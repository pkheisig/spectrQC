import React, { useRef, useEffect, useState, useCallback } from 'react';
import { RefreshCw, Play, Save, Check } from 'lucide-react';

interface MatrixAdjusterProps {
  unmixedData: any[];
  matrix: any[];
  setMatrix: (m: any[]) => void;
  runUnmix: (m: any[]) => void;
  markers: string[];
}

const MatrixAdjuster: React.FC<MatrixAdjusterProps> = ({ unmixedData, matrix, setMatrix, runUnmix, markers }) => {
  const [xMarker, setXMarker] = useState<string>('');
  const [yMarker, setYMarker] = useState<string>('');
  const [isAligning, setIsAligning] = useState(false);
  const [dragStart, setDragStart] = useState<{ x: number, y: number } | null>(null);
  const [dragCurrent, setDragCurrent] = useState<{ x: number, y: number } | null>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);

  // Initialize defaults
  useEffect(() => {
    if (markers.length >= 2 && !xMarker) {
      setXMarker(markers[0]);
      setYMarker(markers[1]);
    }
  }, [markers]);

  const draw = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas || !xMarker || !yMarker || unmixedData.length === 0) return;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const w = canvas.width;
    const h = canvas.height;
    ctx.clearRect(0, 0, w, h);

    // Get Data
    const xVals = unmixedData.map(d => d[xMarker]);
    const yVals = unmixedData.map(d => d[yMarker]);
    
    // Auto-scale (simple min/max with padding)
    // For unmixed data, we often care about scaling around 0.
    // Let's find 1st and 99th percentiles to avoid outliers
    const sortedX = [...xVals].sort((a, b) => a - b);
    const sortedY = [...yVals].sort((a, b) => a - b);
    const xMin = sortedX[Math.floor(sortedX.length * 0.01)] || -100;
    const xMax = sortedX[Math.floor(sortedX.length * 0.99)] || 1000;
    const yMin = sortedY[Math.floor(sortedY.length * 0.01)] || -100;
    const yMax = sortedY[Math.floor(sortedY.length * 0.99)] || 1000;
    
    // Expand a bit
    const xRange = xMax - xMin;
    const yRange = yMax - yMin;
    const xPad = xRange * 0.1;
    const yPad = yRange * 0.1;
    const finalXMin = xMin - xPad;
    const finalXMax = xMax + xPad;
    const finalYMin = yMin - yPad;
    const finalYMax = yMax + yPad;

    const toPx = (v: number, min: number, max: number, size: number) => {
      return ((v - min) / (max - min)) * size;
    };

    // Draw axes
    const zeroX = toPx(0, finalXMin, finalXMax, w);
    const zeroY = h - toPx(0, finalYMin, finalYMax, h);
    
    ctx.strokeStyle = '#334155'; // Slate-700
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(zeroX, 0); ctx.lineTo(zeroX, h); // Y-axis
    ctx.moveTo(0, zeroY); ctx.lineTo(w, zeroY); // X-axis
    ctx.stroke();

    // Draw Points
    ctx.fillStyle = '#3b82f6'; // Blue-500
    ctx.globalAlpha = 0.6;
    
    for (let i = 0; i < unmixedData.length; i++) {
        const x = xVals[i];
        const y = yVals[i];
        const px = toPx(x, finalXMin, finalXMax, w);
        const py = h - toPx(y, finalYMin, finalYMax, h);
        ctx.beginPath();
        ctx.arc(px, py, 1.5, 0, Math.PI * 2);
        ctx.fill();
    }
    
    // Draw Drag Line/Interaction
    if (isAligning && dragStart && dragCurrent) {
        ctx.strokeStyle = '#ef4444'; // Red-500
        ctx.lineWidth = 2;
        ctx.setLineDash([5, 5]);
        ctx.beginPath();
        ctx.moveTo(dragStart.x, dragStart.y);
        ctx.lineTo(dragCurrent.x, dragCurrent.y);
        ctx.stroke();
        ctx.setLineDash([]);
        
        // Draw delta text
        const dx = dragCurrent.x - dragStart.x;
        const dy = dragCurrent.y - dragStart.y;
        ctx.fillStyle = '#fff';
        ctx.font = '12px sans-serif';
        ctx.fillText(`Aligning...`, dragCurrent.x + 10, dragCurrent.y);
    }

  }, [unmixedData, xMarker, yMarker, isAligning, dragStart, dragCurrent]);

  // Animation Loop
  useEffect(() => {
    let animationId: number;
    const render = () => {
        draw();
        animationId = requestAnimationFrame(render);
    };
    render();
    return () => cancelAnimationFrame(animationId);
  }, [draw]);


  const handleMouseDown = (e: React.MouseEvent) => {
    if (!isAligning) return;
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;
    setDragStart({ x: e.clientX - rect.left, y: e.clientY - rect.top });
  };

  const handleMouseMove = (e: React.MouseEvent) => {
    if (!isAligning || !dragStart) return;
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;
    setDragCurrent({ x: e.clientX - rect.left, y: e.clientY - rect.top });
  };

  const handleMouseUp = (e: React.MouseEvent) => {
    if (!isAligning || !dragStart) return;
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;
    
    // CALCULATE ADJUSTMENT
    // Visual drag: from Y_start to Y_end at X_pos.
    // We assume the user grabbed the centroid of the population and dragged it to Y=0 (or wherever).
    // The "slope" correction:
    // We want to remove the dependency of Y on X.
    // If we dragged vertically by dy pixels corresponding to dY intensity units,
    // and this happened at X pixels corresponding to dX intensity units.
    
    // 1. Convert pixels back to data units
    // (Need to access scalers from draw function - refactor needed or duplicate logic)
    // Quick Hack: Re-calculate scale here
    const w = rect.width;
    const h = rect.height;
    const xVals = unmixedData.map(d => d[xMarker]);
    const yVals = unmixedData.map(d => d[yMarker]);
    const sortedX = [...xVals].sort((a, b) => a - b);
    const sortedY = [...yVals].sort((a, b) => a - b);
    const xMin = (sortedX[Math.floor(sortedX.length * 0.01)] || -100) - (sortedX[sortedX.length-1]-sortedX[0])*0.1;
    const xMax = (sortedX[Math.floor(sortedX.length * 0.99)] || 1000) + (sortedX[sortedX.length-1]-sortedX[0])*0.1;
    const yMin = (sortedY[Math.floor(sortedY.length * 0.01)] || -100) - (sortedY[sortedY.length-1]-sortedY[0])*0.1;
    const yMax = (sortedY[Math.floor(sortedY.length * 0.99)] || 1000) + (sortedY[sortedY.length-1]-sortedY[0])*0.1;

    const fromPxX = (px: number) => (px / w) * (xMax - xMin) + xMin;
    const fromPxY = (py: number) => ((h - py) / h) * (yMax - yMin) + yMin;

    const startX = fromPxX(dragStart.x);
    const startY = fromPxY(dragStart.y);
    const endY = fromPxY(e.clientY - rect.top);
    
    const dY = endY - startY; // Change in Y wanted
    // Compensation Coefficient Change: alpha = dY / X
    // We assume user clicked ON the population centroid.
    const alpha = dY / (startX || 1); // Avoid div by zero

    // 2. Update Matrix
    // W_new = C * W_old
    // C is identity with C[y, x] = alpha (actually -alpha or +alpha depending on definition)
    // If Y_new = Y_old + alpha * X_old
    // Then we add alpha to the coefficient.
    
    console.log(`Adjusting ${yMarker} by ${alpha} * ${xMarker}`);

    // Apply to matrix
    // We need to find the ROWS corresponding to our markers in W (if W is Marker x Detector)
    // Actually, "Unmixing Matrix" W usually maps Detectors -> Markers.
    // But in our API, we treat M (Markers x Detectors) and W (also Markers x Detectors for saving).
    // Let's assume 'matrix' state is [ {Marker: "A", FL1: ...}, ... ]
    
    // If we are modifying W directly (rows are Markers, cols are Detectors? No, usually cols are Detectors)
    // Post-multiplication compensation:
    // W_new_row_Y = W_old_row_Y + alpha * W_old_row_X
    
    const newMatrix = matrix.map(row => {
        if (row.Marker === yMarker) {
            // This is the row for Y. We add alpha * row_X to it.
            const rowX = matrix.find(r => r.Marker === xMarker);
            if (!rowX) return row;
            
            const newRow = { ...row };
            Object.keys(row).forEach(key => {
                if (key !== 'Marker' && typeof row[key] === 'number') {
                    newRow[key] = Number(row[key]) + alpha * Number(rowX[key]);
                }
            });
            return newRow;
        }
        return row;
    });

    setMatrix(newMatrix);
    runUnmix(newMatrix);
    
    // Reset
    setIsAligning(false);
    setDragStart(null);
    setDragCurrent(null);
  };

  return (
    <div className="flex flex-col h-full">
        {/* Controls */}
        <div className="flex items-center gap-4 mb-2 p-2 bg-slate-800 rounded border border-slate-700">
            <div className="flex items-center gap-2">
                <span className="text-xs text-slate-400">X Axis:</span>
                <select value={xMarker} onChange={e => setXMarker(e.target.value)} className="bg-slate-900 text-xs px-2 py-1 rounded border border-slate-600">
                    {markers.map(m => <option key={m} value={m}>{m}</option>)}
                </select>
            </div>
            <div className="flex items-center gap-2">
                <span className="text-xs text-slate-400">Y Axis:</span>
                <select value={yMarker} onChange={e => setYMarker(e.target.value)} className="bg-slate-900 text-xs px-2 py-1 rounded border border-slate-600">
                    {markers.map(m => <option key={m} value={m}>{m}</option>)}
                </select>
            </div>
            
            <button 
                onClick={() => setIsAligning(!isAligning)} 
                className={`flex items-center gap-2 px-3 py-1 rounded text-xs font-bold transition-colors ${isAligning ? 'bg-red-600 text-white' : 'bg-blue-600 text-white hover:bg-blue-500'}`}
            >
                {isAligning ? 'Cancel' : 'Align Population'}
            </button>
            <div className="text-[10px] text-slate-500 ml-auto">
                {isAligning ? 'Click population center & drag to desired Y' : 'Select axes to visualize'}
            </div>
        </div>

        {/* Canvas */}
        <div className="flex-1 bg-black/20 rounded-lg border border-slate-700 relative overflow-hidden">
            <canvas 
                ref={canvasRef}
                width={600}
                height={400}
                className={`w-full h-full ${isAligning ? 'cursor-crosshair' : 'cursor-default'}`}
                onMouseDown={handleMouseDown}
                onMouseMove={handleMouseMove}
                onMouseUp={handleMouseUp}
                onMouseLeave={() => { if(isAligning) { setDragStart(null); setIsAligning(false); } }}
            />
        </div>
    </div>
  );
};

export default MatrixAdjuster;
