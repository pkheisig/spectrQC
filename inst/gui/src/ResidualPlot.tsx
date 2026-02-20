import React, { useRef, useEffect, useState, useCallback } from 'react';

interface ResidualPlotProps {
  xKey: string;
  yKey: string;
  data: any[];
  size: number;
  pointColor: string;
  pointOpacity: number;
  pointSize: number;
  sensitivity: number;
  onAdjust: (xKey: string, yKey: string, alpha: number) => void;
}

const ResidualPlot: React.FC<ResidualPlotProps> = ({
  xKey,
  yKey,
  data,
  size,
  pointColor,
  pointOpacity,
  pointSize,
  sensitivity,
  onAdjust
}) => {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [isDragging, setIsDragging] = useState(false);
  const [dragStart, setDragStart] = useState<{ x: number, y: number, dataX: number, dataY: number } | null>(null);
  const [dragCurrent, setDragCurrent] = useState<{ x: number, y: number } | null>(null);

  // Helper to get ranges (memoize in real app, but fast enough for 5k points here)
  const getRanges = () => {
    if (data.length === 0) return { xMin: 0, xMax: 1, yMin: 0, yMax: 1 };

    // We want a stable range, but for residuals it might change.
    // Ideally, for residuals (compensated data), we expect clusters around 0.
    // SpectroFlo usually auto-scales.
    // Let's use robust min/max.
    let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;

    // Use a subset for speed if needed, but 5k is fine
    for (let i = 0; i < data.length; i++) {
      const x = Number(data[i][xKey]);
      const y = Number(data[i][yKey]);
      if (isNaN(x) || isNaN(y)) continue;
      if (x < xMin) xMin = x;
      if (x > xMax) xMax = x;
      if (y < yMin) yMin = y;
      if (y > yMax) yMax = y;
    }

    // Add padding
    const xPad = (xMax - xMin) * 0.1 || 100;
    const yPad = (yMax - yMin) * 0.1 || 100;

    return {
      xMin: xMin - xPad,
      xMax: xMax + xPad,
      yMin: yMin - yPad,
      yMax: yMax + yPad
    };
  };

  const draw = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas || data.length === 0) return;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const w = canvas.width;
    const h = canvas.height;
    ctx.clearRect(0, 0, w, h);

    const { xMin, xMax, yMin, yMax } = getRanges();
    const xRange = xMax - xMin || 1;
    const yRange = yMax - yMin || 1;

    // Draw Axes (Zero lines)
    if (xMin < 0 && xMax > 0) {
      const x0 = ((0 - xMin) / xRange) * w;
      ctx.strokeStyle = '#e2e8f0'; // slate-200
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(x0, 0); ctx.lineTo(x0, h);
      ctx.stroke();
    }
    if (yMin < 0 && yMax > 0) {
      const y0 = h - ((0 - yMin) / yRange) * h;
      ctx.strokeStyle = '#e2e8f0';
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(0, y0); ctx.lineTo(w, y0);
      ctx.stroke();
    }

    // Draw Points
    ctx.fillStyle = pointColor;
    ctx.globalAlpha = pointOpacity;

    for (let i = 0; i < data.length; i++) {
      const x = Number(data[i][xKey]);
      const y = Number(data[i][yKey]);
      if (isNaN(x) || isNaN(y)) continue;

      const px = ((x - xMin) / xRange) * w;
      const py = h - ((y - yMin) / yRange) * h;

      ctx.beginPath();
      ctx.arc(px, py, pointSize, 0, Math.PI * 2);
      ctx.fill();
    }

    // Draw Drag Overlay
    if (isDragging && dragStart && dragCurrent) {
      ctx.strokeStyle = '#ef4444'; // red-500
      ctx.lineWidth = 2;
      ctx.setLineDash([4, 2]);
      ctx.beginPath();
      ctx.moveTo(dragStart.x, dragStart.y);
      ctx.lineTo(dragCurrent.x, dragCurrent.y);
      ctx.stroke();
      ctx.setLineDash([]);

      // Calculate and draw delta
      // const dx = dragCurrent.x - dragStart.x;
      // const dy = dragCurrent.y - dragStart.y;

      // Map back to data space to show alpha?
      // Maybe simpler: just show visual line
    }

  }, [data, xKey, yKey, pointColor, pointOpacity, pointSize, isDragging, dragStart, dragCurrent]);

  useEffect(() => {
    draw();
  }, [draw]);

  const handleMouseDown = (e: React.MouseEvent) => {
    // Only allow dragging if X and Y are different
    if (xKey === yKey) return;

    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;
    const ex = e.clientX - rect.left;
    const ey = e.clientY - rect.top;

    // Map click to data space
    const { xMin, xMax, yMin, yMax } = getRanges();
    const w = rect.width;
    const h = rect.height;

    const dataX = (ex / w) * (xMax - xMin) + xMin;
    const dataY = ((h - ey) / h) * (yMax - yMin) + yMin;

    setIsDragging(true);
    setDragStart({ x: ex, y: ey, dataX, dataY });
    setDragCurrent({ x: ex, y: ey });
  };

  const handleMouseMove = (e: React.MouseEvent) => {
    if (!isDragging) return;
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;
    setDragCurrent({ x: e.clientX - rect.left, y: e.clientY - rect.top });
    draw(); // Re-draw to show line
  };

  const handleMouseUp = (e: React.MouseEvent) => {
    if (!isDragging || !dragStart) return;
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;

    const endX = e.clientX - rect.left;
    const endY = e.clientY - rect.top;

    // Calculate pixel deltas
    const pixelDeltaX = endX - dragStart.x;
    const pixelDeltaY = endY - dragStart.y;

    // Map to data space
    const { xMin, xMax, yMin, yMax } = getRanges();
    const w = rect.width;
    const h = rect.height;

    const dataDeltaX = (pixelDeltaX / w) * (xMax - xMin);
    const dataDeltaY = -(pixelDeltaY / h) * (yMax - yMin); // Negative because canvas Y is inverted

    // Determine drag direction:
    // horizontal drag should feel like Y-axis correction,
    // vertical drag should feel like X-axis correction.

    const absX = Math.abs(pixelDeltaX);
    const absY = Math.abs(pixelDeltaY);

    // Minimum drag threshold
    if (absX < 5 && absY < 5) {
      setIsDragging(false);
      setDragStart(null);
      setDragCurrent(null);
      return;
    }

    // Apply user sensitivity multiplier.
    const baseAlpha = 0.00005 * sensitivity;

    if (absX > absY) {
      // Horizontal drag: adjust Y based on X (intuitive left-right correction).
      const alpha = -dataDeltaX * baseAlpha;
      onAdjust(xKey, yKey, alpha);
    } else {
      // Vertical drag: adjust X based on Y (intuitive up-down correction).
      const alpha = -dataDeltaY * baseAlpha;
      onAdjust(yKey, xKey, alpha);
    }

    setIsDragging(false);
    setDragStart(null);
    setDragCurrent(null);
    draw();
  };

  return (
    <canvas
      ref={canvasRef}
      width={size}
      height={size}
      className={`w-full h-full touch-none ${xKey !== yKey ? 'cursor-crosshair hover:ring-1 hover:ring-blue-400' : 'cursor-default'}`}
      onMouseDown={handleMouseDown}
      onMouseMove={handleMouseMove}
      onMouseUp={handleMouseUp}
      onMouseLeave={() => { if (isDragging) { setIsDragging(false); setDragStart(null); draw(); } }}
    />
  );
};

export default ResidualPlot;
