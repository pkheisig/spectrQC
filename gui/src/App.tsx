import { useState, useEffect, useRef } from 'react';
import { Settings2, Activity, Save, RefreshCw, FileText, Check, Sliders, Palette, Circle, Sun, Moon, Maximize2 } from 'lucide-react';
import axios from 'axios';
import ResidualPlot from './ResidualPlot';

const API_BASE = 'http://localhost:8001';

interface MatrixRow {
  Marker: string;
  [key: string]: any;
}

const COLOR_PALETTES = {
  default: ['#60a5fa', '#f472b6', '#34d399', '#fbbf24', '#a78bfa', '#fb7185', '#38bdf8', '#4ade80', '#facc15', '#c084fc'],
  neon: ['#00ff87', '#ff00ff', '#00ffff', '#ff6b6b', '#ffd93d', '#6bcb77', '#4d96ff', '#ff9f45', '#ee6ff8', '#45b7d1'],
  pastel: ['#a8d8ea', '#aa96da', '#fcbad3', '#ffffd2', '#b5eaea', '#ffc8dd', '#bde0fe', '#cdb4db', '#fde4cf', '#caffbf'],
  warm: ['#ff6b35', '#f7c59f', '#efa00b', '#d65108', '#ff8243', '#ffb627', '#ff7b54', '#ffa07a', '#f4845f', '#e76f51'],
  cool: ['#5390d9', '#7400b8', '#6930c3', '#5e60ce', '#4ea8de', '#48bfe3', '#56cfe1', '#64dfdf', '#72efdd', '#80ffdb']
};

const SCATTER_COLORS = ['#3b82f6', '#ef4444', '#22c55e', '#f59e0b', '#8b5cf6', '#ec4899', '#14b8a6', '#f97316'];

const THEMES = {
  dark: {
    bg: 'bg-slate-950',
    navBg: 'bg-slate-900/40',
    sideBg: 'bg-slate-900/20',
    cardBg: 'bg-slate-900/30',
    border: 'border-white/5',
    text: 'text-slate-200',
    textMuted: 'text-slate-500',
    textDim: 'text-slate-400',
    inputBg: 'bg-slate-800',
    cellBg: 'bg-slate-900/50',
    gridLine: '#1e293b',
    labelBg: 'bg-slate-950/50',
  },
  light: {
    bg: 'bg-gray-50',
    navBg: 'bg-white/80',
    sideBg: 'bg-white/60',
    cardBg: 'bg-white',
    border: 'border-gray-200',
    text: 'text-gray-900',
    textMuted: 'text-gray-500',
    textDim: 'text-gray-600',
    inputBg: 'bg-gray-100',
    cellBg: 'bg-gray-50',
    gridLine: '#e5e7eb',
    labelBg: 'bg-gray-100',
  }
};

const App = () => {
  const [matrices, setMatrices] = useState<string[]>([]);
  const [currentFile, setCurrentFile] = useState('refined_reference_matrix.csv');
  const [matrix, setMatrix] = useState<MatrixRow[]>([]);
  const [loading, setLoading] = useState(true);
  const [detectors, setDetectors] = useState<string[]>([]);
  const [detectorLabels, setDetectorLabels] = useState<string[]>([]);
  const [selectedMarkers, setSelectedMarkers] = useState<string[]>([]);
  const [unmixedData, setUnmixedData] = useState<any[]>([]);
  const [rawData, setRawData] = useState<any[]>([]);
  const [isDragging, setIsDragging] = useState(false);
  const [isUnmixingMatrix, setIsUnmixingMatrix] = useState(false);

  const [lineWidth, setLineWidth] = useState(0.8);
  const [lineOpacity, setLineOpacity] = useState(0.85);
  const [colorPalette, setColorPalette] = useState<keyof typeof COLOR_PALETTES>('default');
  const [showControls, setShowControls] = useState(true);

  const [signatureHeight, setSignatureHeight] = useState(300);
  const [signatureDetWidth, setSignatureDetWidth] = useState(26);
  const [residualCellSize, setResidualCellSize] = useState(130);
  const [pointSize, setPointSize] = useState(1.5);
  const [pointOpacity, setPointOpacity] = useState(0.5);
  const [pointColor, setPointColor] = useState('#3b82f6');
  const [dragSensitivity, setDragSensitivity] = useState(0.5);

  const [theme, setTheme] = useState<'dark' | 'light'>('light');
  const [pageScroll, setPageScroll] = useState(true);
  const t = THEMES[theme];

  useEffect(() => {
    fetchMatrices();
    fetchData();
  }, []);

  const fetchMatrices = async () => {
    const res = await axios.get(`${API_BASE}/matrices`);
    setMatrices(Array.isArray(res.data) ? res.data : []);
  };

  const fetchData = async (filename = currentFile) => {
    setLoading(true);
    const resMatrix = await axios.get(`${API_BASE}/load_matrix?filename=${filename}`);
    if (!resMatrix.data.error) {
      const matrixData = resMatrix.data;
      setMatrix(matrixData);

      const isUnmix = filename.toLowerCase().includes('unmixing') || filename.toLowerCase().includes('w_');
      setIsUnmixingMatrix(isUnmix);

      const detNames = Object.keys(matrixData[0]).filter(k => k !== 'Marker');
      setDetectors(detNames);
      const allMarkers = matrixData.map((r: any) => r.Marker);
      setSelectedMarkers(allMarkers);
      const resData = await axios.get(`${API_BASE}/data`);
      if (!resData.data.error) {
        const raw = resData.data.raw_data;
        setRawData(raw);
        setDetectorLabels(resData.data.detector_labels || detNames);
        await runUnmix(matrixData, raw, filename);
      }
    }
    setCurrentFile(filename);
    setLoading(false);
  };

  const runUnmix = async (currentM: any[], currentRaw: any[], filename = currentFile) => {
    const M_obj: any = {};
    currentM.forEach(row => {
      M_obj[row.Marker] = { ...row };
      delete M_obj[row.Marker].Marker;
    });

    // Check type based on filename if string, or boolean if provided as override
    // Here we use state or passed filename
    let useType = 'reference';
    if (typeof filename === 'string') {
      if (filename.toLowerCase().includes('unmixing') || filename.toLowerCase().includes('w_')) useType = 'unmixing';
    } else if (typeof filename === 'boolean') {
      useType = filename ? 'unmixing' : 'reference';
    }

    const res = await axios.post(`${API_BASE}/unmix`, {
      matrix_json: M_obj,
      raw_data_json: currentRaw,
      type: useType
    });
    setUnmixedData(res.data);
  };

  const handleResidualAdjust = (xMarker: string, yMarker: string, alpha: number) => {
    const newMatrix = matrix.map(row => {
      if (row.Marker === yMarker) {
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
    runUnmix(newMatrix, rawData, isUnmixingMatrix as any);
  };

  const saveMatrix = async () => {
    const newName = currentFile.replace('.csv', '_adjusted.csv');
    await axios.post(`${API_BASE}/save_matrix`, {
      filename: newName,
      matrix_json: matrix
    });
    alert(`Saved to ${newName}`);
    fetchMatrices();
  };

  const svgRef = useRef<SVGSVGElement>(null);
  const handleMouseMove = (e: React.MouseEvent) => {
    if (!isDragging || selectedMarkers.length !== 1) return;
    const marker = selectedMarkers[0];
    const svg = svgRef.current;
    if (!svg) return;
    const rect = svg.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;
    const detIdx = Math.round((x / rect.width) * (detectors.length - 1));
    const detName = detectors[detIdx];
    const newVal = Math.max(0, Math.min(1, 1 - (y / rect.height)));
    const nextMatrix = matrix.map(row => {
      if (row.Marker === marker) return { ...row, [detName]: newVal };
      return row;
    });
    setMatrix(nextMatrix);
    runUnmix(nextMatrix, rawData, isUnmixingMatrix as any);
  };

  const colors = COLOR_PALETTES[colorPalette];
  const markerNames = matrix.map(m => m.Marker);
  const chartWidth = detectors.length * signatureDetWidth;

  if (loading) return (
    <div className={`h-screen w-screen flex items-center justify-center ${t.bg} text-blue-400 font-medium`}>
      <RefreshCw className="animate-spin mr-3" size={24} /> Initializing Workspace...
    </div>
  );

  const lowerTriangleCells = selectedMarkers.length * (selectedMarkers.length - 1) / 2;

  return (
    <div className={`min-h-screen ${pageScroll ? '' : 'h-screen'} w-screen flex flex-col ${t.bg} ${t.text} font-sans ${pageScroll ? 'overflow-auto' : 'overflow-hidden'}`}>
      {/* Navbar */}
      <header className={`h-14 shrink-0 border-b ${t.border} flex items-center px-4 justify-between ${t.navBg} backdrop-blur-xl sticky top-0 z-50`}>
        <div className="flex items-center gap-3">
          <div className="w-10 h-10 bg-blue-600 rounded-lg flex items-center justify-center">
            <Activity className="text-white" size={20} />
          </div>
          <div>
            <h1 className="font-bold text-base leading-tight">spectrQC</h1>
            <p className={`text-[10px] ${t.textMuted} uppercase tracking-widest`}>Interactive Spectral Tuner</p>
          </div>
        </div>
        <div className="flex items-center gap-3">
          <button onClick={() => setPageScroll(!pageScroll)} className={`p-2 rounded-lg ${pageScroll ? 'bg-blue-600/20 text-blue-400' : `hover:bg-white/5 ${t.textDim}`}`} title="Page scroll mode">
            <Maximize2 size={18} />
          </button>
          <button onClick={() => setTheme(theme === 'dark' ? 'light' : 'dark')} className={`p-2 rounded-lg hover:bg-white/10 ${t.textDim}`}>
            {theme === 'dark' ? <Sun size={18} /> : <Moon size={18} />}
          </button>
          <button onClick={() => setShowControls(!showControls)} className={`p-2 rounded-lg ${showControls ? 'bg-blue-600/20 text-blue-400' : `hover:bg-white/5 ${t.textDim}`}`}>
            <Sliders size={18} />
          </button>
          <button onClick={() => fetchData()} className={`p-2 hover:bg-white/5 rounded-lg ${t.textDim}`}>
            <RefreshCw size={18} />
          </button>
          <button onClick={saveMatrix} className="bg-blue-600 hover:bg-blue-500 text-white px-4 py-1.5 rounded-lg text-sm font-semibold flex items-center gap-2">
            <Save size={14} /> Save
          </button>
        </div>
      </header>

      <div className={`flex-1 flex ${pageScroll ? '' : 'overflow-hidden'}`}>
        {/* Left Sidebar */}
        <aside className={`w-52 shrink-0 border-r ${t.border} ${t.sideBg} flex flex-col ${pageScroll ? '' : 'overflow-auto'} sticky top-14 h-[calc(100vh-56px)]`}>
          <div className="p-3">
            <h3 className={`text-[10px] font-bold ${t.textMuted} uppercase tracking-widest mb-2`}>Matrices</h3>
            <div className="space-y-1">
              {matrices.map(m => (
                <button key={m} onClick={() => fetchData(m)} className={`w-full flex items-center gap-2 px-3 py-1.5 rounded text-xs ${currentFile === m ? 'bg-blue-600/10 text-blue-400' : `hover:bg-white/5 ${t.textDim}`}`}>
                  <FileText size={12} />
                  <span className="truncate flex-1 text-left">{m}</span>
                </button>
              ))}
            </div>
          </div>
          <div className="flex-1 overflow-y-auto px-3 pb-3">
            <h3 className={`text-[10px] font-bold ${t.textMuted} uppercase tracking-widest mb-2`}>Signatures</h3>
            <div className="space-y-0.5">
              {matrix.map((m, idx) => (
                <label key={m.Marker} className={`flex items-center justify-between px-3 py-1 rounded cursor-pointer text-xs ${selectedMarkers.includes(m.Marker) ? `${t.cardBg} ${t.text}` : `hover:bg-white/5 ${t.textMuted}`}`}>
                  <div className="flex items-center gap-2">
                    <input type="checkbox" className="hidden" checked={selectedMarkers.includes(m.Marker)} onChange={() => {
                      setSelectedMarkers(prev => prev.includes(m.Marker) ? prev.filter(x => x !== m.Marker) : [...prev, m.Marker]);
                    }} />
                    <div className="w-2 h-2 rounded-full" style={{ backgroundColor: colors[idx % colors.length], opacity: selectedMarkers.includes(m.Marker) ? 1 : 0.3 }} />
                    <span>{m.Marker}</span>
                  </div>
                  {selectedMarkers.includes(m.Marker) && <Check size={10} className="text-blue-400" />}
                </label>
              ))}
            </div>
          </div>
        </aside>

        {/* Main Content */}
        <div className={`flex-1 flex flex-col p-3 gap-3 ${pageScroll ? '' : 'overflow-auto'}`}>
          {/* Controls */}
          {showControls && (
            <div className={`shrink-0 ${t.cardBg} rounded-lg border ${t.border} p-3 text-xs`}>
              <div className="flex flex-wrap gap-x-6 gap-y-2">
                <div className="flex items-center gap-3">
                  <span className={`font-bold ${t.textMuted} uppercase`}>Sig:</span>
                  <span className={t.textMuted}>W</span><input type="range" min="0.2" max="2" step="0.1" value={lineWidth} onChange={e => setLineWidth(Number(e.target.value))} className="w-12 h-1 rounded cursor-pointer accent-blue-500" /><span className={`${t.textMuted} w-4`}>{lineWidth}</span>
                  <span className={t.textMuted}>Op</span><input type="range" min="0.1" max="1" step="0.05" value={lineOpacity} onChange={e => setLineOpacity(Number(e.target.value))} className="w-12 h-1 rounded cursor-pointer accent-blue-500" /><span className={`${t.textMuted} w-6`}>{Math.round(lineOpacity * 100)}%</span>
                  <Palette size={12} className={t.textMuted} /><select value={colorPalette} onChange={e => setColorPalette(e.target.value as keyof typeof COLOR_PALETTES)} className={`${t.inputBg} ${t.text} px-2 py-1 rounded border ${t.border} outline-none text-xs`}><option value="default">Default</option><option value="neon">Neon</option><option value="pastel">Pastel</option><option value="warm">Warm</option><option value="cool">Cool</option></select>
                  <span className={t.textMuted}>H</span><input type="range" min="200" max="500" step="25" value={signatureHeight} onChange={e => setSignatureHeight(Number(e.target.value))} className="w-12 h-1 rounded cursor-pointer accent-blue-500" /><span className={`${t.textMuted} w-6`}>{signatureHeight}</span>
                  <span className={t.textMuted}>DetW</span><input type="range" min="10" max="30" step="2" value={signatureDetWidth} onChange={e => setSignatureDetWidth(Number(e.target.value))} className="w-10 h-1 rounded cursor-pointer accent-blue-500" /><span className={`${t.textMuted} w-4`}>{signatureDetWidth}</span>
                </div>
                <div className="flex items-center gap-3">
                  <span className={`font-bold ${t.textMuted} uppercase`}>Res:</span>
                  <span className={t.textMuted}>Pt</span><input type="range" min="0.5" max="4" step="0.25" value={pointSize} onChange={e => setPointSize(Number(e.target.value))} className="w-10 h-1 rounded cursor-pointer accent-emerald-500" /><span className={`${t.textMuted} w-4`}>{pointSize}</span>
                  <span className={t.textMuted}>Op</span><input type="range" min="0.1" max="1" step="0.05" value={pointOpacity} onChange={e => setPointOpacity(Number(e.target.value))} className="w-10 h-1 rounded cursor-pointer accent-emerald-500" /><span className={`${t.textMuted} w-6`}>{Math.round(pointOpacity * 100)}%</span>
                  <Circle size={10} className={t.textMuted} /><div className="flex gap-1">{SCATTER_COLORS.map(c => (<button key={c} onClick={() => setPointColor(c)} className={`w-4 h-4 rounded-full border ${pointColor === c ? 'border-white' : 'border-transparent'}`} style={{ backgroundColor: c }} />))}</div>
                  <span className={t.textMuted}>Cell</span><input type="range" min="80" max="200" step="10" value={residualCellSize} onChange={e => setResidualCellSize(Number(e.target.value))} className="w-10 h-1 rounded cursor-pointer accent-emerald-500" /><span className={`${t.textMuted} w-6`}>{residualCellSize}</span>
                </div>
              </div>
            </div>
          )}

          {/* Signature Editor */}
          <div className="shrink-0">
            <div className="flex items-center gap-2 mb-1">
              <div className="p-1.5 bg-indigo-500/10 rounded"><Settings2 size={14} className="text-indigo-400" /></div>
              <h2 className="font-bold text-sm">Signature Editor</h2>
              <span className={`text-xs ${t.textMuted}`}>({selectedMarkers.length} / {detectors.length})</span>
            </div>
            <div className={`${t.cardBg} rounded-lg border ${t.border} overflow-auto`} style={{ maxHeight: signatureHeight + 70 }}>
              <div style={{ width: chartWidth }}>
                <svg
                  ref={svgRef}
                  className="cursor-crosshair block"
                  width={chartWidth}
                  height={signatureHeight}
                  viewBox={`0 0 ${chartWidth} 100`}
                  preserveAspectRatio="none"
                  onMouseDown={() => setIsDragging(true)}
                  onMouseUp={() => setIsDragging(false)}
                  onMouseLeave={() => setIsDragging(false)}
                  onMouseMove={handleMouseMove}
                >
                  {[0, 0.25, 0.5, 0.75, 1].map(tick => (
                    <line key={tick} x1="0" y1={(1 - tick) * 100} x2={chartWidth} y2={(1 - tick) * 100} stroke={t.gridLine} strokeWidth="0.2" strokeDasharray="2 2" />
                  ))}
                  {selectedMarkers.map((m, mIdx) => {
                    const row = matrix.find(r => r.Marker === m);
                    if (!row) return null;
                    const points = detectors.map((det, dIdx) => {
                      const px = dIdx * signatureDetWidth + signatureDetWidth / 2;
                      const py = (1 - Number(row[det])) * 100;
                      return `${px},${py}`;
                    }).join(' ');
                    return (
                      <polyline key={m} points={points} fill="none" stroke={colors[mIdx % colors.length]} strokeWidth={lineWidth} strokeLinejoin="round" strokeLinecap="round" style={{ opacity: lineOpacity }} />
                    );
                  })}
                </svg>
                <div className={`flex border-t ${t.border} ${t.labelBg}`} style={{ width: chartWidth }}>
                  {detectorLabels.map((label, idx) => (
                    <div key={idx} className="relative" style={{ width: signatureDetWidth, height: 60 }}>
                      <span className={`absolute text-[8px] ${t.textMuted} whitespace-nowrap origin-top-left`} style={{ left: signatureDetWidth / 2, top: 4, transform: 'rotate(90deg) translateX(0) translateY(-50%)' }}>
                        {label}
                      </span>
                    </div>
                  ))}
                </div>
              </div>
            </div>
          </div>

          {/* Residual Monitor - Lower Triangle Only */}
          <div className={pageScroll ? '' : 'flex-1 min-h-0'}>
            <div className="flex items-center gap-2 mb-1">
              <div className="p-1.5 bg-emerald-500/10 rounded"><Activity size={14} className="text-emerald-400" /></div>
              <h2 className="font-bold text-sm">Residual Monitor</h2>
              <span className={`text-xs ${t.textMuted}`}>({lowerTriangleCells} cells, {unmixedData.length} events)</span>
              <span className="text-[10px] text-red-400 ml-4 font-bold flex items-center gap-1">
                <Sliders size={10} /> Drag on plots to adjust crosstalk
              </span>
              <span className={`text-xs ${t.textMuted} ml-2`}>Sens:</span>
              <input type="range" min="0.01" max="1" step="0.01" value={dragSensitivity} onChange={e => setDragSensitivity(Number(e.target.value))} className="w-16 h-1 rounded cursor-pointer accent-red-500" />
              <input type="number" min="0.01" max="1" step="0.01" value={dragSensitivity} onChange={e => setDragSensitivity(Math.max(0.01, Math.min(1, Number(e.target.value))))} className={`w-12 text-xs ${t.inputBg} ${t.text} border ${t.border} rounded px-1 py-0.5 ml-1`} />
            </div>
            <div className={`${t.cardBg} rounded-lg border ${t.border} p-2 overflow-auto`} style={{ maxHeight: pageScroll ? 'none' : 'calc(100% - 28px)' }}>
              <div style={{ display: 'inline-block' }}>
                {/* Header row for first x-axis label */}
                <div className="flex items-center" style={{ display: selectedMarkers.includes(markerNames[0]) ? 'flex' : 'none' }}>
                  <div style={{ width: 60 }}></div>
                  <div className={`flex items-center justify-center text-[9px] ${t.textMuted} font-medium`} style={{ width: residualCellSize, height: 20 }}>
                    {markerNames[0]}
                  </div>
                </div>
                {/* Lower triangle rows with stair-style x-axis labels */}
                {markerNames.map((rowName, rowIdx) => {
                  if (rowIdx === 0) return null;
                  const rowSelected = selectedMarkers.includes(rowName);
                  const isLastRow = rowIdx === markerNames.length - 1;
                  return (
                    <div key={`row-${rowIdx}`} className="flex items-center" style={{ display: rowSelected ? 'flex' : 'none' }}>
                      <div className={`text-[9px] ${t.textMuted} text-right font-medium pr-2 truncate`} style={{ width: 60 }}>{rowName}</div>
                      {markerNames.slice(0, rowIdx).map((colName, colIdx) => {
                        const colSelected = selectedMarkers.includes(colName);
                        return (
                          <div key={`cell-${rowIdx}-${colIdx}`} className={`border ${t.border} ${t.cellBg} p-0.5`} style={{ width: residualCellSize, height: residualCellSize, display: colSelected ? 'block' : 'none' }}>
                            <ResidualPlot
                              xKey={colName}
                              yKey={rowName}
                              data={unmixedData}
                              size={residualCellSize - 6}
                              pointColor={pointColor}
                              pointOpacity={pointOpacity}
                              pointSize={pointSize}
                              sensitivity={dragSensitivity}
                              onAdjust={handleResidualAdjust}
                            />
                          </div>
                        );
                      })}
                      {/* Diagonal label cell - x-axis label, but skip for last row to avoid redundancy */}
                      {!isLastRow && (
                        <div className={`flex items-center justify-center text-[9px] ${t.textMuted} font-medium`} style={{ width: residualCellSize, height: residualCellSize }}>
                          {markerNames[rowIdx]}
                        </div>
                      )}
                    </div>
                  );
                })}
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default App;