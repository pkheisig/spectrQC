#!/bin/bash

# Kill any existing processes on these ports
lsof -ti:8000 | xargs kill -9 2>/dev/null
lsof -ti:5173 | xargs kill -9 2>/dev/null

echo "Starting spectrQC Backend (R Plumber)..."
# Log R output to a file so we can debug
Rscript run_gui.R > backend.log 2>&1 &

echo "Starting spectrQC Frontend (Vite)..."
cd gui
npm run dev