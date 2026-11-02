@echo off
title CaM-RyR2 Protein Redesign Portal Launcher
echo Activating pymol-env and launching portal...
echo.

:: 0. Run from this folder so uploaded_pdbs\ etc. resolve correctly
cd /d "%~dp0"

:: 1. Initialize Conda environment hooks
call "C:\Users\charh\miniconda3\Scripts\activate.bat" "C:\Users\charh\miniconda3"

:: 2. Activate pymol-env
call conda activate pymol-env

:: 3. Launch Streamlit using the active environment's executable
python -m streamlit run gui.py

pause