@echo off
SETLOCAL EnableDelayedExpansion
TITLE HPLC Chromatogram Viewer

REM Navigate to the project root directory containing this script
CD /D "%~dp0"

ECHO ==============================================================
ECHO   HPLC Chromatogram Viewer - Local Scientific Desktop App
ECHO ==============================================================
ECHO.

REM Activate local virtual environment if one exists
IF EXIST ".venv\Scripts\activate.bat" (
    ECHO Activating virtual environment (.venv)...
    CALL ".venv\Scripts\activate.bat"
) ELSE IF EXIST "venv\Scripts\activate.bat" (
    ECHO Activating virtual environment (venv)...
    CALL "venv\Scripts\activate.bat"
)

ECHO Starting local Streamlit server on http://localhost:8501 ...
ECHO (All data remains strictly local; telemetry is disabled.)
ECHO.

python -m streamlit run app.py --server.address=localhost --browser.gatherUsageStats=false

IF %ERRORLEVEL% NEQ 0 (
    ECHO.
    ECHO [ERROR] Failed to launch HPLC Viewer. Ensure Python and requirements.txt are installed:
    ECHO   pip install -r requirements.txt
    PAUSE
)
ENDLOCAL
