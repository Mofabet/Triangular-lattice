@echo off
cd /d "%~dp0"
echo === trilattice : full interactive dashboard ===
python -m pip install -e ".[fast]" -q
if errorlevel 1 (
    echo.
    echo Install failed -- see error above.
    pause
    exit /b 1
)
python -m trilattice.cli animate examples\config.toml -T 1800
if errorlevel 1 pause
