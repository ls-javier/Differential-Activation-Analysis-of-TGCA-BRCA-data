@echo off
echo Starting TCGA-BRCA Differential Activation Analysis Shiny App...
echo.
echo The app will open in your default browser.
echo Press Ctrl+C to stop the server.
echo.
cd /d "%~dp0"
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" -e "shiny::runApp('app', port = 7890, launch.browser = TRUE)"
pause
