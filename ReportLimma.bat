@echo off

rem Initialize the variable for the argument
set "config="

rem Check if an argument was provided
if not "%~1"=="" (
    set "config=%~1"
    echo Configuration file path: %~1
) else (
    echo Configuration file was not provided.
    exit /b 1
)

for %%I in ("%~dp0") do set "base_path=%%~fI"

set "Rexe_relative=R-Portable\App\R-Portable\bin\x64\Rscript.exe"
set "Rexe_path=%base_path%%Rexe_relative%"
set "Rscript=%base_path%app_wo_GUI.R"

echo Executing command: "%Rexe_path%" "--vanilla %Rscript% %config%"

"%Rexe_path%" "--vanilla %Rscript% %config%"