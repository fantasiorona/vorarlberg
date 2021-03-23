@echo off
start /B ./build/Debug/SimpleGA.exe > ./simple_ga_test.txt 2>&1
if errorlevel 1 (
   echo Failure Reason Given is %errorlevel%
   exit /b %errorlevel%
)