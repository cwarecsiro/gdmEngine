@echo off
cls
setlocal enabledelayedexpansion

:: collect filepaths in variable
set ARGS=%*

:: used in trying to work out how to pass a white-space separated string as an arg across these envs
::echo %ARGS% > \\ces-10-cdc\\OSM_CDC_MMRG_work\\users\\bitbucket\\gdm_workflow\\argtest.txt

"C:\Program Files\Git\bin\sh" %ARGS%
