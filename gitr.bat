@echo off
cls
setlocal enabledelayedexpansion

:: collect filepaths in variable
set ARGS=%*

:: echo %ARGS% > \\ces-10-cdc\\OSM_CDC_MMRG_work\\users\\bitbucket\\gdm_workflow\\argtest.txt

"C:\Program Files\Git\bin\sh" "//ces-10-cdc//OSM_CDC_MMRG_work//users//bitbucket//gdm_workflow//gitr.sh" %ARGS%
