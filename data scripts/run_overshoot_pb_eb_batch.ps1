# Runs overshoot_pb_eb_script.m headless on every trial folder of the flies
# collected on a given date under the overshoot data root, one fresh MATLAB
# process per trial (keeps memory bounded; one failure doesn't stop the rest).
# Trials with no tif are skipped; trials that already have bump_results_EB.mat
# are skipped unless -Force. Fly -99 folders (test recordings) are ignored.
#   powershell -File run_overshoot_pb_eb_batch.ps1 [-Date 20261012] [-Force]
param(
    [string]$Date = (Get-Date -Format 'yyyyMMdd'),
    [string]$Root = 'Z:\noah_np123\Data\flyg\overshoot',
    [switch]$Force,
    [string]$Script = 'overshoot_pb_eb_script.m',      # single-structure PB volumes (13 frames/vol): overshoot_bar_stim_script.m
    [string]$Var = 'obs2_trialDir',                    # ... with obs_trialDir
    [string]$DoneFile = 'bump_results_EB.mat'          # ... bump_results.mat
)
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$logDir = Join-Path $scriptDir "overshoot_batch_logs_$Date"
New-Item -ItemType Directory -Force $logDir | Out-Null
$trials = Get-ChildItem $Root -Directory -Filter "$Date-*" | Where-Object { $_.Name -notmatch "^$Date-99_" } |
    ForEach-Object { $fly = $_; Get-ChildItem $fly.FullName -Directory -Filter 'trial_*' | ForEach-Object { $_ } }
foreach ($t in $trials) {
    $tag = "$($t.Parent.Name)__$($t.Name)"
    if (-not (Get-ChildItem $t.FullName -Filter *.tif)) { "SKIP (no tif): $tag" ; continue }
    if (-not $Force -and (Test-Path (Join-Path $t.FullName $DoneFile))) { "SKIP (done): $tag"; continue }
    "RUN: $tag  $(Get-Date -Format HH:mm:ss)"
    $dir = ($t.FullName -replace "'", "''") + '\'
    $cmd = "$Var='$dir'; try, run('$($scriptDir -replace "'", "''")\$Script'); catch ME, disp(getReport(ME)); exit(1); end"
    & matlab -batch $cmd *> (Join-Path $logDir "$tag.log")
    "  exit $LASTEXITCODE  $(Get-Date -Format HH:mm:ss)"
}
"ALL DONE"
