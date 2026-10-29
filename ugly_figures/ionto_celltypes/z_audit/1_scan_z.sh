#!/bin/sh
# Step 1: list every file that marks a trial folder under the five Z:\pablo\lpsp_p2x2* roots.
# Run from a scratch/work directory (Git Bash); writes z_scan.tsv there (~3 min over the network).
cd /z/pablo && find lpsp_p2x2 lpsp_p2x2_redo lpsp_p2x2_reredo lpsp_p2x2_syt8s lpsp_p2x2_walking \
  \( -iname '*ficTracData_DAQ.mat' -o -iname 'imagingData*.mat' -o -iname 'mask.mat' -o -type d -iname 'registration*' \
     -o -iname 'trialSettings.csv' -o -iname '*.tif' \) -printf '%p\t%s\n' > "$OLDPWD/z_scan.tsv"
