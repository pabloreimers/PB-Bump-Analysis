# Step 5: current Z: status + genotype per trial -> ../tables/p2x2_z_trial_status.csv
# (read by fig_ionto_celltypes.m). Run from the work dir holding z_trials.json:
#   python -I <path>/5_z_status.py
import csv,json,os
B=chr(92)
rows=json.load(open('z_trials.json'))
out=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','tables','p2x2_z_trial_status.csv')
with open(out,'w',newline='') as f:
    w=csv.writer(f); w.writerow(['key','z_status','z_geno','z_path'])
    for r in rows:
        g=r['geno'] or ('lpsp' if 'blind' in r['tname'].lower() else '')
        st=r['status']
        if st=='ok' and r['date'] and r['date'].startswith('2024') and '_low' in r['tname'].lower(): st='2024 low ATP'
        w.writerow([r['key'],st,g,'Z:'+B+'pablo'+B+r['path'].replace('/',B)])
print('wrote',len(rows))
