# Step 6 (after fig_ionto_celltypes.m has run): merge the Z: scan with tables/p2x2_trials_used.csv and
# tables/per_fly_ejection_features.csv into tables/p2x2_trial_ledger.csv + p2x2_fly_ledger.csv.
#   python -I 6_final_ledger.py <work dir> <path to ../tables>
import csv, re, collections, json, sys
B = chr(92)
SP = sys.argv[1]; OUT = sys.argv[2]
z = list(csv.DictReader(open(SP + '/p2x2_trial_ledger.csv')))
used = list(csv.DictReader(open(OUT + '/p2x2_trials_used.csv')))
norm = lambda t: re.sub(r'^(\d{8}-\d+)[-_]', lambda m: m.group(1) + '_', t.lower())
U = {u['key']: u for u in used}
# .mat copies per trial (old inventory + the new build file)
inv = list(csv.DictReader(open(SP + '/p2x2_inventory.csv')))
ncopy = collections.Counter(); npath = collections.defaultdict(set)
for r in inv:
    p = [x for x in r['meta'].split(B) if x]
    if p[-1].startswith('registration'): p = p[:-1]
    k = p[2] + '|' + norm(p[-1]); ncopy[k] += 1; npath[k].add(B.join(p))
CALL={}
for x in csv.DictReader(open(OUT+'/per_fly_ejection_features.csv')): CALL[x['fly']]=(x['call'],x['conf'])
rows = []
for r in z:
    k = r['root'] + '|' + norm(r['tname'])
    u = U.get(k)
    copies = ncopy[k] + (1 if (u and 'redo_syt8s' in u['file']) else 0)
    why = (u['why'] if u else '') or ''
    if not u:
        if r['status'] != 'ok': why = r['status']
        elif not r['reg'] or 'ch1' in r['reg'] or r['reg'].startswith('imagingData_2024') or re.match(r'imagingData_trial0', r['reg']):
            why = 'never saved; old single-channel registration'
        else: why = 'never saved'
    m = re.search(r'fly\s*(\d+)', r['flyfolder'] or '', re.I)
    fly = f"{r['root']}|{r['date']}|fly {m.group(1)}" if m else f"{r['root']}|{r['tname']}"
    rows.append(dict(root=r['root'], date=r['date'], fly=fly, flyfolder=r['flyfolder'], trial=r['tname'],
                     genotype=(u['geno'] if u and u['geno'] else r['geno'] or '?'), indicator=r['ind'],
                     z_status=r['status'], mat_copies=copies, stored_paths=len(npath[k]),
                     used=(u['used'] in ('1','true','True') if u else False), why_not_used=why if not (u and u['used'] in ('1','true','True')) else '',
                     stratum=(u['cond'] if u else ''), n_stims=(u['n_stims'] if u else ''), n_ej=(u['n_ej'] if u else ''),
                     z_path=r['path']))
with open(OUT + '/p2x2_trial_ledger.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
# fly ledger
F = collections.defaultdict(list)
for r in rows: F[r['fly']].append(r)
fl = []
for k, v in F.items():
    fl.append(dict(fly=k, substance=CALL.get(k,('',''))[0], confidence=CALL.get(k,('',''))[1], root=v[0]['root'], genotype='/'.join(sorted(set(r['genotype'] for r in v))),
                   indicator='/'.join(sorted(set(r['indicator'] for r in v))), n_trials=len(v),
                   n_saved=sum(r['mat_copies'] > 0 for r in v), n_dup_copies=sum(max(r['mat_copies'] - 1, 0) for r in v),
                   n_used=sum(r['used'] for r in v),
                   strata=';'.join(sorted(set(s for r in v for s in r['stratum'].split(';') if s))),
                   why_not_used=';'.join(sorted(set(r['why_not_used'] for r in v if r['why_not_used'])))))
with open(OUT + '/p2x2_fly_ledger.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(fl[0].keys())); w.writeheader(); w.writerows(fl)
# summary
print('trials on Z:', len(rows), ' flies on Z:', len(fl))
print('trials with >=1 .mat copy:', sum(r['mat_copies'] > 0 for r in rows),
      ' with >1 copy:', sum(r['mat_copies'] > 1 for r in rows), ' stored under 2 paths:', sum(r['stored_paths'] > 1 for r in rows))
print('\nFLIES by root x genotype: total / used / reason not used')
T = collections.defaultdict(lambda: [0, 0, collections.Counter()])
for x in fl:
    t = T[(x['root'], x['genotype'])]; t[0] += 1
    if x['n_used']: t[1] += 1
    else: t[2][x['why_not_used'] or '?'] += 1
for k in sorted(T): print('  %-22s %-6s %3d total  %3d used   unused: %s' % (k[0], k[1], T[k][0], T[k][1], dict(T[k][2])))
print('\nflies used per stratum (lpsp / empty):')
S = collections.defaultdict(lambda: collections.defaultdict(set))
for r in rows:
    if r['used']:
        for s in r['stratum'].split(';'): S[s][r['genotype']].add(r['fly'])
for s in sorted(S): print('  %-34s %3d / %3d' % (s, len(S[s]['lpsp']), len(S[s]['empty'])))
