# Step 4: match p2x2_inventory.csv (.mat copies) to z_trials.json by root + normalized trial name;
# prints coverage / duplicates / missing and writes p2x2_trial_ledger.csv (first pass). Run from the work dir.
import csv,re,collections,json
B=chr(92)
def norm(t): return re.sub(r'^(\d{8}-\d+)[-_]', lambda m: m.group(1)+'_', t.lower())
rows=json.load(open('z_trials.json'))
for r in rows:
    if r['geno']=='' and 'blind' in r['tname'].lower(): r['geno']='lpsp'; r['geno_rule']='blinded, untagged fly'
    else: r['geno_rule']='tag'
    r['low']='low' in r['tname'].lower()
    r['dark']='dark' in r['tname'].lower()
Z={r['key']:r for r in rows}
inv=list(csv.DictReader(open('p2x2_inventory.csv')))
copies=collections.defaultdict(list); unmatched=[]
for r in inv:
    p=[x for x in r['meta'].split(B) if x]
    if p[-1].startswith('registration'): p=p[:-1]
    key=p[2]+'|'+norm(p[-1])
    if key in Z: copies[key].append({'file':r['file'],'path':B.join(p[2:]),'nstims':r['nstims'],'fr':r['fr']})
    else: unmatched.append(B.join(p))
print('mat trial copies:',len(inv),' matched:',len(inv)-len(unmatched),' unmatched:',unmatched[:3])
pv=collections.Counter(len(set(c['path'] for c in v)) for v in copies.values())
print('Z trials with >=1 .mat copy: %d / %d ; distinct stored paths per trial: %s'%(len(copies),len(rows),dict(pv)))
nf=collections.Counter(len(set(c['file'] for c in v)) for v in copies.values()); print('number of .mat files holding each trial:',dict(sorted(nf.items())))
for r in rows:
    r['in_mat']=r['key'] in copies; r['mat_files']=';'.join(sorted(set(c['file'] for c in copies.get(r['key'],[]))))
    r['stored_paths']=len(set(c['path'] for c in copies.get(r['key'],[])))
    r['flykey']=r['root']+'|'+r['date']+'|'+('fly '+r['fly'] if r['fly'] else r['tname'])
# experiment label
def exp(r):
    if r['root']=='lpsp_p2x2': return '2024 EPG 7f' + (' (low ATP)' if r['low'] else '')
    if r['root']=='lpsp_p2x2_redo': return '2026 redo EPG syt8s'
    if r['root']=='lpsp_p2x2_reredo': return '2026 reredo EPG '+r['ind']
    if r['root']=='lpsp_p2x2_syt8s': return '2026 syt8s EPG syt8s'
    return '2026 walking EPG '+r['ind']
for r in rows: r['exp']=exp(r)
fl=collections.defaultdict(list)
for r in rows: fl[r['flykey']].append(r)
T=collections.defaultdict(lambda: collections.Counter())
for k,v in fl.items():
    st=v[0]['status']; key=(v[0]['exp'],st)
    g='/'.join(sorted(set(r['geno'] for r in v)))
    nin=sum(r['in_mat'] for r in v)
    cov='full' if nin==len(v) else ('none' if nin==0 else 'partial')
    T[key][(g,cov)]+=1
print('\nexperiment | folder status : (genotype, .mat coverage) -> flies')
for k in sorted(T): print(' ',k,dict(T[k]))
print('\nTRIALS per experiment|status: in .mat / total; trials stored under 2 paths')
c=collections.Counter(); t=collections.Counter(); d=collections.Counter()
for r in rows: k=(r['exp'],r['status']); t[k]+=1; c[k]+=r['in_mat']; d[k]+=r['stored_paths']>1
for k in sorted(t): print('  ',k,c[k],'/',t[k],' dup-path:',d[k])
with open('p2x2_trial_ledger.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['exp','status','root','date','flykey','flyfolder','tname','geno','geno_rule','ind','low','dark','tif','daq','mask','reg','in_mat','mat_files','stored_paths','path'])
    w.writeheader(); [w.writerow({k:r[k] for k in w.fieldnames}) for r in sorted(rows,key=lambda r:(r['root'],r['date'],r['flykey'],r['tname']))]
missing=[r for r in rows if not r['in_mat'] and r['status']=='ok']
print('\nok-status trials with NO .mat copy:',len(missing))
mc=collections.Counter((r['exp'],r['flykey'].split('|',1)[1].rsplit('|',1)[0]) for r in missing)
for k,v in sorted(mc.items()): print('  ',v,k)
print('  of these, have registration output:',sum(bool(r['reg']) for r in missing),' have mask:',sum(r['mask'] for r in missing))
