# Step 3: z_scan.tsv -> z_trials.json, one record per Z: trial folder (root, date, fly, genotype tag,
# indicator, folder status ok/open loop/exclusions/pilot (dont use)). Run from the work dir.
import csv,re,collections,json
B=chr(92)
trials={}  # z trial dir -> info
for line in open('z_scan.tsv',encoding='utf-8',errors='replace'):
    p,sz=line.rstrip('\n').split('\t'); parts=p.split('/')
    name=parts[-1]
    if name=='trialSettings.csv' and parts[-2]=='csv': d='/'.join(parts[:-2])
    elif name.lower().endswith('ficTracData_DAQ.mat'.lower()) or name=='mask.mat' or name.endswith('.tif'): d='/'.join(parts[:-1])
    elif name.startswith('registration'): d='/'.join(parts[:-1])
    elif parts[-2].startswith('registration'): d='/'.join(parts[:-2])
    else: continue
    t=trials.setdefault(d,{'tif':0,'daq':0,'reg':set(),'mask':0})
    if name.endswith('.tif'): t['tif']+=1
    elif 'ficTracData_DAQ' in name: t['daq']=1
    elif name=='mask.mat': t['mask']=1
    elif parts[-2].startswith('registration') and name.endswith('.mat'): t['reg'].add(name)
# keep real trial folders: has daq or tif
trials={d:t for d,t in trials.items() if t['daq'] or t['tif']}
def info(d):
    parts=d.split('/'); root=parts[0]; tname=parts[-1]
    m=re.match(r'(\d{8})-(\d+)',tname); date=m.group(1) if m else None
    flyf=[x for x in parts[1:-1] if x.lower().startswith('fly')]
    fly=re.match(r'fly\s*(\d+)',flyf[-1],re.I).group(1) if flyf else None
    low=d.lower()
    status='open loop' if 'open loop' in low else ('exclusions' if 'exclusion' in low else ('pilot (dont use)' if 'dont use' in low else 'ok'))
    tn=tname.lower(); ff=(flyf[-1].lower() if flyf else '')
    g=''
    if 'empty' in tn or 'empty' in ff or 'control' in ff: g='empty'
    elif 'lpsp' in tn or 'lpsp' in ff or 'th_p2x2' in tn: g='lpsp'
    elif 'blind' in tn and fly and date and date.startswith('2024') and '/todo/blind/' in d+'/': g='lpsp'
    ind='syt8s' if 'syt8s' in tn else ('8m' if re.search(r'8m',tn) else ('7f' if '7f' in tn else ('7f?' if date and date.startswith('2024') else '?')))
    return dict(root=root,tname=tname,date=date,fly=fly,flyfolder=flyf[-1] if flyf else '',status=status,geno=g,ind=ind)
rows=[]
for d,t in trials.items():
    i=info(d); i.update(path=d,tif=t['tif'],daq=t['daq'],mask=t['mask'],reg=';'.join(sorted(t['reg'])))
    if not i['date']:  # e.g. "20260708-bailed_fly_...": take the date folder
        i['date']=next(x for x in d.split('/') if re.fullmatch(r'\d{8}',x))
    # root + trial name, "-4-epg" normalized to "-4_epg": the 2024 blinded trials reuse
    # trial numbers within a date (20240624-1_blind_3 / _blind_4), so the number alone isn't unique
    i['key']=i['root']+'|'+re.sub(r'^(\d{8}-\d+)[-_]',lambda m: m.group(1)+'_',i['tname'].lower())
    rows.append(i)
json.dump(rows,open('z_trials.json','w'))
print('trial folders on Z:',len(rows))
c=collections.Counter((r['root'],r['status']) for r in rows); [print(' ',v,k) for k,v in sorted(c.items())]
keys=collections.defaultdict(list)
for r in rows: keys[r['key']].append(r)
dup=[v for v in keys.values() if len(v)>1]
print('unique trial keys:',len(keys),' keys with >1 folder on Z:',len(dup))
cc=collections.Counter(tuple(sorted(x['root']+':'+x['status']+':'+('/'.join(x['path'].split('/')[1:-1][:2])) for x in v)) for v in dup)
for k,v in cc.most_common(12): print('  ',v,k)
nod=[r for r in rows if not r['date']]; print('no date in trial name:',len(nod),[r['path'] for r in nod[:5]])
