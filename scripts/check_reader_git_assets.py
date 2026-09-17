#!/usr/bin/env python3
"""Ensure every downloadable asset really exists at its immutable Git revision."""
import argparse,csv,hashlib,json,re,subprocess
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--project-root',type=Path,default=Path('.'))
p.add_argument('--git-repository',type=Path,default=Path('.'))
p.add_argument('--asset-ref',required=True)
a=p.parse_args();root=a.project_root.resolve()
if not re.fullmatch('[a-f0-9]{40}',a.asset_ref):raise SystemExit('A full immutable commit is required')
tree={}
for line in subprocess.check_output(['git','ls-tree','-rz',a.asset_ref],cwd=a.git_repository).split(b'\0'):
    if not line:continue
    meta,path=line.split(b'\t',1);tree[path.decode()]=meta.split()[2].decode()
seen=set();rows=0;inventories=sorted((root/'examples').glob('*/files.tsv'))
def verify_blob(path):
    rel=path.relative_to(root).as_posix()
    if rel not in tree:raise SystemExit('Missing public asset: '+rel)
    data=path.read_bytes()
    blob=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    if blob!=tree[rel]:raise SystemExit('Public asset differs from inventory version: '+rel)
    return data
for inv in inventories:
    verify_blob(inv)
    with inv.open() as f:
        for item in csv.DictReader(f,delimiter='\t'):
            rows+=1;rel=item['source']
            if rel in seen:continue
            data=verify_blob(root/rel)
            if hashlib.sha256(data).hexdigest()!=item['sha256']:raise SystemExit('Incorrect inventory hash: '+rel)
            seen.add(rel)
print(json.dumps({'status':'passed','immutable_asset_commit':a.asset_ref,'inventories':len(inventories),'inventory_rows':rows,'unique_git_assets':len(seen)}))
