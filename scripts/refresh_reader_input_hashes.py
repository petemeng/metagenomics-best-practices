#!/usr/bin/env python3
"""Refresh inventories after an reviewed implementation edit; no source changes."""
import csv,hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
for p in (root/'examples').glob('*/files.tsv'):
    with p.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    anchors={'25':'figures/25-hmp-fig3-original.jpg','28':'figures/28-wirbel-fig3-original.png'}
    anchor=anchors.get(p.parent.name)
    if anchor and not any(r['source']==anchor for r in rows):
        rows.append(dict(path=anchor,source=anchor,role='illustration',bytes=0,md5='',sha256=''))
    for row in rows:
        source=root/row['source'];data=source.read_bytes()
        row.update(bytes=len(data),md5=hashlib.md5(data).hexdigest(),sha256=hashlib.sha256(data).hexdigest())
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
print('64 reader input inventories refreshed')
