#!/usr/bin/env python3
"""Reader-delivery gate: source, data inventories, and optional public payloads.

Static delivery validation is deliberately not described as a cold analysis
rerun. Scientific runs retain separate per-analysis reports.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
import yaml
from lxml import html

FORBIDDEN=re.compile(r'真实论文为[^\n]*证据起点|旧草稿|整仓库离线渲染|这里不(?:复制|转载)论文|附录 A：Methods / Results 模板')
PRIVATE=re.compile(r'/media/desk16/|/home/tly9658/|wx[0-9a-f]{16}|(?:access_token|appsecret)\s*[:=]\s*["\'][A-Za-z0-9]')

def digest(p,kind='sha256'):
    h=hashlib.new(kind)
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def has_decision_recommendation(body):
    prose = re.sub(r'^```.*?^```[ \t]*$', '', body, flags=re.M | re.S)
    recommendations = re.findall(
        r'::: \{\.callout-(?:tip|important)[^\n]*\}\n(.*?)\n:::',
        prose, flags=re.S)
    return any(len(re.sub(r'\s+', '', text)) >= 20 for text in recommendations)

def audit(root,bundle=None,check_hashes=True):
    manifest=yaml.safe_load((root/'tutorial.yaml').read_text());errors=[];reports=[]
    chapters=manifest['series']['chapters']
    if [int(x['number']) for x in chapters]!=list(range(1,78)):errors.append('chapter sequence must be 1..77')
    for c in chapters:
        n=int(c['number']);p=root/c['file'];text=p.read_text();front,body=text[4:].split('\n---\n',1);fm=yaml.safe_load(front);bad=[]
        if fm.get('author')!='Peter':bad.append('author')
        if fm.get('format',{}).get('html',{}).get('number-sections'):bad.append('section numbering')
        if re.match(r'第\s*\d+\s*篇',fm['title']):bad.append('episode prefix in source title')
        if FORBIDDEN.search(body):bad.append('author/template narration')
        if PRIVATE.search(body):bad.append('private content')
        # A topic-specific recommendation callout, not a mandatory stock title.
        # Inspect prose (excluding code) so a literal in a script cannot satisfy it.
        if not has_decision_recommendation(body):
            bad.append('missing decision recommendation')
        if fm.get('reader-mode') not in {'analysis','upstream','evidence'}:bad.append('reader category')
        inv=root/f'examples/{n:02}/files.tsv';files=[]
        if fm.get('reader-mode')!='evidence':
            if '#sec-reader-inputs}' not in body or 'download.file(' not in body:bad.append('missing input download')
            if not re.search(r'base_url <- "https://raw\.githubusercontent\.com/petemeng/metagenomics-best-practices/[0-9a-f]{40}/"',body):bad.append('non-immutable data entry')
            if not re.search(r'/blob/[0-9a-f]{40}/examples/'+f'{n:02}/files.tsv',body):bad.append('non-immutable input inventory link')
            if not re.search(r'/[0-9a-f]{40}/examples/'+f'{n:02}/download.R',body):bad.append('non-immutable verified downloader link')
            if not inv.is_file():bad.append('missing input inventory')
            else:
                with inv.open() as f:files=list(csv.DictReader(f,delimiter='\t'))
                for item in files:
                    src=root/item['source'];dest=Path(item['path'])
                    if dest.is_absolute() or '..' in dest.parts:bad.append('unsafe input destination');continue
                    if not src.is_file():bad.append('missing asset: '+item['source']);continue
                    if check_hashes and (digest(src)!=item['sha256'] or digest(src,'md5')!=item['md5']):bad.append('asset changed: '+item['source'])
                    if item['role'] not in {'input','saved_result','implementation','environment','illustration'}:bad.append('unlabelled asset role')
        if n in {1,2} and ('install.packages(' in body or '#sec-methods}' in body):bad.append('forced overview bootstrap/methods')
        if bundle:
            payload=json.loads((bundle/f'{n:02}/draft.json').read_text());content=payload['content'];doc=html.fromstring(content)
            if not payload['title'].startswith(f'宏基因组最佳实践｜{n}. '):bad.append('draft title')
            if payload['author']!='Peter':bad.append('draft author')
            if not re.search(r'/blob/[0-9a-f]{40}/'+re.escape(c['file'])+'$',payload['content_source_url']):bad.append('non-immutable/non-chapter source URL')
            if FORBIDDEN.search(doc.text_content()) or PRIVATE.search(content):bad.append('public boundary')
            code='\n'.join(x.text_content() for x in doc.xpath('//pre'))
            if fm.get('reader-mode')!='evidence' and 'download.file(' not in code:bad.append('download removed by renderer')
            if n==28 and not re.search(r'cohort28\s*<-',code):bad.append('undefined cohort28')
            if n==71 and not (re.search(r'frozen\s*<-',code) and re.search(r'read_frozen\s*<-\s*function',code)):bad.append('undefined SEM inputs')
            if re.search(r'(theme_pub|save_pub)\s*<-\s*function|install\.packages\(',code):bad.append('generic bootstrap in WeChat')
            if not doc.xpath('//img'):bad.append('missing result figure')
        errors.extend(f'{n:02}: {e}' for e in bad)
        reports.append({'chapter':n,'mode':fm.get('reader-mode'),'input_files':len(files),'passed':not bad})
    return {'status':'failed' if errors else 'passed','scope':'editorial-and-input-delivery; not a full scientific rerun','chapters':len(reports),'passed_chapters':sum(r['passed'] for r in reports),'errors':errors,'reports':reports}

def main():
    p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,default=Path('.'));p.add_argument('--bundle',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--skip-hashes',action='store_true');a=p.parse_args()
    r=audit(a.project_root.resolve(),a.bundle,not a.skip_hashes);a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k!='reports'},ensure_ascii=False));return 0 if r['status']=='passed' else 1
if __name__=='__main__':raise SystemExit(main())
