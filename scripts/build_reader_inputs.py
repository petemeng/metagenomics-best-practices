#!/usr/bin/env python3
"""Build chapter-scoped downloadable input inventories, without running tools.

Inventories distinguish input tables, saved analysis results, implementations,
and environment files. Large raw data/databases are never implicitly fetched.
This command does not claim that downloading results recomputes an analysis.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import re
import shutil
from pathlib import Path
from urllib.parse import quote
import yaml

from revise_reader_narrative import split, CONCEPT, INTERPRETATION

ROOT=Path(__file__).resolve().parents[1]
REPO='https://raw.githubusercontent.com/petemeng/metagenomics-best-practices/'
PUBLIC='https://github.com/petemeng/metagenomics-best-practices/'
PRIVATE=re.compile(rb'/media/desk16/|/home/tly9658/|access_token|AppSecret|appsecret|thumb_media_id|refresh_token')
NON_DISTRIBUTED={
    'data/small/13-qc-frozen/multiqc/13-multiqc-report.html',
    'data/small/74-nfcore-mag-workflows-frozen/source/mag-changelog.md',
}
PATH_RE=re.compile(r'(?:data/small|results|env|scripts)/[A-Za-z0-9_./+\-]+')


def sha(p,kind='sha256'):
    h=hashlib.new(kind)
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def files_below(path):
    if path.is_file():return [path]
    if path.is_dir():return sorted(p for p in path.rglob('*') if p.is_file() and not p.is_symlink())
    return []


def candidates(root,number,text):
    selected=set()
    anchors={25:'25-hmp-fig3-original.jpg',28:'28-wirbel-fig3-original.png'}
    if number in anchors:
        selected.add(root/'figures'/anchors[number])
    for pattern in [f'{number:02}-*', f'{number:02}.*']:
        for p in (root/'data/small').glob(pattern):selected.update(files_below(p))
    # Explicit cross-chapter inputs are downloaded too: a path reference is
    # not permission to require readers to execute the previous chapter.
    for token in PATH_RE.findall(text):
        path=root/token.rstrip('.')
        if path==root/'data/small' or path==root/'results':continue
        selected.update(files_below(path))
    for p in (root/'scripts').glob(f'*article{number:02}*'):
        if p.suffix in {'.R','.py','.sh'}:selected.add(p)
    # A few older scripts use unpadded article numbers.
    if number<10:
        for p in (root/'scripts').glob(f'*article{number}_*'):
            if p.suffix in {'.R','.py','.sh'}:selected.add(p)
    for p in (root/'results').glob(f'{number:02}-*'):
        selected.update(files_below(p))
    for _ in range(4):
        before=len(selected)
        for p in list(selected):
            if p.suffix not in {'.R','.py','.sh'}:continue
            content=p.read_text(errors='replace')
            for name in PATH_RE.findall(content):
                child=root/name.rstrip('.')
                if child.is_file() and str(child.relative_to(root)).startswith(('scripts/','env/')):selected.add(child)
            for name in re.findall(r'^from ([A-Za-z0-9_]+) import',content,re.M):
                child=root/'scripts'/f'{name}.py'
                if child.is_file():selected.add(child)
        if before==len(selected):break
    return selected


def build(root,ref):
    manifest=yaml.safe_load((root/'tutorial.yaml').read_text());report=[]
    for chapter in manifest['series']['chapters']:
        n=int(chapter['number']);qmd=root/chapter['file'];front,body=split(qmd.read_text())
        if n in CONCEPT|INTERPRETATION:continue
        folder=root/f'examples/{n:02}';folder.mkdir(parents=True,exist_ok=True)
        rows=[];omitted=[]
        for p in sorted(candidates(root,n,body)):
            rel=p.relative_to(root).as_posix()
            if rel in NON_DISTRIBUTED:
                omitted.append({'path':rel,'reason':'third-party report excluded from public source'})
                continue
            if '__pycache__' in rel or (p.suffix in {'.png','.jpg','.jpeg','.pdf','.tiff','.log','.pyc'} and not rel.startswith('figures/')):continue
            if rel.startswith('results/') and p.suffix not in {'.tsv','.csv','.gz'}:continue
            if p.stat().st_size>95*1024*1024:
                omitted.append({'path':rel,'reason':'large upstream artifact; use original data instructions'});continue
            # Text receipts containing workstation identity are not public
            # teaching inputs. Do not silently redact scientific data bytes.
            if p.suffix in {'.json','.tsv','.txt','.yaml','.yml','.csv'} and PRIVATE.search(p.read_bytes()):
                omitted.append({'path':rel,'reason':'private execution metadata'});continue
            role='input'
            src=rel
            if rel.startswith('scripts/'):role='implementation'
            elif rel.startswith('figures/'):role='illustration'
            elif rel.startswith('env/'):role='environment'
            elif rel.startswith('results/'):
                role='saved_result';src=f'data/reader-results/{n:02}/'+rel.removeprefix('results/')
                dest=root/src;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
            rows.append({'path':rel,'source':src,'role':role,'bytes':p.stat().st_size,'md5':sha(p,'md5'),'sha256':sha(p)})
        if not rows:raise RuntimeError(f'No inputs for {n}')
        with (folder/'files.tsv').open('w') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
        size=sum(r['bytes'] for r in rows)/1024/1024
        readme=f'''# {chapter['title']}

数据与代码清单：[files.tsv](files.tsv)。约 {size:.1f} MiB；不含原始 FASTQ 或大型参考数据库。

清单中 `input` 是用于本例的公开输入或有来源记录的处理后表格；`saved_result` 是已经计算的结果，便于核对图表，不代表重新拟合模型；`implementation` 是分析脚本；`environment` 是环境记录。

下载后按[本章完整说明]({PUBLIC}blob/{ref}/{chapter['file']})运行。原始数据的研究出处、样本筛选、统计单位和局限见正文。重计算流程还需要正文指定的软件、原始输入和数据库，下载这份清单不等于完成上游重跑。

在空文件夹中下载并运行 `download.R`，文件会保留正文所用的相对目录。已存在且与清单不一致的文件会报错，不会覆盖。
'''
        (folder/'README.md').write_text(readme)
        code=f'''options(timeout = max(600, getOption("timeout")))
base_url <- "{REPO}{ref}/"
files <- read.delim(paste0(base_url, "examples/{n:02}/files.tsv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(!grepl("(^/|(^|/)\\\\.\\\\.(/|$))", files$path)))
for (i in seq_len(nrow(files))) {{
  target <- files$path[i]
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(target)) {{
    temporary <- tempfile(tmpdir = dirname(target))
    download.file(paste0(base_url, files$source[i]), temporary, mode = "wb", quiet = TRUE)
    if (unname(tools::md5sum(temporary)) != files$md5[i]) {{
      unlink(temporary)
      stop("Downloaded file checksum mismatch: ", target)
    }}
    if (!file.rename(temporary, target)) stop("Cannot save: ", target)
  }}
  if (unname(tools::md5sum(target)) != files$md5[i]) {{
    stop("Existing input differs from this version: ", target,
         ". Use a new folder; do not overwrite your own data.")
  }}
}}
message("Inputs downloaded; saved results have NOT been recomputed.")
'''
        (folder/'download.R').write_text(code)
        # Short reader-visible download; verification is in the complete
        # download script, not duplicated as generic WeChat plumbing.
        entry=f'''## 下载本例数据与分析代码 {{#sec-reader-inputs}}

[数据与代码清单]({PUBLIC}tree/{ref}/examples/{n:02})列出了本例的输入表、分析脚本和环境文件（约 {size:.1f} MiB）。原始 FASTQ 和大型数据库按下文的数据来源另行获取。清单中的已计算结果仅供核对；读取结果表不等于重新运行统计模型或上游工具。

在新建文件夹中运行下面的 R 代码，下载完成后即可按正文中的相对路径读取文件。需要核验文件版本时，可使用[带完整性检查的下载脚本]({REPO}{ref}/examples/{n:02}/download.R)。

```{{r}}
#| label: reader-input-download-{n:02}
#| eval: false
options(timeout = 600)
base_url <- "{REPO}{ref}/"
files <- read.delim(paste0(base_url, "examples/{n:02}/files.tsv"))
for (i in seq_len(nrow(files))) {{
  path <- files$path[i]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) {{
    download.file(paste0(base_url, files$source[i]), path, mode = "wb")
  }}
}}
```

'''
        # Input setup must precede even the formerly hidden loader chunks.
        body=re.sub(r'## 下载本例数据与分析代码 \{#sec-reader-inputs\}.*?(?=^## |^```\{r\}\n#\| label: (?:load-|setup-)|\Z)','',body,flags=re.M|re.S) if '#sec-reader-inputs}' in body else body
        # Insert after the opening recommendation but before the first code.
        m=re.search(r'^```\{r\}',body,re.M)
        h=re.search(r'^## ',body,re.M)
        pos=min(x.start() for x in [m,h] if x is not None)
        body=body[:pos]+entry+body[pos:]
        qmd.write_text(front+body)
        report.append({'chapter':n,'files':len(rows),'size_mib':round(size,2),'omitted':omitted})
    out=root/'qa/reader-revision';out.mkdir(parents=True,exist_ok=True)
    (out/'input-inventories.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'chapters':len(report),'files':sum(x['files'] for x in report),'total_mib':round(sum(x['size_mib'] for x in report),1)}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,default=ROOT);p.add_argument('--ref',required=True);a=p.parse_args();build(a.project_root.resolve(),a.ref)
