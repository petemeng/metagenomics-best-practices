#!/usr/bin/env python3
"""Run selected complete analyses from a fresh directory and input inventory.

Use --base-url for an immutable public commit; by default test the identical
files through file URLs before publication. Download and analysis results are
reported separately. No heavy upstream workflow is launched.
"""
import argparse,json,os,re,shutil,subprocess,tempfile,time
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,default=Path('.'));p.add_argument('--chapter',type=int,choices=[25,28,71],required=True);p.add_argument('--base-url');p.add_argument('--output',type=Path,required=True);a=p.parse_args()
root=a.project_root.resolve();work=Path(tempfile.mkdtemp(prefix=f'metagenome-reader-cold-{a.chapter:02}-'))
base=a.base_url or root.as_uri()+'/'
code=(root/f'examples/{a.chapter:02}/download.R').read_text()
code=re.sub(r'^base_url <- .*$',f'base_url <- "{base}"',code,flags=re.M)
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
started=time.monotonic();log=work/'download.log'
with log.open('w') as out:
    download=subprocess.run(['Rscript','-'],input=code,text=True,cwd=work,env=env,stdout=out,stderr=subprocess.STDOUT,timeout=600)
result={'chapter':a.chapter,'work_dir':str(work),'input_transport':'public HTTPS' if a.base_url else 'local file URL','download_returncode':download.returncode,'download_log':str(log),'analysis_recomputed':False}
if download.returncode==0:
    if a.chapter==25:
        cmd=['Rscript','scripts/validate_article25_composition_core.R','--project-root','.','--input-dir','data/small/25-composition-core','--output-dir','results/cold-analysis','--figure-dir','figures']
    elif a.chapter==28:
        # The manuscript is needed only by the report's editorial assertions;
        # all scientific inputs come from the chapter's download inventory.
        qmd=next((root/'chapters').glob('28-*.qmd'));shutil.copy2(qmd,work/'chapter.qmd')
        cmd=['Rscript','scripts/validate_article28_cross_cohort.R','--project-root','.','--input-dir','data/small/28-cross-cohort','--output-dir','results/cold-analysis','--figure-dir','figures','--chapter','chapter.qmd']
    else:
        cmd=['Rscript','scripts/run_article71_sem_models.R','--input-dir','data/small/71-sem-frozen','--output-dir','results/cold-analysis']
    with (work/'analysis.log').open('w') as out:
        run=subprocess.run(cmd,cwd=work,env=env,stdout=out,stderr=subprocess.STDOUT,timeout=900)
    result.update(analysis_returncode=run.returncode,analysis_recomputed=run.returncode==0,analysis_log=str(work/'analysis.log'))
result['seconds']=round(time.monotonic()-started,1);result['status']='passed' if result['analysis_recomputed'] else 'failed'
a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));raise SystemExit(0 if result['status']=='passed' else 1)
