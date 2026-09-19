"""Reject known reader-facing prose regressions without scanning code examples."""
import argparse
import json
from pathlib import Path
import re

def prose_errors(text):
    text=re.sub(r'\A---\n.*?\n---\n','',text,flags=re.S)
    text=re.sub(r'^(`{3,}|~{3,})[^\n]*\n.*?^\1\s*$','',text,flags=re.M|re.S)
    text=re.sub(r'`+[^`\n]*`+','',text)
    rules={
        'legacy_math_delimiters':r'\\[\[\]()]',
        'wrong_wetland_groups':r'conventional\s*[/、,]\s*intensive|intensive\s*[/、,]\s*temporal',
        'maintainer_prose':r'教程合同|这里不复制论文原图|不是伪装成|可复制阅读框架|审计与升级',
        'legacy_public_episode':r'16S最佳实践[（(]\s*\d+\s*/\s*55',
        'formula_typesetting_leak':r'\(sum\(completeness-5\\times contamination\)\)',
    }
    return [name for name,pattern in rules.items() if re.search(pattern,text)]

def main():
    p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,default=Path('.'));p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    files=[a.project_root/'index.qmd',*sorted((a.project_root/'chapters').glob('*.qmd'))]
    failures=[{'source':str(f.relative_to(a.project_root)),'errors':errors} for f in files if (errors:=prose_errors(f.read_text()))]
    result={'status':'failed' if failures else 'passed','chapters':len(files),'failures':failures}
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');print(json.dumps(result,ensure_ascii=False));return bool(failures)

if __name__=='__main__':raise SystemExit(main())
