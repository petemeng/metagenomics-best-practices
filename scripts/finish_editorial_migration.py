#!/usr/bin/env python3
"""Category-specific presentation fixes; never recompute biological pipelines."""
import re
from pathlib import Path
from revise_reader_narrative import split, sections, FENCE, CONCEPT, INTERPRETATION

ROOT=Path(__file__).resolve().parents[1]

def remove_sections(body, anchors):
    for a,b,h in reversed(sections(body)):
        if any('#'+x+'}' in h for x in anchors):body=body[:a]+body[b:]
    return body

for p in [ROOT/'index.qmd']+sorted((ROOT/'chapters').glob('*.qmd')):
    n=1 if p.name=='index.qmd' else int(p.name[:2]);front,body=split(p.read_text())
    if n in {1,2}:
        body=remove_sections(body,{'sec-setup','sec-methods','sec-polish'})
        body=FENCE.sub(lambda m:'' if re.match(r'```\{?(?:r|bash)',m[0]) else m[0],body)
        body=re.sub(r'^> \*\*(?:学习目标|真实数据|真实锚点|预计资源)：.*\n','',body,flags=re.M)
        body=body.replace('## 读取作者公开元数据生成的队列表','## 八个队列，怎样区分发现与验证')
        body=body.replace('## 读取测量边界证据表','## 需要绝对数量或功能活动时，还缺什么？')
        body=body.replace('## 读取并审计两个输入表','## 怎样读这张路线选择表')
        body=body.replace('## 转成长表，避免把层级混成一个变量','## 样本、基因和基因组，不是同一种分析单位')
        body=body.replace('## 从真实指标构建路线锚点','## 三条路线分别适合哪些问题？')
        body=body.replace('这里重绘两张决策图：','两个研究对应的分析对象如下：')
        body=re.sub(r'\n{4,}','\n\n\n',body)
        front=front.replace('eval: true','eval: false')
    elif n in CONCEPT:
        # Retain actual computed evidence tables, not the author-facing R
        # plumbing that happens to produce their layout.
        hidden=[]
        for a,b,h in reversed(sections(body)):
            if '#sec-preparation}' in h or '#sec-setup}' in h:
                for m in FENCE.finditer(body[a:b]):
                    if m[0].startswith('```{r}'):
                        block=m[0]
                        if not re.search(r'#\|\s*eval:\s*false',block):
                            hidden.append(block.replace('```{r}\n','```{r}\n#| include: false\n',1))
                body=body[:a]+body[b:]
        body='\n'+'\n\n'.join(reversed(hidden))+'\n\n'+body
        def evidence_code(m):
            block=m[0]
            if block.startswith('```{r}') and '#| include: false' not in block:
                if not re.search(r'^#\| echo:', block, re.M):
                    block=block.replace('```{r}\n','```{r}\n#| echo: false\n',1)
            return block
        body=FENCE.sub(evidence_code,body)
    if n in INTERPRETATION:
        # This revision is confined to interpreting existing aggregate
        # evidence. Do not add executable discovery, host-inference,
        # virulence or mobilization workflows or targeted sequence inputs.
        body=remove_sections(body,{'sec-preparation','sec-setup'})
        body=FENCE.sub(lambda m:'' if re.match(r'```\{?(?:bash|sh|python|r)(?:\W|$)',m[0]) else m[0],body)
        body=re.sub(r'\n{4,}','\n\n\n',body)
    # A hidden input-definition chunk is unsafe to copy from the rendered
    # article: make the small loader visible, keeping formatting-only code
    # out of sight where it is genuinely not part of the analysis.
    if n in {27,28,29,63,64,71}:
        body=body.replace('#| include: false','#| echo: true\n#| results: hide')
    p.write_text(front+body)

# In the composition tutorial the main figures formerly preceded every data
# read. Place each one next to the operation and explanation it summarizes.
p=ROOT/'chapters/25-composition-core-microbiome.qmd';t=p.read_text()
targets={
 '25-multirank-mean-composition.png':r'^## 把 prevalence',
 '25-individual-stool-composition.png':r'^## 把 prevalence',
 '25-prevalence-abundance.png':r'^## bootstrap',
 '25-core-membership-sensitivity.png':r'^## 锁定结果',
}
for name,target in targets.items():
    m=re.search(r'^!\[[^\n]+'+re.escape(name)+r'[^\n]*\n',t,re.M)
    if not m:raise RuntimeError(name)
    block=m[0];t=t[:m.start()]+t[m.end():]
    h=re.search(target,t,re.M)
    if not h:raise RuntimeError(target)
    t=t[:h.start()]+block+'\n'+t[h.start():]
t=t.replace('完整 92 个成员保存在结果表。','共有 92 个成员满足主核心定义。')
p.write_text(t)
print('category-specific revisions complete')
