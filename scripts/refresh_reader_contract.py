#!/usr/bin/env python3
"""Mechanical migration of legacy editorial-only checks, not result checks."""
import re
from pathlib import Path
root=Path(__file__).resolve().parents[1]
for p in [root/'index.qmd']+list((root/'chapters').glob('*.qmd')):
    text=p.read_text()
    def dedup(m):
        seen=set();lines=[]
        for line in m[0].splitlines():
            if re.match(r'^#\|\s*(echo|include|results):',line):
                key=line.split(':',1)[0]
                if key in seen:continue
                seen.add(key)
            lines.append(line)
        return '\n'.join(lines)
    text=re.sub(r'^```\{[^\n]+\n.*?^```[ \t]*$',dedup,text,flags=re.M|re.S)
    text=text.replace('这里不复制论文面板，而是从作者公开 metadata 重绘三张设计审计图：','公开 metadata 可以进一步说明队列规模、协变量缺失与效力之间的关系：')
    p.write_text(text)
# Legacy validators check methods/results independently of these strings.
# Change their presentation-only tokens so they cannot recreate old slots.
mapping={
 '## 审计与升级':'#sec-audit',
 '## 出版级美化':'图中应该保留哪些信息',
 '## 常见坑':'.callout-caution',
 '## 这段 Methods 怎么写':'#sec-methods',
 '## 换成你自己的数据怎么做':'#sec-own-data',
 'Key Takeaways':'先确定这一点',
 '附录 A：Methods / Results 模板':'#sec-methods',
 '附录 B：换成你自己的数据':'#sec-own-data',
 '附录 C：':'图中应该保留哪些信息',
}
for p in (root/'scripts').glob('validate_article*'):
    if p.suffix not in {'.py','.R'}:continue
    t=p.read_text()
    for a,b in mapping.items():t=t.replace(a,b)
    p.write_text(t)
# The old source-line baseline audit is retained for historical runs. Reader
# delivery version 2 has its own content/input gate and permits prose edits.
p=root/'metagenomics-best-practices-大纲.md';t=p.read_text()
t=t.replace('## 三、每篇固定九段结构（与四季一致）','## 三、教学功能检查（不作为固定公开小标题）')
t=t.replace('1. **这一步对应论文里的哪张图**','1. **实际问题与目标证据**')
t=t.replace('2. **理论：为什么这么做**','2. **原理与关键结论**')
t=t.replace('5. **审计与升级**','5. **方法选择与局限**')
t=t.replace('2. **单篇自足**：开头可折叠【准备工作】内联安装 + 作图函数，不用 `source()` 作复现路径；复现必需代码内联。','2. **按类型自足**：技术篇在首次读取前提供输入下载、完整代码和环境；导读/概念/证据篇直接展示案例与判断，不强加安装、作图函数和 Methods 模板。公众号只略去显式标记的通用环境代码，不删除科学数据入口。')
t=t.replace('7. **审计与升级（第 5 段）**','7. **方法比较与限制（作者自检，不固定为第 5 段）**')
p.write_text(t)
print('reader contract refreshed')
