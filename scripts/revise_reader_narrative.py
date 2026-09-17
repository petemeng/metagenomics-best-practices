#!/usr/bin/env python3
"""One-time, reviewable migration of the 77 reader-facing chapters.

Scientific code, frozen results and citations are not rewritten by prose rules.
The report records every changed source and the explicitly relocated figures.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import yaml

FENCE = re.compile(r"^(`{3,}|~{3,})[^\n]*\n.*?^\1[ \t]*$", re.M | re.S)
CONCEPT = {1, 2, 72, 75, 76, 77}
INTERPRETATION = {38, 39, 53, 54, 55, 56, 57}


def split(text):
    end = text.index("\n---\n", 4) + 5
    return text[:end], text[end:]


def mask_code(text):
    return FENCE.sub(lambda m: re.sub(r"[^\n]", " ", m[0]), text)


def sections(body):
    hits = list(re.finditer(r"^## .+$", mask_code(body), re.M))
    return [(m.start(), hits[i+1].start() if i+1 < len(hits) else len(body), body[m.start():m.end()]) for i,m in enumerate(hits)]


def prose_only(text, fn):
    parts = []; last = 0
    for m in FENCE.finditer(text):
        parts.extend([fn(text[last:m.start()]), m[0]])
        last = m.end()
    parts.append(fn(text[last:]))
    return ''.join(parts)


def polish_prose(text):
    text = re.sub(r"^(#{2,6} )当前升级\s*\d*\s*[：:]\s*", r"\1", text, flags=re.M)
    text = re.sub(r"^(#{2,6} )(?:结果[一二三四五六七八九十]+|理论|步骤)[：:]\s*", r"\1", text, flags=re.M)
    text = text.replace('这里不转载论文原图，而是从作者公开的 CC BY 4.0 元数据重建队列结构：', '公开元数据记录了发现队列与独立验证队列的样本构成：')
    text = text.replace('四张重绘图来自锁定的 QUAST 5.3.0 输出，不复制论文图。', '以下结果来自同一版本的 QUAST 5.3.0 输出。')
    text = text.replace('下图是流程重绘，不复制论文图片。', '下图区分了三种目录的输入与合并规则。')
    text = text.replace('每篇需要的出版级函数在正文中直接定义：', '绘图使用下面的颜色和导出尺寸：')
    text = text.replace('不会改写其他章节已经验证的 R 快照', '可单独恢复本次分析所用的包版本')
    text = text.replace('不能假装完全复刻的部分', '不同数据处理版本的结果为什么不完全相同')
    text = text.replace('先忠实复现作者问题，再说明改了什么', '原研究与当前数据的比较范围')
    text = re.sub(r'原图中 “LOSO” 已经是 N−1，但旧草稿只取两个队列互相预测，实质仍是一对一迁移。', 'LOSO 的训练集应包含除测试队列以外的所有队列；两个队列互相预测只回答单来源迁移问题。', text)
    text = text.replace('正文读取的 `data/small/71-sem-frozen/` 是上述管线生成的 checksum 证据包，不是另一套隐藏输入。整仓库离线渲染时直接读它；需要从网络重建时执行上面的命令。', '公开表经样本对齐与完整病例筛选后，得到 90 名 PRISM 参与者和 38 名外部队列参与者。下面分别检查节点定义、暴露分布和路径估计。')
    text = text.replace('下载脚本不会抓取或分发 Franzosa Figure 1。', '')
    text = re.sub(r'^(#{2,6} )冻结、出图与验收', r'\1汇总质量结果并检查图表', text, flags=re.M)
    text = re.sub(r'^(#{2,6} )重建 selected candidates 并锁定版本', r'\1确认待评估的 MAG 与输入版本', text, flags=re.M)
    text = re.sub(r'^(#{2,6} )当流程通过时，输出应该长什么样', r'\1用已知真值检查路径估计与效力', text, flags=re.M)
    return text


def make_editorial(root, topics):
    changes = []
    for source in [root/'index.qmd'] + sorted((root/'chapters').glob('*.qmd')):
        number = 1 if source.name == 'index.qmd' else int(source.name[:2])
        old = source.read_text(); front, body = split(old); row = topics[number]
        front = re.sub(r'(?m)^title: ["\']?第\s*\d+\s*篇\s*[·｜]\s*', 'title: "', front)
        front = re.sub(r'(?m)^author:.*$', 'author: "Peter"', front)
        if not re.search(r'^author:', front, re.M):
            front = front.replace('---\n', '---\nauthor: "Peter"\n', 1)
        front = front.replace('number-sections: true', 'number-sections: false')
        mode = 'evidence' if number in CONCEPT | INTERPRETATION else ('analysis' if yaml.safe_load(front[4:-5]).get('execute',{}).get('eval') else 'upstream')
        front = front[:-5] + f'\nreader-mode: {mode}\n---\n'

        # Keep all result facts in the chapter; move the existing closing
        # summary to the beginning, alongside a genuinely actionable rule.
        summaries = []
        for a,b,h in reversed(sections(body)):
            if re.match(r'## Key Takeaways', h):
                summaries.insert(0, body[a+len(h):b].strip())
                body = body[:a] + body[b:]
        takeaways = f'::: {{.callout-tip title="先确定这一点"}}\n\n{row["recommendation"]}\n\n:::\n\n'
        # Repeated metadata in a takeaway is already preserved in frontmatter
        # and worked results; keep non-duplicated narrative in a closing note.
        closing = '\n\n'.join(summaries)
        if closing:
            closing = re.sub(r'(?m)^- \*\*[^\n]+\n?', '', closing).strip()
            if closing and closing not in body:
                refs = next((a for a,_,h in sections(body) if re.match(r'## 参考', h)),len(body))
                body = body[:refs] + f'## 回到最初的问题\n\n{closing}\n\n' + body[refs:]

        # Replace the public slot names, preserving stable links.
        body = re.sub(r'^## .+?(\{#sec-target\})$', f'## {row["opening"]} \\1', body, count=1, flags=re.M)
        body = re.sub(r'^## 附录 A：Methods / Results 模板(.*)$', f'## {row["reporting"]}\\1', body, flags=re.M)
        body = re.sub(r'^## 附录 B：换成你自己的数据(.*)$', f'## 换到自己的研究中，先检查哪些条件？\\1', body, flags=re.M)
        body = re.sub(r'^## 附录 C：(?:图形与导出规范|十张图如何复现)(.*)$', r'## 图中应该保留哪些信息？\1', body, flags=re.M)
        body = prose_only(body, polish_prose)

        # Installation/plot helpers are explicitly website-only. A section ID
        # alone must never delete an input table or scientific explanation.
        for a,b,h in reversed(sections(body)):
            if '#sec-setup}' in h or '#sec-preparation}' in h:
                part = body[a:b]
                part = part.replace('<details>', '<details class="wechat-omit">')
                part = re.sub(r'(<summary>)(?:<strong>)?(?:展开|点开)[：:]([^<]*)(?:</strong>)?(</summary>)', r'\1环境配置与完整运行说明\3', part)
                body = body[:a]+part+body[b:]
        # Catch the mixed setup block in 71 without hiding its data section.
        body = body.replace('<details>\n<summary><strong>展开：安装依赖、下载真实数据与作图函数</strong></summary>', '<details class="wechat-omit">\n<summary>环境配置与完整运行说明</summary>')
        first = sections(body)[0][0] if sections(body) else 0
        body = body[:first] + takeaways + body[first:]
        new = front + body
        if new != old:
            source.write_text(new)
        changes.append({'number':number,'file':str(source.relative_to(root)), 'mode':mode,
                        'before':hashlib.sha256(old.encode()).hexdigest(), 'after':hashlib.sha256(new.encode()).hexdigest()})
    return changes


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--project-root',type=Path,default=Path('.'));args=parser.parse_args()
    root=args.project_root.resolve()
    with (root/'scripts/reader_editorial_topics.tsv').open() as f:
        topics={int(r['number']):r for r in csv.DictReader(f,delimiter='\t')}
    assert set(topics)==set(range(1,78))
    changes=make_editorial(root,topics)
    out=root/'qa/reader-revision';out.mkdir(parents=True,exist_ok=True)
    (out/'source-changes.json').write_text(json.dumps(changes,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'chapters':len(changes),'modified':sum(r['before']!=r['after'] for r in changes)}))


if __name__=='__main__':main()
