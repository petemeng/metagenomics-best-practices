#!/usr/bin/env python3
"""Explicit topic-based placements, not a lexical guess of scientific meaning."""
import re
from pathlib import Path
from revise_reader_narrative import split, sections, prose_only
root=Path(__file__).resolve().parents[1]
placements={
3:{'03-cohort-balance':'把总样本量拆回队列内比较','03-covariate-completeness':'审计协变量能否进入预定模型','03-power-sensitivity':'用两个 alpha 情景计算最小可检测差异'},
4:{'04-crc-library-depth':'读取 768 个 CRC 样本的 reported library size','04-nonpareil-saturation':'读取真实 Lake Lanier coverage 曲线','04-endpoint-depth-anchors':'对照同一研究中的三个 endpoint anchors'},
5:{'05-extraction-bias':'读取并验证 27 个真实 mock profiles','05-protocol-bias-range':'从 relative efficiency 计算成对偏差','05-syndna-quantification':'审计 synDNA 是否恢复已知 cell composition'},
19:{'19-read-flow':'完整运行 nucleotide + translated workflow','19-gene-family-stratification':'先在 native RPK regroup，再分别派生 CPM 与 relative abundance','19-pathway-contributions':'从表读起，不从上一章的内存对象起步','19-abundance-coverage':'从表读起，不从上一章的内存对象起步'},
20:{'20-normalization-denominators':'用真实 HUMAnN 3.9 跑四个代表性分支','20-special-feature-budget':'从审计表读出 denominator，而不是凭感觉解释','20-pathway-contributions':'展开全部 24 个分支','20-prevalence-zero-sensitivity':'在真实队列中比较 prevalence definitions'},
21:{'21-table-unit-map':'用各工具的原生列名读取，不先改成统一的 abundance','21-denominator-closure':'审计 counts、model estimates 与 closure','21-genome-equivalent-calibration':'从同一 sequence universe 计算 genome equivalents 与 RPKG','21-zero-strata-semantics':'分层关系按 output type 审计'},
22:{'22-resolution-boundaries':'读入三种谱表并对齐 metadata','22-alpha-hill-numbers':'统一计算 Hill q=0/1/2','22-beta-ordination':'保留 raw negative eigenvalues，再画 PCoA','22-sensitivity-recruitment':'过滤、零值与 recruitment sensitivity'},
23:{'23-design-permutation-space':'生成 9,999 个 region-restricted 唯一置换','23-pcoa-cap':'PCoA 与 partial CAP 共用同一 Bray-Curtis 几何','23-permanova-dispersion':'PERMDISP 使用 spatial median 和同一置换矩阵','23-pairwise-sensitivity':'不用“更显著”选 sensitivity branch'},
24:{'24-zeller-marker-redraw':'读表并锁定 110 名独立受试者','24-maaslin3-two-part':'MaAsLin3 两部分主模型','24-functional-associations':'MaAsLin3 两部分主模型','24-method-concordance':'读取冻结结果并只汇报可报告行'},
26:{'26-resource-release-audit':'用 `dryrun` 发现资源，再锁定自动选择','26-metadata-completeness':'检查 metadata 完整性与跨 study 主键','26-query-attrition':'把观察单位写进样本查询','26-lineage-compatibility':'在 `mergeData()` 前执行五项合同'},
52:{'52-pangenome-prevalence':'运行主阈值','52-gene-content-pcoa':'由 presence/absence 矩阵计算 prevalence、Jaccard 与 PCoA','52-gene-content-dendrogram':'由 presence/absence 矩阵计算 prevalence、Jaccard 与 PCoA'},
53:{'53-species-sharing-negative-control':'先用 species table 做失败得很有价值的负对照','53-mother-infant-strain-evidence':'把论文报告的 marker-SNV 证据单独登记','53-relatedness-classifier':'审计 SameStr 的相关/无关负对照与 FMT 三角设计'},
66:{'66-metabolite-overlap':'先画证据边界，再解释差异','66-mag-quality-landscape':'以固定 seed 重建 MAG，并保留完整质量证据'},
}
count=0
for p in sorted((root/'chapters').glob('*.qmd')):
    n=int(p.name[:2]);front,body=split(p.read_text())
    for stem,target in placements.get(n,{}).items():
        m=re.search(r'^!\[[^\n]*\([^\n]*'+re.escape(stem)+r'\.[^\n]*\n',body,re.M)
        if not m:raise RuntimeError(stem)
        block=m[0];body=body[:m.start()]+body[m.end():]
        candidates=[(a,b,h) for a,b,h in sections(body) if h.startswith('## '+target)]
        if len(candidates)!=1:raise RuntimeError((n,target))
        _,end,_=candidates[0];body=body[:end].rstrip()+'\n\n'+block+'\n'+body[end:];count+=1
    # The data entry belongs after motivation, not before the first question.
    ss=sections(body)
    if len(ss)>1 and '#sec-reader-inputs}' in ss[0][2] and '#sec-target}' in ss[1][2]:
        a,b,_=ss[0];c,d,_=ss[1]
        if '`r ' not in body[c:d]:body=body[:a]+body[c:d]+body[a:b]+body[d:]
    def headings(text):
        def one(m):
            s=m[0]
            for a,b in [('审计','检查'),('冻结','保存'),('固化','保存'),('验收','结果核对'),('门禁','检查'),('账本','记录表'),('合同','约定')]:s=s.replace(a,b)
            s=s.replace('一次性 builder 与日常离线 QA','完整分析与结果核对')
            s=s.replace('从表读起，不从上一章的内存对象起步','读入功能表，核对丰度与分层关系')
            s=s.replace('自己的 FASTQ 用当前 SameStr 工作流完整生成 SNV profiles','新增样本需要哪些独立的菌株证据？')
            return s
        return re.sub(r'^#{2,6} .+$',one,text,flags=re.M)
    p.write_text(front+prose_only(body,headings))
print(f'relocated {count} result figures')
