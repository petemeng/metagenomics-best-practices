# 宏基因组分析最佳实践 · Shotgun 系列（第二季）

> 仓库名：`metagenomics-best-practices`
> 定位：**以 shotgun metagenomics 为入口，复现高分论文中的物种、功能、MAG、菌株、病毒、代谢与跨队列分析，并把每类结果转化为可投稿证据。**
> 参照物：`sc-best-practices` 教程（网页书籍 + GitHub + 可复现环境）；姊妹季：16S 系列 `microbiome-best-practices`。
> 发布形式：**每篇独立成文**；网页书籍按技术体系归档，**公众号发布顺序先走完一个 read-based 发表闭环，再进 genome-resolved 深度路线**。

---

## 一、核心定位与原则

1. **上游标准化到可审计，下游围绕论文证据链深挖。**（内部简记：上游是过路、下游是肉。对外用前一句，避免被做组装/MAG/算法的人抓措辞。）上游压缩到能跑通即可，但要留全可复现记录；下游每篇锚定一篇真实高分文献 + 其公开代码，复现到发表级美观。
2. **复现驱动教学。** 不讲抽象方法，而是"复现某篇高分论文的 Figure X，作者代码在此，逐行拆解，再教你换成自己的数据"。
3. **区分度不只在图，在严谨性。** 跨队列验证、菌株级传播证据、MAG 质控标准（MIMAG）、绝对定量——这些研究设计的严谨性才是分水岭。
4. **诚实讲边界。** Kraken2 假阳性、MAG 完整度/污染/嵌合权衡、菌株"传播"的证据等级、代谢建模是潜力非实测、病毒—宿主关联是预测、相对丰度不等于绝对数量、纯生信到不了因果。
5. **算力与数据库是第一现实约束（宏基因组特有）。** 组装需大内存（单样本常 100+ GB RAM）、数据库数十至上百 GB、运行动辄数小时。**每篇上游章节开头必须标注硬件门槛（RAM/磁盘/核数/耗时）**，并给"没有服务器怎么办"（HPC/云/子集演示）。
6. **一条主干、三个分析层级（本季骨架，取代旧的"两条主线"）。** 从原始 reads 出发分三层，产出与结论边界各不同，第 01–02 篇必须画清地图，读者才不会把不同层级的表混用：
   ```
   原始 reads
   ├── 1. Reference-based read profiling（免组装谱）
   │      MetaPhlAn / mOTUs（物种）· Kraken2/Bracken · HUMAnN（功能）
   └── 2. Assembly-based（组装）
          ├── 2A. Gene / contig-centric（基因/contig 中心）
          │      非冗余基因目录 · eggNOG/KO/GO · CAZymes · ARG/VF · BGC · 病毒/质粒 contig
          └── 2B. Genome-resolved（基因组解析）
                 MAG · 系统基因组 · 菌株/SNV · 泛基因组 · 传播 · 基因组尺度代谢模型
   ```
   **关键纠正：eggNOG-mapper / dbCAN / CARD-RGI / antiSMASH 属于 2A（输入是预测蛋白/基因目录/contig/MAG），不属于 HUMAnN 所在的 read-based 功能层。**
7. **工具分工：命令行（bash/conda）跑上游与重计算，R 收尾下游统计与出版级作图。**

---

## 二、数据集策略（沿用 16S 季，针对宏基因组调整）

- **不强求单一数据集贯穿全季。** 每章用其锚定论文自己公开的数据 + 代码复现。
- **下游统计（人体）自足弹药库：`curatedMetagenomicData`（Bioconductor）**——大量已发表人群队列的标准化谱表 + 丰富临床元数据直接做成 R 对象，装包即得、无需重跑上游。相当于 16S 季 microeco、转录组季 airway、代谢组季 sacurine。
  - **版本溯源（必须讲清，否则是可复现性层面的错误陈述）**：`curatedMetagenomicData` 的标准化分类谱来自 **MetaPhlAn3**、功能潜力来自 **HUMAnN3**，**不是本教程上游所用的 MetaPhlAn4**。绝不能让读者以为"这些数据就是用本教程 MetaPhlAn4 流程生成的"。每篇用公开数据的章节须附"数据谱系卡片"（见 §六）。
- **环境方向公共数据入口**：`curatedMetagenomicData` 偏人体，环境向单列一篇讲 MGnify / GEM / GTDB / 环境 MAG catalog（第 73 篇）。
- **上游重计算章节无法"装包即得"自足**：用公开小数据子集（已发表研究少量样本或官方教程数据），`data/download.sh` 从 SRA/ENA/Zenodo 拉；关键中间产物（谱表、bins、MAG、CheckM2/GUNC 结果）固化进小数据区，让读者从任意"下游半段"单独开始。
- **数据入库规则：原始 FASTQ 与大数据库不进 Git**（放外部源脚本拉取，写清 release + checksum）；小数据、示例谱表、示例 MAG、中间结果表可进 Git。

---

## 三、教学功能检查（不作为固定公开小标题）

1. **实际问题与目标证据**（展示目标 Figure 的结构与论文链接；许可允许时展示原图，否则重绘分析示意或用自己的复现图——公众号转载顶刊原图需单独考虑图形许可）
2. **原理与关键结论**（先给结论，长推导折叠进可跳过的"深入"段落）
3. **准备工作**（可折叠，内联安装 + 作图函数；上游篇标注硬件门槛与数据库下载）
4. **可复制代码**（read-based/下游从谱表或 curatedMetagenomicData 读起；标注环境/耗时/内存/磁盘；术语首次中英对照）
5. **方法选择与局限**（先忠实复现作者做法，再用当前标准审计：这版方法/参数/验证放到今天，哪里已过时、有何局限、该换成更稳的替代或补什么验证——把"作者怎么做"升级为"今天该怎么做"。本系列区别于普通复现教程的核心一步）
6. **出版级美化**（色盲友好配色、字体、显著性标注；PDF/SVG 矢量输出，TIFF/PNG 按期刊要求 300–600 ppi）
7. **常见坑方框**（审稿人高频攻击点）
8. **这段 Methods 怎么写**（含工具版本 + 数据库 release）
9. **换成你自己的数据怎么做**

---

## 四、方向标注（多标签，取代单一图标）

一个章节可带多个标签，更贴合真实使用场景：

```
受众：[通用] [人体] [环境]
层级：[Read-based] [Gene-centric] [MAG]
读长：[短读长] [长读长]
算力：[低算力] [高算力]
环节：[统计] [组装] [发表]
```

示例：
```
CAZymes            [通用] [Gene-centric] [MAG] [功能]
菌株传播            [人体] [MAG] [菌株] [高算力]
C/N/S 循环          [环境] [MAG] [功能] [代谢]
跨队列验证          [人体] [Read-based] [统计] [发表]
```

---

## 五、全套篇目（约 77 篇，开放式，后续持续更新）

> 网页书籍按下面技术体系归档；**公众号发布顺序**：第 0–6 部分先走完 → 完成"read-based 人体疾病论文"发表闭环 → 再进第 7 部分起的 genome-resolved 深度路线。

### 第 0 部分 · 导论
| # | 标题 | 标签 |
|---|------|------|
| 01 | 开篇：宏基因组能回答什么、回答不了什么（vs 16S / 宏转录 / 其他组学） | 通用 |
| 02 | 一条主干、三个分析层级：read profiling / gene-contig-centric / genome-resolved 如何取舍 | 通用 |

### 第 1 部分 · 实验设计与样本（GPT: 深度与样本量拆开、补建库与绝对定量）
| # | 标题 | 标签 |
|---|------|------|
| 03 | 研究设计、样本量与统计效力（病例对照/配对/纵向/多中心/发现验证队列/功效模拟/批次混杂） | 通用·统计 |
| 04 | 宏基因组到底要测多深：决策框架（科学问题→目标产物→群落复杂度→宿主 DNA 比例→深度）+ 饱和曲线 | 通用 |
| 05 | **从样本到文库：DNA 提取、建库偏差与定量宏基因组**（提取/裂解强度/HMW DNA/PCR-free vs PCR/insert size/GC bias/mock/spike-in；**相对丰度≠绝对数量，shotgun 同样是组成型数据**；qPCR/流式/每克细胞数/绝对 taxon 与 resistome 丰度） | 通用 |
| 06 | 宿主去除与低生物量：湿实验宿主耗竭 + 计算去宿主的取舍 | 通用 |
| 07 | 污染与对照：blank、kitome、index hopping、decontam | 通用 |
| 08 | 认识数据与平台：FASTQ、短读 vs 长读（PacBio HiFi / Nanopore）对全流程的影响 | 通用 |

### 第 2 部分 · 环境搭建与算力
| # | 标题 | 标签 |
|---|------|------|
| 09 | WSL2 + conda/mamba（Windows 用户主线） | 通用 |
| 10 | **算力现实**：内存/磁盘/核数怎么估、HPC/云、**Apptainer/Singularity**（HPC 常禁 Docker daemon）、SLURM/job array/scratch/失败续跑 | 通用·高算力 |
| 11 | 安装 bioBakery（MetaPhlAn/HUMAnN）+ 组装分箱工具链 + 数据库落地（release + checksum） | 通用 |
| 12 | 安装 R + `curatedMetagenomicData` + 关键包生态 | 通用·统计 |

### 第 3 部分 · 上游 QC 与预处理
| # | 标题 | 标签 |
|---|------|------|
| 13 | 读长质控与修剪：FastQC/MultiQC + fastp | 通用 |
| 14 | **去宿主、复杂度评估与重复 reads 的识别与处理**（不把去重写成默认必做：PCR/optical duplicate vs 高丰度真实重复 vs 低复杂度；何时只报 duplicate rate、何时去重、去重是否改变丰度、UMI 数据的区别） | 通用 |

---
**══ 第一闭环：Read-based 快速发表路线（公众号先发到这里，读者即可完成一篇人体疾病宏基因组分析）══**

### 第 4 部分 · Reference-based read profiling（物种）
| # | 标题 | 标签 |
|---|------|------|
| 15 | **MetaPhlAn4：marker 基因、SGB 与物种级组成分析**（已知/未知 SGB、数据库 release、未分类 reads、marker coverage、低丰度阈值、检出 vs 定量；菌株留给 StrainPhlAn；版本+数据库日期必须写死） | 通用·Read-based |
| 16 | Kraken2 + Bracken：k-mer 分类 + 丰度重估 | 通用·Read-based |
| 17 | 数据库与置信阈值：Kraken2 库选择、confidence、假阳性控制 | 通用·Read-based |
| 18 | 谱工具横评：MetaPhlAn vs Kraken2/Bracken vs mOTUs，结果为何不同 | 通用·Read-based |

### 第 5 部分 · Reference-based read profiling（功能）
| # | 标题 | 标签 |
|---|------|------|
| 19 | HUMAnN3：gene family / pathway abundance / pathway coverage、stratified vs unstratified | 通用·Read-based |
| 20 | 功能谱归一化与解读：CPM/relative/coverage 区别、物种对通路的贡献、unmapped/unintegrated/unclassified | 通用·Read-based |

### 第 6 部分 · 下游统计与跨队列（read-based，curatedMetagenomicData 自足）
| # | 标题 | 标签 |
|---|------|------|
| 21 | **宏基因组表格到底是什么**：counts/coverage/CPM/relative abundance/genome equivalents 与组成型问题（哪种可直接求和、哪种需基因长度/深度校正、0 是未检出还是不存在、pathway coverage vs abundance、stratified 怎么处理） | 通用·统计 |
| 22 | 多样性：Alpha / Beta（物种谱、MAG 表、基因谱三种粒度） | 通用·统计 |
| 23 | 排序与检验：PCoA / CAP + PERMANOVA + betadisper | 通用·统计 |
| 24 | 差异分析：**MaAsLin3** / ANCOM-BC2 / ALDEx2（MaAsLin3 区分丰度关联与检出率关联；物种谱 + 功能谱都做；复现原论文→当前升级→敏感性分析→结论稳健性） | 通用·统计 |
| 25 | 组成与多分类级可视化、核心微生物组 | 通用·统计 |
| 26 | `curatedMetagenomicData` 深用 + **数据谱系卡片**（profiler/版本/数据库 release/feature 类型/单位/归一化必须标注；MetaPhlAn3/HUMAnN3 来源） | 人体·统计 |
| 27 | 随机森林 / 梯度提升分类 + ROC + 重要性 + 防过拟合规范 | 通用·统计 |
| 28 | **跨队列 meta 分析 + 外部验证**（SIAMCAT / 批次校正 / leave-one-dataset-out，疾病 signature 可重复性） | 人体·统计·发表 |
| 29 | 共现网络 + 网络进阶（鲁棒性 / keystone / Zi-Pi） | 环境·统计 |

> **▲ 第一闭环完成**：读者到此可独立完成一篇"人体疾病宏基因组 read-based"论文分析（物种+功能+差异+ML+跨队列验证+成套图）。

---
**══ 第二闭环：Genome-resolved 深度路线 ══**

### 第 7 部分 · 组装（补长读长与 hybrid）
| # | 标题 | 标签 |
|---|------|------|
| 30 | 短读长组装：MEGAHIT vs metaSPAdes、单样本 vs 混合组装（co-assembly） | 通用·组装·高算力 |
| 31 | **长读长组装**：Nanopore/PacBio HiFi、HMW DNA/basecalling、metaFlye/hifiasm-meta/metaMDBG、polishing、环状 contig、长读特有错误、methylation（长读改善重复区/rRNA operon/BGC/CRISPR-Cas/ARG-MGE 连接/完整基因组） | 通用·长读长·组装·高算力 |
| 32 | **Hybrid assembly 与 polishing**：短读准确性 + 长读连续性、strain mixture、何时混合 vs 直接 HiFi | 通用·长读长·组装 |
| 33 | 组装质控：QUAST、N50、contig 统计、什么样的组装可用 | 通用·组装 |

### 第 8 部分 · Gene / contig-centric 分析（GPT: 功能注释归到此层）
| # | 标题 | 标签 |
|---|------|------|
| 34 | 非冗余基因目录：prodigal 预测 + MMseqs2/CD-HIT 聚类 | 通用·Gene-centric |
| 35 | 基因丰度谱：reads 回比基因目录、按功能汇总（RPKM/TPM/CPM） | 通用·Gene-centric |
| 36 | Gene-centric 功能注释：eggNOG-mapper（KO/COG/GO）、功能暗物质与未注释比例 | 通用·Gene-centric·功能 |
| 37 | CAZymes 碳水活性酶（dbCAN） | 通用·Gene-centric·MAG·功能 |
| 38 | **Resistome 专题**：CARD-RGI/deepARG，read vs assembly，identity/coverage/alignment length，ARG abundance，抗性基因存在≠表型耐药，宿主归属，ARG 是否位于质粒/MGE | 通用·Gene-centric·MAG |
| 39 | **Virulome 专题**：VFDB / ABRicate | 通用·Gene-centric |
| 40 | **BGC 与天然产物**：antiSMASH/GECCO/DeepBGC、BGC 完整性/新颖性、BiG-SCAPE/BiG-SLiCE 聚类、BGC abundance、MAG-BGC 对应 | 通用·Gene-centric·MAG |

### 第 9 部分 · 分箱与 MAG（本季最成熟板块，补 GUNC/完整 MIMAG/人工校正）
| # | 标题 | 标签 |
|---|------|------|
| 41 | 回比与深度：reads 回比（bowtie2/bwa）、jgi 深度矩阵 | 通用·MAG·高算力 |
| 42 | 分箱器：**MetaBAT2（基线）/ SemiBin2 / TaxVAMB（深度学习）**；三个决策——单样本 vs 多样本、单一 vs ensemble、如何用 CheckM2+GUNC 选结果 | 通用·MAG·高算力 |
| 43 | Bin 精炼与整合：DAS Tool / binette | 通用·MAG |
| 44 | **MAG 质控：CheckM2 + GUNC + 完整 MIMAG + assembly graph**（CheckM2 擅长标志基因完整度/冗余污染，但远缘片段的非冗余污染/嵌合可能逃过它，需 GUNC clade separation；报告 completeness/contamination/strain heterogeneity/GUNC/contig 数/N50/genome size/rRNA/tRNA/coding density/GTDB/coverage/人工检查） | 通用·MAG·高算力 |
| 45 | 去冗余：dRep（跨样本 MAG 聚类到物种级基因组、ANI、representative genome） | 通用·MAG |
| 46 | MAG 分类：GTDB-Tk（GTDB 分类系统、release 写死） | 通用·MAG |
| 47 | 新物种与系统基因组：novel MAG 判定、GToTree/IQ-TREE 建基因组树 | 通用·MAG |
| 48 | MAG 丰度定量：CoverM（覆盖度/breadth/相对丰度） | 通用·MAG·统计 |
| 49 | **MAG 人工校正、命名、数据库提交与补充材料制作**（anvi'o/Bandage 检查、异常 coverage/GC、嵌合 contig、`Candidatus`/novel species 谨慎表述、NCBI/ENA 提交、BioProject/BioSample/WGS-MAG accession、MIMAG metadata、MAG 质量总表） | 通用·MAG·发表 |

### 第 10 部分 · 菌株与种内变异（补泛基因组/附属基因）
| # | 标题 | 标签 |
|---|------|------|
| 50 | inStrain：SNV、微多样性（microdiversity）、nucleotide diversity | 通用·菌株·高算力 |
| 51 | StrainPhlAn：marker 基因菌株系统发育 | 人体·菌株 |
| 52 | **种内功能差异：泛基因组、附属基因与菌株功能型**（PanPhlAn/MIDAS2/anvi'o pangenomics；core vs accessory、gene presence/absence、strain-specific genes、species phylogeny vs gene-content tree；同一物种≠同一 ARG/VF/BGC/质粒/phage defense） | 通用·菌株·MAG |
| 53 | 菌株传播与共享：母婴垂直传播、FMT 定植（engraftment）、传播的证据等级（种水平相似≠传播，需菌株 SNV + 排除共同来源混杂） | 人体·菌株 |

### 第 11 部分 · 病毒 / 质粒 / 真核（补 MIUViG、定量、实验类型差异）
| # | 标题 | 标签 |
|---|------|------|
| 54 | **病毒发现与质量标准**：geNomad/VirSorter2、CheckV、vOTU 聚类、**MIUViG**；**total metagenome vs virus-enriched virome（filtration/nuclease/WGA，两种实验不能直接比）** | 通用·病毒 |
| 55 | 病毒分类与丰度：vOTU abundance、read mapping threshold、breadth of coverage、viral taxonomy、temperate/lytic、prophage vs 游离病毒、跨样本共享 | 通用·病毒 |
| 56 | **病毒—宿主关联的证据等级**：同一 MAG 上 prophage > CRISPR spacer > 长读/Hi-C 连接 > k-mer/组成预测 > 共丰度（各层说明假阳性；iPHoP 输出是预测非证实） | 通用·病毒 |
| 57 | 质粒与可移动元件：geNomad 质粒识别、MGE、ARG 的可移动性 | 通用·Gene-centric |
| 58 | 真核微生物：EukDetect / EukRep（真菌、原生生物） | 环境·人体 |

### 第 12 部分 · 代谢建模（MAG 功能重建，属 genome-resolved）
| # | 标题 | 标签 |
|---|------|------|
| 59 | 基因组解析代谢：DRAM / METABOLIC（MAG 代谢潜力注释） | 通用·MAG·代谢 |
| 60 | 基因组尺度代谢模型（GEM）：gapseq / CarveMe 从 MAG 建模（依赖 gap-filling，结论谨慎） | 通用·MAG·代谢 |
| 61 | 群落代谢与流量：MICOM / SMETANA（互养、竞争、代谢流预测——潜力非实测通量） | 人体·代谢 |
| 62 | 元素循环与功能类群（环境向）：C/N/S 循环通路盘点 | 环境·MAG·代谢 |

### 第 13 部分 · 多组学整合（补宏蛋白组，完善证据层级）
| # | 标题 | 标签 |
|---|------|------|
| 63 | 与代谢组整合：Procrustes / Mantel / HAllA / sPLS-DIABLO | 人体·统计 |
| 64 | 宏转录组联动：功能"有没有"vs"表达没表达"、DNA/RNA 比 | 通用 |
| 65 | **宏蛋白组：从功能潜力到真实蛋白表达**（metagenome-informed 蛋白库、PSM、shared peptide、protein inference、微生物 vs 宿主蛋白、taxon-specific protein、pathway 汇总；DNA→RNA→protein→metabolite 四层对齐） | 人体 |
| 66 | MAG × 代谢物 / 表型：把菌株基因组与功能读出对接 | 人体·MAG |

> **四层证据链**：DNA 能不能做 → RNA 是否转录 → Protein 是否表达 → Metabolite 是否形成可检测生化结果。

### 第 14 部分 · 纵向与因果（与 16S 季因果簇衔接，可复用）
| # | 标题 | 标签 |
|---|------|------|
| 67 | 纵向：线性混合模型 / 轨迹 / 定植稳定性 | 人体·统计 |
| 68 | 微生物组与生存结局：Cox / Kaplan–Meier / 时间依赖 ROC（cutoff 数据泄漏） | 人体·统计 |
| 69 | 中介分析：暴露 → 微生物（含菌株/功能）→ 结局 | 人体·统计 |
| 70 | 孟德尔随机化：工具变量、水平多效性、敏感性分析 | 人体·统计 |
| 71 | 结构方程模型：环境因子 → 微生物 → 表型（piecewiseSEM / PLS-PM） | 环境·统计 |
| 72 | 因果证据阶梯：关联→纵向→中介/MR→干预→培养→FMT/无菌鼠→分子机制 | 通用·发表 |

### 第 15 部分 · 公共数据资源
| # | 标题 | 标签 |
|---|------|------|
| 73 | **公共宏基因组资源：MGnify / GEM / GTDB / 环境 MAG catalog**（补 curatedMetagenomicData 的人体偏向，环境方向的公开数据入口） | 环境·统计 |

### 第 16 部分 · 工程化与发表（GPT: 补真正面向发表的收尾）
| # | 标题 | 标签 |
|---|------|------|
| 74 | 用 nf-core/mag 等 Nextflow/Snakemake 流程一键复现全上游（**固定 `-r <release>`、Apptainer profile、`-params-file`、SLURM、数据库 release+checksum**；nf-core/funcscan 组织 ARG/AMP/BGC） | 通用·组装·发表·高算力 |
| 75 | 一篇宏基因组论文的主图和补图如何组织（研究设计/差异物种功能/MAG-菌株-病毒发现/机制或代谢模型/外部或实验验证/补充 QC+MAG 表+数据库版本+敏感性分析） | 通用·发表 |
| 76 | **报告标准：STORMS（人体）/ STREAMS（环境与非人宿主）/ MIMAG / MIUViG** | 通用·发表 |
| 77 | 数据、代码与基因组提交：SRA/ENA raw reads、BioProject/BioSample、metadata、MAG/vOTU/gene catalog/annotation、Zenodo release、GitHub tag、container digest、database manifest、Methods/Data availability/Code availability | 通用·发表 |

> **重点与边界说明**
> - **三个分析层级（§一.6）是本季骨架。** eggNOG/dbCAN/ARG/BGC 属 gene-contig-centric，不与 HUMAnN 同层——这是相对旧版大纲最重要的结构纠正。
> - **MAG（第 9 部分）** 是本季最大的"16S 做不到"卖点。质控从 CheckM2 升级为 **CheckM2 + GUNC + 完整 MIMAG**（含 5S/16S/23S rRNA 与足够 tRNA），区分 near-complete / MIMAG high-quality / medium-quality / circular complete，别把 `>90%/<5%` 一律写成"高质量完整基因组"。
> - **菌株（第 10 部分）** 补泛基因组/附属基因，形成"物种丰度→菌株系统发育→SNV→附属基因→功能差异"链条，是真正超越 16S 的核心。
> - **病毒（第 11 部分）** 拆成发现+质量标准（MIUViG、total vs enriched virome）、分类+定量、宿主关联证据等级三篇。
> - **绝对定量（第 05 篇）** 补相对丰度陷阱——shotgun 同样是组成型数据。
> - **长读长/hybrid（第 31–32 篇）** 补齐"导论提了却没教"的断裂，可能是本季最有吸引力的内容之一。
> - **curatedMetagenomicData 版本溯源**：MetaPhlAn3/HUMAnN3 来源，非本教程 MetaPhlAn4，每篇附数据谱系卡片。
>
> **额外数据要求（影响"单篇自足"）**：第 68（生存）需临床随访、第 70（MR）需外部 GWAS、第 71（SEM）需实测环境+表型，开头须写清"你还需要什么数据"，MR 用公开示例演示。

---

## 六、GitHub + 网页版技术方案

### 选型：Quarto（同 16S 季）
- 原生跑 R / Python / bash 代码块，渲染书籍型网站，内建搜索/暗色模式/多语言/代码复制/折叠。

### 推荐仓库结构
```
metagenomics-best-practices/
├── _quarto.yml            # 站点配置：约 77 篇章节树、主题、导航
├── index.qmd              # 首页 = 第 01 篇
├── chapters/              # 每篇一个 .qmd，文件名带编号 + 英文 slug
├── R/theme_pub.R          # 共享作图函数（与四季同款，供整仓库 source）
├── env/
│   ├── biobakery.yml      # bioBakery/read-based conda 环境（锁版本）
│   ├── assembly.yml       # 组装/分箱/MAG conda 环境（锁版本）
│   └── renv.lock          # R 包锁版本
├── db/
│   └── download_db.sh     # 大数据库（不进 git，写清 release+checksum）：MetaPhlAn/Kraken2/GTDB/CheckM2/GUNC/CARD/eggNOG…
├── data/
│   ├── download.sh        # 原始/示例数据不进 git，从 SRA/ENA/Zenodo 拉
│   └── small/             # 小数据、示例谱表、示例 MAG、中间结果可进 git
├── figures/
├── Dockerfile             # 一键可复现环境（镜像体积与数据库分离）
├── .github/workflows/publish.yml   # push 后自动 render + 部署
├── LICENSE                # 内容 CC-BY-4.0，代码 MIT
└── README.md
```

### 数据谱系卡片（每篇用公开数据的章节固定附）
```
Study / Accession / Raw reads available /
Profile source / Profiler / Profiler version / Database release /
Feature type / Unit / Normalization /
Original publication / Reprocessed or author-provided
```

### 关键决策（沿用四季，宏基因组特有项加粗）
- **服务器实跑出图 + 图直接进推文**。
- **执行分层（本季尤其重要）**：下游统计章节（第 6、13、14 部分，基于 curatedMetagenomicData/谱表）`eval: true` 真跑 + `freeze: auto`；上游重计算章节（组装/分箱/MAG/病毒/HUMAnN）`eval: false`，只跑一次固化中间产物。
- **`eval:false` 不能只留截图**：必须保存完整命令、`params.yml`、软件版本、数据库版本、运行日志、MultiQC、峰值内存、CPU/wall time、输入与关键输出 checksum、pipeline commit/tag。截图只作教学展示，不作可复现记录。
- **固定工作流版本**：不写 `nextflow run nf-core/mag`，而写 `nextflow run nf-core/mag -r <release> -profile apptainer -params-file params.yml`；具体 release 由 Codex 在服务器核定并锁定（不在文档硬编，因版本随时更新）。
- **图内文字一律英文**（无头服务器无中文字体会出豆腐块；且便于投稿）。
- **确定性**：随机步骤（分箱、抽样、ML、置换）`set.seed`/固定 `--seed`。
- **作图规范集中在 `R/theme_pub.R`**（单篇正文内联作图函数保证自足）。
- **图导出**：PDF/SVG 矢量；TIFF/PNG 按期刊要求 300–600 ppi（矢量图不用 dpi 描述）。
- **数据库与数据分层**；**算力标注制度**（每篇上游章节标 RAM/磁盘/核数/耗时）。
- **部署**：GitHub Pages + Actions；国内镜像 Netlify/Gitee Pages。
- **公众号 ↔ 网页落差**：.qmd 写作 → mdnice/Md2All 转公众号；网页是代码权威源。
- **双语预留**：文件名/URL 用英文/拼音。

---

## 七、写作规范（六条铁律 + 审计升级，与四季一致）

1. **每章从谱表/表格或 curatedMetagenomicData 读起**（下游），上游从明确示例数据/中间产物读起；不隐性依赖别篇。
2. **按类型自足**：技术篇在首次读取前提供输入下载、完整代码和环境；导读/概念/证据篇直接展示案例与判断，不强加安装、作图函数和 Methods 模板。公众号只略去显式标记的通用环境代码，不删除科学数据入口。
3. **用真实、可引用数据**：下游用 curatedMetagenomicData（附数据谱系卡片），上游用公开研究子集 + 注明出处；至少 `head()` 展示结构；章末给参考。
4. **图内全英文**。
5. **正文不写元说明**：不写"本篇可独立跑通""这体现全系列约定"之类；不写"作者代码通常长这样"这种评论别人代码的段落。
6. **确定性 set.seed**。
7. **方法比较与限制（作者自检，不固定为第 5 段）**：忠实复现后用当前标准审计作者做法的局限并给升级——本系列灵魂，也是人工审阅的重点段落。

---

## 八、待办 / 下一步

- [ ] 生成/更新本季 `AGENTS.md`（三层级归类、两闭环发布顺序、CheckM2+GUNC、MaAsLin3、绝对定量、数据谱系卡片、eval:false 完整记录、固定 workflow release、Apptainer）
- [ ] 搭 Quarto 脚手架（`_quarto.yml` 含约 77 篇章节树、双 conda 环境 + renv、CI、README）
- [ ] 写 2 篇样板定两大卖点范式：
      **① 第 44 篇（MAG 质控：CheckM2+GUNC+MIMAG，如何制作审稿人认可的 MAG 质量总表）** — genome-resolved 上游范式
      **② 跨队列验证（用 curatedMetagenomicData 复现疾病 signature + leave-one-dataset-out）** — 人体 read-based 发表闭环范式（比普通 Alpha/Beta 更能体现本季区分度）
- [ ] 为每章确定锚定文献 + 代码出处
- [ ] 关键弹药库：bioBakery（MetaPhlAn4/HUMAnN3/StrainPhlAn/MaAsLin3/PanPhlAn）、MIDAS2、nf-core/mag、nf-core/funcscan、GTDB-Tk、CheckM2、GUNC、dRep、inStrain、SemiBin2/TaxVAMB、DAS Tool/binette、metaFlye/hifiasm-meta/metaMDBG、geNomad/CheckV/iPHoP、antiSMASH/BiG-SCAPE、DRAM/METABOLIC、gapseq/CarveMe/MICOM、SIAMCAT、curatedMetagenomicData、MGnify
