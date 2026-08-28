# metagenomics-best-practices

面向“懂生物、不懂生信”科研读者的 77 篇中文 shotgun metagenomics
教程。网页版使用 Quarto Book；微信公众号文章从同一份 QMD、真实公开数据、
执行结果和原创重绘图派生。

## 发布规则

- `tutorial.yaml` 是内容、执行和发布的唯一契约。
- 网页版是完整代码权威源；公众号版重讲解但保留数据入口、关键参数和结果边界。
- 默认引用论文并从公开数据重算、重绘，不直接复制无明确许可的论文原图。
- 下游章节必须真实执行；上游重型章节必须保留一次性运行的完整命令、参数、
  版本、日志、资源、checksum 和固化产物。
- 只有 `qa_report.json.status == "passed"` 才能生成发布包。
- 公众号发布包先在本地完成严格审计；经明确授权后只创建草稿，不自动发布或群发。
- GitHub Pages deploy job 默认关闭，只有人工设置
  `ENABLE_PRODUCTION_DEPLOY=true` 后才允许部署。
- 公开仓库会把纯文本运行日志中的本地工作区根路径归一为 `/workspace`；
  原始数据、参数和结果数值不改。可由脚本重建的第三方大报告与上游 changelog
  不纳入 Git。

## 当前状态

- 新版大纲固定为 77 篇和 17 个技术部分。
- `tutorial.yaml`、`_quarto.yml` 与 77 个 QMD 路径已经建立。
- 第 28 篇旧样板已重做为真正 N−1 队列训练、留一完整队列外测。
- 第 29 篇已完成 spring-level 条件关联网络、边稳定性、Zi-Pi、拓扑零模型与删除鲁棒性审计。
- 第 30 篇已完成两个真实不均一 DNA mocks 的 MEGAHIT/metaSPAdes single/co 六分支组装、八次 read-back mapping 与资源审计。
- 第 31 篇已完成同一 71-strain MOCK1 DNA 的完整历史 ONT R9/PacBio HiFi 四分支长读组装、base-level read-back、circular-junction 与资源审计。
- 第 32 篇已完成同一 MOCK1 DNA 的 short-only、两条 short-read-first hybrid、两条 long-only 与两条 Polypolish 分支，并以 71 references 审计 recovery、consensus error、低丰度回收和资源。
- 第 33 篇已统一审计第 30–32 篇的 15 套非重复真实 assemblies，以 71/87 exact references、17 条评价记录和两个确定性正控拆开 N50、NA50、recovery、misassembly 与任务可用性。
- 第 34 篇已从第 30 篇六套真实 assemblies 预测 441,407 个 ORFs，按 individual、co-assembly 与 mix 三种策略建立非冗余基因目录，并以 87 个 exact mock genomes、MMseqs2/CD-HIT 与阈值敏感性审计代表序列、成员谱系和真值边界。
- 第 35 篇已将 MOCK1/MOCK2 共 7,999,482 条 clean reads 回比 93,782-gene 主目录，分开 83.6705%/83.3977% raw mapping 与 2,784,234/2,777,443 条主 assigned reads；89,339 个基因获得 UniRef90 best hit，15,392 个基因连到 MetaCyc reaction，并以 CPM/RPKM/TPM closure 和一对多守恒建立功能账本。
- 第 36 篇已用 eggNOG-mapper 2.1.15、eggNOG 5.0.2 与 DIAMOND 2.0.15 全量注释同一 93,782-protein 目录：84,511 条获得 seed ortholog，50,800 条带 KO，non-electronic/all-evidence GO 分别覆盖 8,821/20,471 条；四态 gene/read ledger 和 COG/KO/GO fractional allocation 全部守恒。
- 第 37 篇已用 dbCAN 5.2.9 的 DIAMOND、family HMM 与 dbCAN-sub 实跑同一 93,782-protein 目录：3,899 条候选中 2,050 条达到两工具主共识、1,605 条为三工具共识；*B. thetaiotaomicron* 正控回收 117 个 CGC，其中 25 个带底物预测。
- 第 38 篇已用 RGI 6.0.8 与 CARD 4.0.1 实跑 catalog、co-assembly 和双阳性对照：主集合为 36 条 Perfect/Strict genes，6,763 条 Loose-only genes 只进敏感性分析；ARG sequence evidence 与耐药表型明确分开。
- 第 39 篇已用 ABRicate 1.4.0、BLAST+ 2.17.0 与 VFDB 2026-07-24 实跑 core/full 和双阈值：core 90/80、core 80/80、full 90/80 分别回收 93、123、184 条 catalog genes，并以双物种正控审计 context 与表型边界。
- 第 40 篇已用 antiSMASH 8.0.4 与 GECCO 0.10.3 实跑完整/切碎 *Salinispora*、*Nostoc* 和真实 co-assembly：后者得到 71/21 个 regions，71 个 antiSMASH regions 链接 1,499 个 catalog genes；fragmentation、MIBiG similarity 与化合物新颖性边界分别审计。
- 第 41 篇已把 MOCK1/MOCK2 的 1,999,853/1,999,888 对 clean reads 回比同一 18,354-contig MEGAHIT co-assembly：overall alignment 为 87.48%/87.25%，breadth 为 95.2965%/95.0964%；JGI 与 Samtools 深度使用完整同坐标账本。
- 第 42 篇已在同一 10,203-contig、74.93-Mb 坐标上实跑 MetaBAT2 single/multisample、SemiBin2、VAMB 与 TaxVAMB 五分支，再对 171 个 raw bins 独立运行 CheckM2 1.1.0 和 GUNC 1.1.0；truth 只用于事后 benchmark。
- 第 43 篇已从五套 checksum-locked partitions 实跑 DAS Tool 与 Binette，并对 53 个 refined candidates 独立重跑 CheckM2/GUNC；truth-blind 规则选择 23 个不重叠 Binette candidates，共 57,793,968 bp。
- 第 44 篇已对 23 个 selected MAGs 独立运行 CheckM2/GUNC/CheckM1 与完整 MIMAG 审计：17 个通过 pre-marker HQ 前门，最终仅 4 个 high-quality、19 个 medium-quality；58,070 个 FASTG nodes 以 exact-sequence SHA-256 映回 k141 坐标并解析 5,843 条 edges。
- 第 45–49 篇已完成 dRep 去冗余、GTDB-Tk R232 分类、新物种系统基因组、CoverM 定量以及 MAG 校正/命名/提交，均使用 checksum-locked 真实 MAG 与固化产物。
- 第 50–53 篇已完成 inStrain、StrainPhlAn、泛基因组/附属基因与传播证据等级，并保留同株阈值、输入可调用性和阴性结论边界。
- 第 54–58 篇已完成病毒发现/CheckV/vOTU、病毒分类丰度、病毒—宿主证据、质粒移动元件与真核微生物章节；候选、分类、宿主与表型证据分层报告。
- 第 59 篇已对 24 个真实 PRJEB52977 representative MAGs 加 4 个确定性截短对照实跑 DRAM 1.5.0 / METABOLIC-G 4.0：汇总 66,519 条 annotations，工具间 KO 中位 Jaccard 为 0.853，module agreement 为 98.4%；8,635/8,635 项离线检查通过。
- 第 60 篇已从 8 个真实 MAG 与 4 个确定性截短对照生成 gapseq 2.1.0 / CarveMe 1.6.6 的 48 个 SBML 模型，并完成 168 次 FBA 审计。24/24 个 gap-filled models 在声明的 construction medium 中可行，48/48 个 no-uptake audits 无正 biomass，4/4 个 exact controls 结构一致；403 个冻结 payload 的 676/676 项独立复算检查通过。
- 第 62 篇已用 500 个真实 shotgun metagenomes、56 个 hot springs 与 780 个 MAGs 审计 C/N/S 元素循环潜力。主分析保留 DiTing-style abundance index，同时要求完整 DNF marker gate；MAG carrier 仅按严格完整规则定义。温度、pH 与 community–MAG concordance 分别有 3、7、14 个 FDR < 0.05 关联，818 个 strict process–MAG calls 中只有 1 个 MAG 覆盖完整四步 denitrification chain；37 个冻结 payload 的 251/251 项独立复算检查通过。
- 第 01 篇使用 Wirbel CRC 公开元数据建立测量边界和跨队列研究设计。
- 第 02 篇用 Wirbel CRC 与 UHGG/UHGP 真实指标建立三层分析路线决策。
- 第 03 篇从同一 CRC metadata 审计 16 个设计臂、协变量缺失和效力规划边界。
- 第 04 篇联合 CRC reported library-size、Lake Lanier Nonpareil coverage
  与环境 AMR 深测序锚点，建立 endpoint-specific 深度停止规则。
- 第 05 篇用 Costea Phase III mock community、McLaren 偏差分解与
  synDNA 定量基准，建立提取/建库偏差、全过程 spike-in 和绝对定量工作流。
- 第 06 篇用 Marotz 唾液配对实验与 Longhi 痰液剂量系列，区分湿实验宿主
  耗竭、计算去宿主、靶标保留和 residual-human privacy audit。
- 第 07 篇用 Salter 真实 shotgun 稀释实验、`decontam` 官方 MUClite
  对照对象和 UDI/index-pair 证据，拆开流程污染、逐样本 burden 与 index hopping。
- 第 08 篇用同一 71-strain MOCK1 的 Illumina、historical ONT R9 与
  PacBio CCS 数据，拆开 FASTQ Q、mapped identity、read span 与 assembly outcome。
- 第 09 篇把 Windows/WSL2 控制面与 Linux 执行面分开审计，并在专用
  conda 环境中核对 71 个已解析包、4 个直接依赖和第 08 篇真实数据 lineage。
  Windows 命令按 Microsoft 官方文档给出，但未在这台原生 Linux 主机上冒充执行。
- 第 10 篇把 4 个 ENA FASTQ 的 10.945 GB 压缩输入盘点、三个确定性本地任务、
  checksum sentinel、Apptainer 1.5.2 容器合同和 SLURM array 模板接成资源闭环；
  本机没有执行 Apptainer 或 SLURM，因此相应证据保持 `NOT_RUN`。
- 第 11 篇实际安装 bioBakery 与 assembly/binning 两个 native-Linux 环境，
  冻结 554/239 个包并验证 16 个入口、3 个内置测试；CONCOCT 的 user-site/
  `pkg_resources` 故障及 HUMAnN 覆盖 Bowtie2/DIAMOND 入口的故障均保留并修复。
  第 11 篇当时的五库未下载边界仍保留；MetaPhlAn 与 HUMAnN 数据库随后分别由
  第 15、19 篇另立清单真实验证；CheckM2 与 GUNC 后续由第 42–44 篇以锁版
  数据库真实验证，GTDB-Tk 与 dRep 则由第 45–49 篇另立锁版证据验证。
- 第 12 篇固定 R 4.4.1 / Bioconductor 3.19 与 185-record `renv.lock`，
  精确核对 17 个核心包；一次性真实取回 curatedMetagenomicData 3.12.0 的
  ExperimentHub `EH7091`，冻结 298 features × 24 samples 的完整
  `TreeSummarizedExperiment`，并以 121/121 检查证明离线复用、对象方向、
  相对丰度单位和 phyloseq 转换守恒。
- 第 13 篇从 Meslier 71-strain MOCK1 的官方 Illumina run `ERR9765746`
  确定性读取前 100,000 个同步 pairs，实际运行 FastQC 0.12.1、fastp
  1.3.6 与 MultiQC 1.35。保守过滤保留 99,991 pairs，九对仅因最短
  50 bp 门槛被移除，adapter-trimmed reads 为零；离线验证 74/74 PASS，
  57 个冻结文件 checksum 全部通过。
- 第 14 篇用 Hostile 论文中的公开人源 WGS run `ERR194147` 与同一 MOCK1
  `ERR9765746` 建立分离的阳性/保留双控制，各取前 20,000 个同步 pairs。
  实际运行 Hostile 2.0.2、Bowtie2 2.5.4、Samtools 1.21 与完整
  `human-t2t-hla` 2023-07 索引；人源 pairs 删除 99.750%，微生物 pairs
  保留 100.000%。低复杂度过滤和序列精确去重仅作敏感性分支，主流程均关闭；
  离线验证 52/52 PASS，47 个冻结 payload checksum 全部通过。
- 第 15 篇从第 13 篇的 99,991 个 clean pairs 运行 MetaPhlAn 4.2.5、
  Bowtie2 2.5.5 与完整 `mpa_vJan26_CHOCOPhlAnSGB_202605` 数据库。199,929 条
  reads 进入 profiling，12,739 条直接 marker-map records 支持 31 个 species / 31 个
  known SGB；terminal-clade estimated reads 为 261,721，因此默认 unclassified
  被截断为 0.000%，但不解释为每条 read 都已分类。一次性初始化 105/105、
  离线验证 100/100、23 个冻结 payload checksum 全部通过。
- 第 16 篇复用同一 99,991 个 clean pairs，实际运行 Kraken2 2.17.1、
  Bracken package 3.1p1（wrapper CLI 3.0.1）与 checksum-locked
  Standard-8 2026-06-26 数据库。Kraken2 分类 86,373 个 paired fragments，
  13,618 个未分类；species/150 bp/threshold 10 的 Bracken 主模型保留
  65 个物种，最终新增 24,156 个 fragments 并报告 84,147 个 species-level
  estimates。threshold 0、100 bp 和 genus 分支分别审计 feature count、
  read-length 与 rank 敏感性；一次性初始化 86/86、离线验证 73/73、
  33 个冻结 payload checksum 全部通过。
- 第 17 篇在同一 2026-06-26 release 下实际比较 Standard-8、Standard-16
  与 PlusPF-8，使用 99,991-pair MOCK1 阳性真值和 20,000-pair 人类方法对照，
  完成 30 个 database × confidence × control 主分支及 6 个新增 hit-group
  分支。三库均覆盖 63/69 个 current truth species；操作性比较点
  Standard-16、confidence 0.10、hit groups 2 分类 65,882 个 fragments，
  恢复 55/63 个 reference-present species，人类 lineage-external burden 为
  0.5/10k。一次性初始化 404/404、离线验证 74/74、111 个冻结 checksum
  全部通过；该点明确不写成通用最优阈值。
- 第 18 篇在同一 99,991-pair MOCK1 输入上横评 MetaPhlAn 4.2.5、
  Kraken2 2.17.1 + Bracken 3.1p1 与 mOTUs 4.1.0。mOTUs database 4.1
  实际声明 GTDB R226；实跑获得 1,122 条 marker-aligned reads 和 673 个
  aligned inserts。显式 NCBI species—SGB—mOTU crosswalk 将 69 个 current
  truth species 收敛为 52 个三方严格一对一物种，覆盖 70.738% expected mass；
  默认分支分别恢复 31/69、53/63、29/64 个 reference-present species。
  一次性初始化 83/83、离线验证 101/101、28 个冻结 payload checksum 全部通过。
- 第 19 篇在同一 99,991-pair clean MOCK1 上真实运行 HUMAnN 3.9 的 nucleotide
  与 translated 两级检索。MetaPhlAn vJun23 prescreen 选择 32 个 profile species，
  HUMAnN alias mapping 对应 89 个 database names 和 37 个 ChocoPhlAn files；
  prescreen 未被绕过。199,982 条 reads 分为 63,968 nucleotide mapped、3,033
  translated mapped 和 132,981 unmapped；得到 11,776 个 ordinary gene families
  与 147 个 pathways。native gene-family RPK 另按内置 `uniref90_rxn` mapping
  regroup 为 1,402 个 ordinary MetaCyc reactions；reaction strata 0 个违反加和
  约束，`UNGROUPED` 保留 113,812.639 RPK，one-to-many mapped expansion 为
  1.935×。147/147 个 pathway abundance 体现独立重建的非加和语义。主流程实测
  1,861.48 秒、6.62 GiB peak RSS；一次性初始化冻结 34 个 payload，离线验证
  202/202 全部通过。
- 第 20 篇把第 19 篇冻结的 gene-family、reaction、pathway abundance 与 coverage
  表拆成 3 个 feature spaces、2 种单位、2 种归一化模式和 2 种 special-row
  policy，共审计 24 个分支，closure failure 为 0；其中 4 个 pathway 分支由
  HUMAnN 3.9 `humann_renorm_table` 实际生成并与独立重算一致。真实队列部分固定
  curatedMetagenomicData 3.12.0 的 AsnicarF_2017 pathway abundance/coverage：
  11,173 行、24 个 profiles、15 名受试者和 445 个 ordinary unstratified pathways，
  仅用于阈值、重复测量和零值敏感性示范，不做疾病组推断；离线验证 110/110 通过。
- 第 21 篇联合第 13、15、16、19、20 篇的真实冻结证据与同一 MOCK1 的
  MicrobeCensus 校准，把 observed fragments、model estimates、RPK、CoPM、
  relative abundance、coverage、average genome size、genome equivalents 与 RPKG
  写成 16 行逐列语义合同，并预先裁决 18 种下游变换。150 bp 主分支在
  29,809,773 bp 中估计 AGS 2,576,205.26 bp 和 11.5712 genome equivalents；
  100 bp 敏感性分支为 2,469,448.80 bp 和 12.0714 genome equivalents。
  真实 AsnicarF_2017 的 10,680 个 pathway-profile cells 被拆成四种 abundance ×
  coverage zero 状态；不做生物分组检验。离线验证 252/252 通过，四张英文图均已检查。
- 第 22 篇用同批 24 个 AsnicarF_2017 profiles / 15 名受试者的 298-species
  MetaPhlAn 表与 HUMAnN gene-family 表，外加独立的 500-sample / 56-spring /
  780-MAG 温泉 catalog，统一审计 Hill q=0/1/2、Bray–Curtis、binary Jaccard、
  Aitchison 与 Cailliez-corrected PCoA。gene-family 主分析预先固定 prevalence
  ≥5/24（178,928 features），MAG 组成只在 recovered catalog 内闭合，并把平均
  12.9828% reads recruitment 单列为覆盖边界；不执行生物分组检验或置换。
  离线验证 79/79 通过，论文 CC BY 原图和四张英文重绘图均已检查。
- 第 23 篇把同一温泉数据的 500 个局部样本等权聚合为 56 个 hot-spring
  推断单位，在 780-MAG recovered catalog 上固定 Bray–Curtis、marginal
  PERMANOVA、BroadRegion 内 9,999 个唯一非原位置换和 spatial-median PERMDISP。
  Temperature regime 的 total/partial R² 为 0.1434/0.2114，restricted
  *p* = 0.0001；PERMDISP *F*(4, 51) = 4.296、*p* = 0.0043，因此不写成纯
  centroid shift。十个 pairwise 对比中 6 个可估、4 个 not estimable、2 个
  Holm rejection；九个预声明敏感性分支全部保留。离线验证 123/123 通过，
  论文 CC BY 原图和四张英文重绘图均已检查。
- 第 24 篇在 110 名独立 ZellerG_2014 CRC/Control 受试者上分开建模 abundance
  与 prevalence，并平行审计 MaAsLin3、ANCOM-BC2 和 ALDEx2。物种与通路主
  feature spaces 分别为 212 和 394 项；离线验证 100/100 通过。
- 第 25 篇从 HMP_2012 的 748 个 profiles 固定 490 个 subject-habitat
  representatives，按 rank 分别闭合，并以 Wilson 区间、1,000 次受试者 bootstrap
  和 112 个阈值分支区分 point core 与 stable core；离线验证 97/97 通过。
- 第 26 篇审计 curatedMetagenomicData 3.12.0 的 package snapshot、论文集合、
  732 个 date-stamped resources、177 个跨研究重名 sample IDs，并真实合并三队列
  923 species × 261 biological units；离线验证 157/157 通过。
- 第 27 篇在 114 名 ZellerG_2014 独立受试者上执行 5×5 外层、4 折内层的
  Random forest/XGBoost nested CV；安全置换 median AUROC 为 0.4963，而故意全局
  选特征的泄漏分支为 0.7083；离线验证 122/122 通过。
- 第 28 篇使用同一 cMD 2021-03-31 lineage 的 8 个 CRC 队列、771 名独立受试者
  和 897 个联合物种。五队列随机效应 meta-analysis 与三队列方向验证独立于模型；
  8 轮真正 LODO 的 macro AUROC 为 0.7786（层级 bootstrap 95% CI 0.7107–0.8324），
  random-effects pooled AUROC 为 0.7859、I² 为 52.09%；离线验证 112/112 通过。
- 第 29 篇复用第 23 篇 56 个 spring-level 单位与 780-MAG catalog，主筛选保留
  93 个 MAG；CLR 后校正 BroadRegion、temperature 与 pH，再以 huge graphical
  lasso + StARS 得到 70 条 conditional associations。1,000 次分层 bootstrap
  支持其中 30 条边达到 0.70；7 个 Louvain modules 中没有达到预声明 Zi hub
  阈值的节点，因此仅报告 5 个 topology-priority follow-up candidates，不称为
  ecological keystones。度保持零模型与 1,000 次随机删除对照均已完成；离线验证
  82/82 通过，论文 CC BY 原图和五张英文重绘图均已检查。
- 第 30 篇从 PRJEB52977 的 MOCK1/2 官方 Illumina runs 各精确抽取 2,000,000
  个同步 pairs，fastp 后保留 1,999,853/1,999,888 pairs；实际运行 MEGAHIT
  1.2.9 与 metaSPAdes 4.3.0 的四个 single branches 和两个 co-assembly branches，
  再用 Bowtie2 2.5.5 完成八次 read-back mapping。本次 6 个组装的 ≥1 kb
  assembled sequence 为 54.48–88.17 Mb、N50 为 5,906–17,124 bp，回帖率为
  83.79%–96.58%；这些指标只描述本次工作流的表示度与资源权衡，不替代
  reference-aware correctness。22 条任务全部正常退出，66 个冻结 payload
  checksum 与 194/194 项离线检查全部通过，四张英文图均已检查。
- 第 31 篇使用同一 MOCK1 DNA 的完整 `ERR9765780`（696,944 条历史 ONT R9
  reads；3.126 Gbp）与 `ERR9765783`（524,805 条 PacBio HiFi reads；5.400 Gbp），
  实跑 Flye 两分支、hifiasm-meta `--force-rs` 参数敏感性分支与 metaMDBG。
  四套 `≥1 kb` assemblies 为 143.77–179.43 Mb，N50 为 0.441–2.155 Mb；
  minimap2 2.31 `-c` read-back 的 aligned input bases 为 71.82%–99.04%，
  base-level weighted identity 为 89.57%–99.799%。152 条 `≥10 kb` circular
  candidates 中 118 条得到至少 3 条跨首尾 reads 支持；hifiasm-meta 的峰值
  RSS 为 135.72 GiB。68 个冻结 payload checksum 与 333/333 项离线检查全部
  通过，四张英文图均已人工检查；这些结果不构成平台单因素或普适工具排名。
- 第 32 篇从 `ERR9765746` 无放回同步抽取恰好 10,000,000 对 Illumina reads，
  与完整 `ERR9765780` ONT R9、`ERR9765783` HiFi reads 建立七分支比较。
  MetaQUAST 5.3.0 对 71 个 truth genomes 产生 497 行标准化结果；HiFi-only
  的 N50 为 2.014 Mb、genome fraction 为 70.125%、36 个 genomes 达到
  `>=99%`，mismatch/indel 为 8.11/10.42 per 100 kbp。ONT draft 经 default
  Polypolish 后 mismatch 与 indel 分别下降 20.1% 与 42.7%，但 `>=99%`
  回收数仍为 24。209 个固化 payload checksum 与 447/447 项离线检查全部通过，
  四张英文图已人工检查；结论限定为本 mock 与锁定输入下的 strategy comparison。
- 第 33 篇复用第 30–32 篇 15 套 checksum-identified biological assemblies，
  对 MOCK1 的 11 个工作流和 2 个诊断正控使用 71-genome truth，对 MOCK2 single
  与 MOCK1+MOCK2 co branches 使用同一 87-genome denominator。QUAST/MetaQUAST
  5.3.0 共生成 17 条统一指标和 1,271 行逐 reference 记录。Fragmentation 正控
  在保留全部 bases 与 43/36 个 `>=90%`/`>=99%` genomes 时把 N50 从 2.014 Mb
  降至 50 kb；block-rotation 正控保持 N50/L50 不变，却把 misassemblies 从 284
  提高到 323、NA50 降低 13.3%。45 个冻结 payload checksum 与 224/224 项离线
  检查全部通过，四张英文图已人工检查；正文不声明普适阈值或 N50 correctness。
- 第 34 篇复用第 30 篇 PRJEB52977 MOCK1/MOCK2 的六套 `>=1 kb` assemblies，
  用 Prodigal 2.6.3 预测 441,407 个 assembly ORFs，并把 87 个 exact mock genomes
  的 270,679 个 callable ORFs 压缩为 260,868 个 95%/95% truth clusters。主分析
  采用 MMseqs2 9.d36de 的两阶段 individual-plus-co membership 展开，MEGAHIT mix
  目录用 93,782 个 representatives 表示 216,191 个 raw ORFs；truth recovery 为
  24.741%，catalog support 为 92.633%。CD-HIT 4.8.1 与三组阈值分支作为方法/
  参数敏感性，不把目录大小写成基因丰富度。44 个冻结 payload checksum 与
  132/132 项离线检查全部通过，四张英文图已人工检查。
- 全书已在固定 Quarto 环境中完成 77/77 页面本地渲染。
- 正式全流程 QA 已通过：16 个来源均完成下载与 checksum 核对，51 个执行步骤和
  756 条断言全部通过；隔离 staging 中生成并核对了 77/77 个章节页面。
- `tutorial.yaml` 当前列出 77/77 篇已验证文章；第 41–44 篇也已纳入正式执行步骤，
  不再依赖发布目录中的历史产物。
- 77 篇本地公众号审阅包均已生成并通过严格公开内容审计；生成器不调用公众号
  发布或群发接口，账号侧草稿操作仍是单独、需明确授权的发布阶段。
- bioBakery、assembly/binning、CheckM2、GUNC、CheckM1、GTDB-Tk、dRep、
  病毒组、DRAM/METABOLIC、gapseq/CarveMe，以及第 61–77 篇所需的统计与整合
  环境均保存锁版证据、固化产物或章节级验证记录。

## 第 01–44 篇真实数据和执行证据

- 公开记录：[Zenodo 3517209](https://doi.org/10.5281/zenodo.3517209)
- 文件：`meta_all.tsv`
- 版本：1.0
- MD5：`da18b10fdabae6308329e80b73991f84`
- 许可：CC BY 4.0
- 来源论文：[Wirbel et al., *Nature Medicine*, 2019](https://doi.org/10.1038/s41591-019-0406-6)
- 文章使用范围：5 个发现/meta-analysis 队列、3 个独立验证队列，
  共 768 个 CRC/control 样本

`data/small/01-crc-cohort-summary.tsv` 是可直接渲染的小表；
`scripts/prepare_intro_data.R` 可从固定源重新生成并核对全部计数。
`scripts/prepare_article03_data.R` 从同一源生成
`data/small/03-crc-design-audit.tsv`，QA 会重新生成并逐字节核对。
`scripts/prepare_article04_data.R` 从同一源生成 768 行
`data/small/04-crc-library-size.tsv`，并从 Nonpareil `v3.5.5`
包内未修改的 `LakeLanier.npo` 重建
`data/small/04-lake-lanier-coverage.tsv`。原始 `.npo` SHA-256 为
`8db00c35eadb64c36eaccad35007c9e7a5d553efb0d5c12b6f859ab1097f52f4`，
第三方数据许可与 immutable source 见
`data/small/04-nonpareil-data-NOTICE.txt`。

第 05 篇从 Costea 等人的公开仓库固定提交
`24c0edd557f18f0b62b8fd46a60fe3819515227a` 读取 Phase III
样本信息、MetaPhlAn2 谱和理论 mock composition，三个源文件均由
`tutorial.yaml` 锁定 SHA-256。`scripts/prepare_article05_data.R`
生成 27 个文库、10 个物种的偏差矩阵，并把 Zaramela 等人的 synDNA
Table 1 数值固化为定量基准；来源、许可与再生边界见
`data/small/05-data-NOTICE.txt`。

第 06 篇把 Marotz 等 Results/Figure 2 的 6 种唾液处理，以及 Longhi 等
main Table 2/Supplementary Table 6 的 7 个 saponin 条件固化为 CC BY
来源表。`scripts/prepare_article06_data.R` 重算 read budget 与组成权衡，
`scripts/render_article06_figures.R` 一次生成三张 PDF/PNG/TIFF 图。Longhi
剂量系列来自 1 份痰液，只作机制性风险提示；字段映射、补充文件 checksum、
accession 与许可见 `data/small/06-data-NOTICE.txt`。

第 07 篇把 `decontam::MUClite` 真实 16S V4 对照对象导出为 byte-locked
三件套，仅用于执行 frequency/prevalence classifier；shotgun 机制结论另由
Salter 等人的真实 *S. bongori* dilution experiment 锚定，不把 ASV 表冒充
shotgun 数据。`scripts/prepare_article07_data.R` 固定 `decontam 1.24.0`
并导出 feature 判定、plate/threshold sensitivity、逐样本 burden 和过滤后
三件套；Costello 与 Illumina 的 index-hopping 数值保持情境边界，不作为统一
合格阈值。hash、许可与转录范围见 `data/small/07-data-NOTICE.txt`。

第 08 篇固定 Meslier 等 CC BY 4.0 benchmark 的 project `PRJEB52977`、
sample `SAMEA14435832`，以及 Illumina `ERR9765746`、historical ONT R9
`ERR9765780`、PacBio CCS `ERR9765783`。一次性 builder 从四个
checksum-identified ENA FASTQ 流式读取各平台前 5,000 条完整 records，只提交
15,000 行 per-read metrics、四条 60-base display prefixes 与 lineage hash，
不保存约 10.8 GB 原始 FASTQ。论文全量指标与非随机 prefix 分表保存；
paper Table 4 的 `ERR9765446` 与 official ENA `ERR9765746` discrepancy、
historical ONT 边界和原生 POD5/BAM 保留要求见
`data/small/08-data-NOTICE.txt`。

第 09 篇在 Ubuntu 22.04.5 原生 Linux 主机上运行
`scripts/validate_article09_environment.py`，使用
`env/platform-smoke.yml` 创建
`metagenome-platform-smoke-2026.07` 专用环境，并逐项核对 Python 3.12.13、
SeqKit 2.10.0、pigz 2.8、matplotlib-base 3.10.5 的来源与 build。真实数据
smoke test 读取第 08 篇已提交的 15,000 行 metrics 和 lineage summary；
Windows/WSL2 控制面保留为 `NOT_RUN`，需要读者在 Windows PowerShell 中另行留证。

第 10 篇复用同一 MOCK1 的四条 ENA source records 与 15,000 行真实 read-prefix
metrics。`scripts/validate_article10_compute.py` 在本地首轮完成 Illumina、ONT、
PacBio 三个 5,000-row task，第二轮核对 JSON SHA-256 后跳过 3/3；这只验证
scratch、原子发布和恢复语义，不用于外推完整 FASTQ 或 assembler 资源。
容器合同固定 Apptainer 1.5.2 与 linux/amd64 OCI manifest digest
`sha256:eafc1edb577d2e9b458664a15f23ea1c370214193226069eb22921169fc7e43f`。
本机 `apptainer`、`sbatch`、`sacct` 均不存在，legacy Singularity 3.7.2 的
SIF exec 失败记为 `EXPECTED_FAIL`，没有冒充目标集群实测。

第 11 篇在 native Linux 上创建 `metagenome-biobakery-2026.07` 与
`metagenome-assembly-2026.07`。完整 explicit locks 分别包含 554 与 239 个
SHA-256 标识的归档；`data/small/11-environment-evidence.tsv` 固化 15 个真实
入口和 3 个安装测试，`data/small/11-install-self-tests.log` 保留两类入口回归。
MetaPhlAn vJan26、CheckM2 v3 与 GTDB R232 具有发布方 MD5，HUMAnN
ChocoPhlAn v201901_v31 与 UniRef90 v201901b 因发布方未给 checksum 而默认阻断。
所有数据库状态仍为 `NOT_DOWNLOADED`，因此没有声称 profile、`checkm2 testrun`
或 `gtdbtk check_install` 已通过。

第 12 篇使用 Bioconductor 3.19 的 curatedMetagenomicData 3.12.0。一次性
retrieval 以 pattern `AsnicarF_2017.relative_abundance` 找到两个日期资源，
固定最新的 `2021-10-14.AsnicarF_2017.relative_abundance`（ExperimentHub
`EH7091`）。原始 Hub cache 为 18,187 bytes，SHA-256 为
`ad631532fbbab39dfc3550a676a58310887d317e947392dcae9e0d4e4d69db27`；
完整 RDS 快照为 199,904 bytes，SHA-256 为
`2952774730bff2af9e13c9c40058320aed524dc2dc7408b44dc6e697c06564b2`。
对象包含 298 个分类 features、24 个样本、15 位受试者、一个
`relative_abundance` assay、7 个 taxonomy 字段和 22 个样本 metadata 字段。
在线对象与 `LOCAL = TRUE` cache replay 完全一致；日常 QA 只读 SHA-256 锁定
快照，不访问 ExperimentHub。该百分比谱来自 cMD3 的 MetaPhlAn 3 lineage，
不冒充本系列 MetaPhlAn 4 上游输出。

第 13 篇复用同一 MOCK1 的 Illumina paired-end run `ERR9765746`，streaming
builder 从 ENA checksum-identified R1/R2 各取前 100,000 条完整 records，
逐条核对 normalized ID 与顺序。原始和 clean FASTQ 只留在 ignored scratch；
`data/small/13-qc-frozen/` 固化 raw/clean FastQC zip/HTML、fastp JSON/HTML、
MultiQC report/data、完整命令、资源日志和 57-file SHA-256 manifest。实测
combined Q30 从 89.9032% 变为 89.9034%，四份 FastQC 合计 36 PASS、8 WARN、
0 FAIL；GC 与 length-distribution WARN 保留为复杂 mock 与真实读长结构的诊断
信号，不自动解释为污染。数据、许可、hash、完整 ENA MD5 的证据边界与资源使用见
`data/small/13-data-NOTICE.txt`。

第 14 篇把公开 NA12878 WGS run `ERR194147` 仅作为去宿主阳性方法控制，
并把 `ERR9765746` MOCK1 作为独立的微生物保留控制；两者各使用前 20,000 个
同步 pairs，不把非随机前缀冒充全量估计。一次性运行使用 Hostile 2.0.2、
Bowtie2 2.5.4、Samtools 1.21 和官方 `human-t2t-hla` 2023-07 索引
（归档 3,934,284,979 bytes，SHA-256
`5b584f5c28abeec5dba78bd37b53fa476dd42af57051d2fb7d2f2098e3a2df13`）。
`data/small/14-host-removal-frozen/` 保存路径归一化的命令、资源、双控制守恒账本、
复杂度阈值敏感性和序列精确重复审计，checksum manifest 覆盖 47 个 payload；
不包含 FASTQ 或原始人源 read ID。主流程只保留双端均未比对到人参考的 pair，
低复杂度过滤和去重默认关闭。完整来源、许可、索引文件 hash 与隐私边界见
`data/small/14-data-NOTICE.txt`。

第 15 篇复用第 13 篇 checksum-locked clean FASTQ，固定完整 MetaPhlAn vJan26
metadata 与 Bowtie2 archives 的字节数、MD5、SHA-256 和提取文件清单。数据库本体、
FASTQ 和 mapout 留在 ignored storage；`data/small/15-metaphlan-frozen/` 只固化
路径归一化命令/日志、profile、species/SGB/marker/depth/threshold 表与 23-file
checksum manifest。20k、50k 和 99,991-pair 分支分别检出 13、25 和 31 个
species/SGB；0.1% 与 1% post-profile 阈值分别保留 29/99.832% 和
20/94.375% classified clades/abundance。VDB targets 与 PKL 排除的 SGB6011
index targets 独立记账，不进入 SGB marker breadth。完整来源、数据库 hash、
许可和 0% unclassified 解释边界见 `data/small/15-data-NOTICE.txt`。

第 16 篇继续使用同一 clean paired-FASTQ lineage，固定 Kraken2 Standard-8
2026-06-26 archive 的 release-specific URL、5,946,578,575-byte 大小、publisher
MD5 `7685f43cce057c2ca18511c925399b72`、本地 SHA-256 和 17 个解包文件 MD5。
5.54-GiB archive、8.64-GB index 与 per-fragment Kraken output 均留在 Git 外；
`data/small/16-kraken-bracken-frozen/` 只保存 standard report、四组 Bracken
outputs、分类/重分配/参数敏感性账本、路径归一化日志与 33-entry checksum
manifest。主分析对 99,991 个 pairs 报告 86.3808% classified；Bracken
species/150 bp/threshold 10 从 214 个非零 Kraken species rows 中保留 65 个，
而 threshold 0 保留 214 个。100 bp 分支相对主模型 total variation 为 0.00134，
genus 分支报告 54 个 taxa。数据库 cap、paired-fragment 单位、within-rank
fraction、未分类与无法分配质量边界见 `data/small/16-data-NOTICE.txt`。

第 17 篇固定三份 2026-06-26 release-specific archives 及 publisher MD5；
压缩包合计 22.24 GiB，解包索引合计 31.59 GiB，全部 17-file internal MD5
均通过。71 个 Meslier Supplementary Table S3 assemblies 经 2026-07-21
NCBI Datasets 18.33.1 快照压缩为 69 个 current species；每个库的 crosswalk
均为 62 个 exact assemblies、3 个 alternate same-species references 和 6 个
no same-species references。`data/small/17-kraken-database-confidence-frozen/`
只保存 standard reports、Bracken tables、truth/reference/control/resource audits
与 111-entry checksum manifest；FASTQ、archives、indexes 和 per-fragment outputs
留在 Git 外。human-control ledger 明确排除 *Homo sapiens* clade 与 human-lineage
ancestor direct calls 后才计算 unsupported burden；空 Bracken 分支保留真实 exit 1
和 header-only table，不插补 composition。完整边界见
`data/small/17-data-NOTICE.txt`。

第 18 篇复用第 13 篇 checksum-locked clean FASTQ、第 15 篇 MetaPhlAn
profile 与第 16–17 篇 Kraken/Bracken 证据，并实际运行 mOTUs 4.1.0。
mOTUs database 4.1 固定为 Zenodo record `20322482`：archive 为
5,552,983,255 bytes，publisher MD5 为
`471ea128f0c0839f5c4629b949ea5f8a`，本地 SHA-256 为
`7d2d6382ecf766b23ef362311715cd612243af82c00c936da4597afb4e4df375`；
解包后占 9,699,381,148 bytes，taxonomy 文件实际声明 GTDB R226。
一次映射后从同一 MGC table 计算 `g=1/3/6`，默认 `g=3`；BAM、MGC、FASTQ
和大型数据库留在 Git 外。`data/small/18-profiler-benchmark-frozen/` 保存
三套原生 profile、reference-aware crosswalk、共同域组成、资源账本、路径归一化
日志和 28-entry checksum manifest。严格共同一对一域为 52 species / 70.737791%
expected mass；默认分支 common-domain TV 分别为 0.405、0.363、0.323。
来源、许可和数据库计数差异见 `data/small/18-data-NOTICE.txt`。

第 19 篇继续复用第 13 篇的 99,991 个同步 clean pairs，但为 HUMAnN 3.9
重新生成兼容的 `mpa_vJun23_CHOCOPhlAnSGB_202403` prescreen；第 15 篇 vJan26
profile 的预期非兼容退出被作为 release sensitivity evidence 保留。两份
MetaPhlAn archives 使用 publisher MD5 验证；ChocoPhlAn full v201901_v31 与
UniRef90 annotated v201901b 的上游未发布 checksum，因此只声明 final HTTPS URL、
Content-Length、tar traversal 和本次 retrieval SHA-256，不冒充 publisher checksum。
`data/small/19-humann3-frozen/` 保存 native RPK/RPK-derived tables、community CPM/
relative-abundance branches、1,402-reaction `uniref90_rxn` regroup 表及 mapping/
expansion audit、未经归一化的 0–1 pathway coverage、read-flow、prescreen、database
inventory、资源账本、标准化日志、version-output compatibility wrapper 与 34-entry
checksum manifest；FASTQ、archives、database indexes、sample-specific index 和
per-read intermediates 留在 Git 外。reaction mapping 固定为 57,511,575 bytes，
SHA-256 `8419ce78a62ca9130914f2c347a9708111cedc7de52ba274659ce51ec7de7752`。
完整边界见 `data/small/19-data-NOTICE.txt`。

第 20 篇同时使用两条明确分开的证据线：方法审计复用第 19 篇 checksum-locked
MOCK1 HUMAnN 表；真实队列敏感性审计使用 curatedMetagenomicData 3.12.0 的
AsnicarF_2017 `pathway_abundance`（EH7089）和 `pathway_coverage`（EH7090）。两份
ExperimentHub 资源各为 11,173 × 24，映射到 15 名受试者；导出前按 pathway ID
显式对齐 abundance 与 coverage 的行顺序。`data/small/20-cmd-pathway/` 保存原始
`.rda`、派生表、metadata、资源 manifest、许可说明和 checksum；
`data/small/20-functional-normalization-frozen/` 保存 4 个真实 HUMAnN 3.9
`humann_renorm_table` pathway 输出、命令、日志、版本和 checksum。验证器独立重算
24 个归一化分支，审计 special-feature denominator、pathway contribution、profile/
subject prevalence、coverage threshold 和三个 pseudocount；所有图内文字均为英文。

第 21 篇复用第 13、15、16、19 篇与 AsnicarF_2017 队列的五组既有冻结证据，
再加入 `data/small/21-table-semantics-frozen/` 的九项 MicrobeCensus 运行证据；
验证时对六份 SHA-256 manifest 共 163 个 payload 逐项核对。一次性运行固定
MicrobeCensus tag v1.1.1 / commit
`dfc42d356bfd7943633cde6c0fbfc0b116f29ae2`（源码内部仍报告 1.1.0），只对
RAPsearch2 2.15 的 Python 3 bytes preflight 做兼容解码，不改 marker、阈值、权重
或 AGS 算法。99,991 pairs / 199,982 reads / 29,809,773 bp 的同一 clean sequence
universe 用于 150 bp 主分支和 100 bp 敏感性分支；FASTQ、数据库与逐 read 结果不进
常规 QA。`results/21-table-semantics/` 保存 16 行 table-semantics contract、18 项
transformation legality、closure/zero/strata/RPKG 审计和 252/252 validation summary。

第 22 篇的人体分支固定 curatedMetagenomicData 3.12.0 的 EH7091 species 表与
EH7086 gene-family 表；二者严格样本匹配为 24 profiles / 15 subjects。MAG 分支是
独立的 Korchagina et al. 2026 美国西部温泉数据，固定 Figshare record
`30284068` v2 的 metadata、BIOM 与 recruitment 三个文件及 publisher MD5。
`data/small/22-diversity-inputs/` 保存 9 项 checksum-locked 冻结输入；gene-family
原对象的 2,704,846 行先去除 taxon-stratified / special rows，再冻结 prevalence
≥3/24 的 415,581 行，主分析使用 ≥5/24 的 178,928 行，≥12/24 的 1,896 行仅作
敏感性分支。温泉表为 780 MAGs × 500 samples，独立 recruitment 平均
0.12982826194。`results/22-alpha-beta-diversity/` 保存数据谱系、alpha、beta、PCoA、
重复测量、filter/pseudocount/recruitment 敏感性和 79/79 validation summary；
所有 QA 复算均从冻结表开始且不访问网络。

第 23 篇复用上述 Figshare v2 数据谱系，但不把 500 个局部样本当成独立生态
重复。`data/small/23-ordination-permanova/` 保存 56 × 780 等权 spring-level
主矩阵、catalog-read-weighted 敏感性矩阵、spring metadata、sample-to-spring
ledger、温度分组和分析合同以及 7 项 SHA-256；`results/23-ordination-permanova/`
保存 aggregation、design balance、permutation space、PCoA/CAP、marginal/sequential
PERMANOVA、PERMDISP、pairwise estimability、九个 sensitivity branches 与 123/123
validation summary。主置换矩阵固定 seed 20260723，SHA-256 为
`94119a47186a7a786bc389917c0c57661a2a7b3488a6d2386b918bc71a4951b6`。

第 24–27 篇分别在 `data/small/24-differential-abundance/`、
`25-composition-core/`、`26-cmd-lineage/` 和 `27-machine-learning/` 保存真实
cMD 输入、分析合同与 SHA-256；对应 `results/` 目录保存逐方法、逐阈值、逐折
prediction、校准、重要性和完整验证记录。第 24 篇的方法环境另由
`env/differential-abundance.yml` 与 Linux-64 lock 固定，第 27 篇由
`env/machine-learning-renv.lock` 固定。

第 28 篇的八份官方 ExperimentHub RDA 均锁定 URL、release、bytes 与 SHA-256，
原始文件留在 `downloads/article28/`、不进 Git；
`scripts/prepare_article28_cross_cohort.R` 重建 897 × 771 fraction table、metadata、
cohort summary、feature audit 与 `analysis-contract.tsv`。完整 LODO、Hedges g/REML、
层级 bootstrap、8×8 transfer、coefficient stability 和 112 项检查由
`scripts/validate_article28_cross_cohort.R` 生成，R 依赖由
`env/cross-cohort-renv.lock` 固定。

第 29 篇从第 23 篇已校验的 spring-level 表重建独立输入；
`scripts/prepare_article29_network_data.R` 固化 56 × 780 组成表、环境元数据、
feature-filter audit、27 行分析合同和 5 项 SHA-256。完整 graphical-lasso/StARS、
1,000 次 BroadRegion 分层 bootstrap、Zi-Pi、1,000 个 degree-preserving nulls、
targeted-vs-random deletion 与敏感性分支由
`scripts/validate_article29_network.R` 执行；R 依赖由 `env/network-renv.lock` 固定。

第 30 篇固定 Meslier 等 PRJEB52977 的两个真实不均一 DNA mocks：
`ERR9765746`（MOCK1）与 `ERR9765747`（MOCK2）。
`scripts/download_article30_assembly_reads.sh` 按 ENA byte count 与 MD5 校验四个
FASTQ archives；`scripts/select_article30_read_pairs.py` 用 seed 20260730 做一次遍历、
无放回的精确同步抽样；原始和清洁 FASTQ 均不进 Git。
`scripts/run_article30_short_read_assembly.sh` 实跑 QC、六个组装、六个索引和八次
mapping，`scripts/freeze_article30_short_read_assembly.py` 固化版本、命令、资源、
日志、统计表与六套 ≥1 kb contigs，`scripts/validate_article30_short_read_assembly.py`
离线重算指标、核对 66 项 checksum 并生成四张三格式图。

第 31 篇固定同一项目、同一 `SAMEA14435832` MOCK1 DNA 的完整历史 ONT R9
`ERR9765780` 与 PacBio HiFi `ERR9765783` archives；bytes、MD5、SHA-256、FASTQ
grammar、read/base counts 均独立核对。`scripts/run_article31_long_read_assembly.sh`
运行 Flye、hifiasm-meta 和 metaMDBG 四分支，并以 minimap2 `-c` 完成四次
base-level read-back 与四次 junction mapping；首次 ONT/Flye 的本地 IPC 失败及
同参数 resume 均保留在 attempts 账本。`scripts/freeze_article31_long_read_assembly.py`
固化 68 个 checksum-covered payloads，其中四套 `≥1 kb` assemblies 共 193 MB；
原始 reads、PAF 与 16 GB 工作目录留在 Git 外。专用验证器离线复算组装指标、
identity、junction、资源和解释边界，并生成四张 PDF/PNG/350-dpi LZW TIFF 图。

第 32 篇在同一 MOCK1 DNA 上固定 10,000,000 对 Illumina reads 与完整
ONT R9/HiFi reads；SPAdes 4.3.0 建立 short-only 和两条 short-read-first
hybrid，Flye 2.9.6 基线建立两条 long-only，ONT draft 再进入
BWA-MEM 0.7.19 `-a` 与 Polypolish 0.6.1 default/careful。七套 `>=1 kb`
assemblies 用 MetaQUAST 5.3.0 在 71 references、500 bp minimum alignment 与
97% identity 下统一评价。`data/small/32-hybrid-assembly-polishing-frozen/`
保留 209 个 checksum-covered payloads，其中包含七套 assemblies、71 个独立
references、63 份实际生成的 per-reference reports 和 497 行标准化结果；
原始 FASTQ、SAM 与 31 GB work 目录留在 Git 外。专用验证器通过
447/447 项检查并生成四张 PDF/PNG/350-dpi LZW TIFF 图。

第 33 篇不再读取 FASTQ，而是逐个核对第 30–32 篇 15 套非重复压缩 FASTA 的
上游 SHA-256，再在统一 `>=1 kb` 空间运行一次 reference-free QUAST 和三次
MetaQUAST。MOCK1、MOCK2、MOCK1+MOCK2 分别使用 71、87、87 个 exact
references；两个从真实 Flye HiFi assembly 派生的确定性正控分别保持 base
multiset 或完整 contig-length multiset，用于验证 N50 与 structure 指标语义。
`data/small/33-assembly-qc-frozen/` 保留 45 个 checksum-covered 小表、报告、日志、
资源、环境和脚本，不包含 FASTQ 或 assembly FASTA。专用验证器通过 224/224 项
检查并生成四张 PDF/PNG/350-dpi LZW TIFF 图。

第 34 篇逐个核对第 30 篇六套 assembly FASTA 的上游 SHA-256，用 Prodigal
2.6.3 分别执行 metagenome-mode ORF prediction，并以 MMseqs2 9.d36de 建立
individual、co-assembly 与两阶段 mix catalogs。真实 truth 由同一 benchmark
repository 的 87 个 exact MOCK2 genomes 预测并按相同 95% identity / 95%
coverage 合并；CD-HIT 4.8.1 及 MMseqs2 90%/80%、99%/95% 分支只用于敏感性。
`data/small/34-nonredundant-gene-catalog-frozen/` 保留 44 个 checksum-covered
payloads，包括主目录 93,782 条 paired protein/nucleotide representatives、完整
expanded membership、小表、日志、资源、环境与脚本；11 GB 一次性工作目录留在
Git 外。专用验证器通过 132/132 项检查并生成四张 PDF/PNG/350-dpi LZW TIFF 图。

第 35 篇逐个核对第 30 篇 MOCK1/MOCK2 四份 clean FASTQ 与第 34 篇 paired
catalog 的上游 SHA-256。论文兼容分支固定抽取每样本 10,000 条 forward reads；
主分支用 Bowtie2 `--very-sensitive-local -k 2` 把 R1/R2 分别作为 reads，并以
MAPQ 10、95% identity 和 80% query coverage 建立 assigned ledger。Raw counts
另派生 CPM/RPKM/TPM；DIAMOND 2.2.4 以 iterative sensitive 模式搜索 36.3-GB
UniRef90 v201901b，再用 HUMAnN 3.9 reaction mapping 比较 equal-split 守恒与
copy-full 膨胀。89,339 个 catalog genes 获得 best hit，15,392 个 genes 连到
MetaCyc reaction；主 assigned reads 的 reaction-linked 比例为 21.4299%/21.4047%，
copy-full CPM 总量膨胀到 1.212374×/1.212982×。冻结目录包含 45 个
checksum-covered payloads 和 1 个 manifest，专项验证器通过 225/225 项检查并
生成四张 PDF/PNG/350-dpi LZW TIFF 图。FASTQ、索引、SAM/BAM 和数据库留在 Git 外。

第 36 篇逐个核对第 34 篇 93,782 条 primary protein representatives 与第 35 篇
两套 raw-count ledger 的 SHA-256 和 ID universe。主分支用 eggNOG-mapper 2.1.15、
eggNOG 5.0.2、DIAMOND 2.0.15、sensitive iterative search、automatic tax scope、
non-electronic GO evidence 与 ortholog-derived PFAM transfer；同一 84,511 条 seed
hits 再以 `all` GO evidence 重注释作敏感性分析。四个互斥状态分别包含 9,271、
6,243、24,965 和 53,303 条 genes；no-seed 的 gene fraction 为 9.8857%，而
MOCK1/MOCK2 assigned-read fractions 为 3.3867%/3.3049%。主分支耗时 11:58.62、
峰值 RAM 7.14 GiB，GO sensitivity 耗时 1:09.89。冻结目录含 45 个
checksum-covered payloads 加一个 manifest，专项验证器通过 136/136 项检查并生成
四张 PDF/PNG/350-dpi LZW TIFF 图；约 50.94 GB 的 installed database assets 与
12.06 GB archives 留在 Git 外。

第 37 篇复用同一 93,782-protein catalog 和两套 assigned-read ledger，固定
dbCAN 5.2.9、database `db_v5-2-9_5-5-2026`、DIAMOND/family HMM/dbCAN-sub
三条证据以及两工具主门槛。主集合含 2,050 条 genes，三工具共识 1,605 条，
730 条获得 substrate candidates；完整 *Bacteroides thetaiotaomicron* VPI-5482
正控检出 117 个 CGC。一次性 catalog/正控运行分别耗时 20:36/7:41、峰值 RAM
5.52/2.98 GiB；39 个 payload 由 checksum 覆盖，离线验证通过 84/84 项。

第 38 篇固定 RGI 6.0.8、CARD 4.0.1、DIAMOND 2.2.4、Perfect/Strict 主集合、
Loose sensitivity 与关闭 nudge 的合同。catalog 主集合含 36 条 genes，MOCK1/
MOCK2 分别承载 1,682/1,778 条 raw reads；co-assembly、*P. aeruginosa* 和
*S. aureus* 正控分别检出 34、32、21 条主 ORFs。40 个冻结 payload 由 checksum
覆盖，离线验证通过 85/85 项；CARD 数据库和中间 JSON 留在 Git 外。

第 39 篇固定 VFDB core/full 的 2026-07-24 snapshot、ABRicate 1.4.0、BLAST+
2.17.0 与 core 90% identity/80% reference coverage 主门槛。主 catalog 93 条
genes 在 MOCK1/MOCK2 承载 3,469/3,303 条 raw reads；core 80/80 和 full 90/80
分别回收 123/184 条 genes，阈值与数据库范围敏感性保持分离。38 个冻结 payload
由 checksum 覆盖，离线验证通过 83/83 项。

第 40 篇固定 antiSMASH 8.0.4、GECCO 0.10.3、Prodigal 2.6.3、PFAM 35.0、
MIBiG 4.0、MITE 1.3 与 25% reciprocal coordinate overlap。真实 co-assembly 的
antiSMASH/GECCO 分别检出 71/21 个 regions，各有 17 个得到另一工具坐标支持；
antiSMASH regions 链接 1,499 个 catalog genes，membership missing 为零，并在
MOCK1/MOCK2 承载 69,790/68,816 条 raw reads。八个工具分支均完成，最重分支
耗时 22:25、峰值 RAM 1.77 GiB；57 个冻结 payload 由 checksum 覆盖，离线验证
通过 115/115 项。

第 41 篇逐个核对第 30 篇 MEGAHIT co-assembly 和两套 clean FASTQ 的 SHA-256，
用 Bowtie2 2.5.5、SAMtools 1.23.1 与 MetaBAT2 2.18 的
`jgi_summarize_bam_contig_depths` 建立 full-coordinate depth ledger。MOCK1/
MOCK2 的 length-weighted JGI mean depth 为 5.812685/5.786746，18,353 条 contigs
在两个样本均检出；log-depth Pearson/Spearman 为 0.791678/0.483371。冻结目录
包含 69 个文件，专项验证器通过 215/215 项检查并重绘三套三格式发表图。

第 42 篇把同一 co-assembly 固定到 `>=1.5 kb` 的 10,203 条 contigs，并把两列
深度与 Kraken Standard-8 taxonomy 锁定后再启动五个 binner branches。MetaBAT2
single/multisample、SemiBin2、VAMB、TaxVAMB 分别产生 33/34/43/26/35 个 bins；
独立 CheckM2+GUNC gate 分别保留 22/22/24/15/18 个。两个 MetaBAT2 分支共享
8,052 条已分箱 contigs，ARI 为 0.981018。冻结目录包含 71 个文件，专项验证器
通过 361/361 项检查；mock truth 和 Kraken taxonomy 的使用边界分别记录。

第 43 篇从 171 个源 bins 重建 exact memberships，用 DAS Tool 1.1.7 与 Binette
1.2.1 搜索 non-overlapping candidate sets，并对全部 53 个 candidates 独立重跑
CheckM2/GUNC。DAS Tool 与 Binette 分别产生 26/27 个 candidates，其中 21/23 个
通过预注册 minimum gate；因此 truth-blind 主规则选择 Binette 的 23 个 MAG，
passing score sum 为 2,054.76。冻结目录包含 61 个文件，专项验证器通过 208/208
项检查并重绘四套发表图。

第 44 篇对这 23 个 truth-blind candidates 独立运行 CheckM2 1.1.0、GUNC 1.1.0、
CheckM 1.2.5、barrnap 1.10.5、tRNAscan-SE 2.0.13 与 Prodigal 2.6.3。全部通过
GUNC，17 个通过 CheckM2 `>90%/<5%` 前门，但只有 5 个含完整 rRNA set、21 个
达到至少 18 个 tRNA isotypes，最终按完整 MIMAG 分为 4 个 high-quality 与
19 个 medium-quality。CheckM marker lineage 决定 16 个 bacterial、7 个 archaeal
annotation modes；58,070 个 FASTG nodes 全部以 exact-sequence SHA-256 映回
k141 坐标，共解析 5,843 条 edges。冻结目录含 351 个文件、350 个 checksum
payload，专项验证器通过 581/581 项检查并重绘四套三格式发表图。

## 固定环境

- R 4.4.1
- Quarto 1.9.38
- `env/renv.lock`：R 4.4.1 / Bioconductor 3.19 的 185 条精确包记录，覆盖
  第 01–08 篇绘图与第 12 篇 curatedMetagenomicData 对象生态
- `env/platform-smoke.yml`：第 09–10 篇固定 Linux 环境审计与本地 smoke 的直接依赖
- `env/biobakery.yml` + `env/biobakery-linux-64.lock`：MetaPhlAn 4.2.5、
  HUMAnN 3.9、Bowtie2 2.5.5 与 DIAMOND 2.2.4
- `env/assembly.yml` + `env/assembly-linux-64.lock`：Python 3.10.20、MEGAHIT
  1.2.9、SPAdes/metaSPAdes 4.3.0、Bowtie2 2.5.5，以及后续章节使用的
  MetaBAT2、MaxBin2、CONCOCT、DAS Tool 与 CoverM
- `env/long-read-assembly.yml` + `env/long-read-assembly-linux-64.lock`：
  Python 3.10.20、Flye 2.9.6、metaMDBG 1.4、minimap2 2.31、Samtools 1.23.1、
  SeqKit 2.11.0 与 104-package Linux-64 transaction；hifiasm-meta 0.3.5-r81
  另由 tag/commit/source SHA-256 锁定后编译
- `env/hybrid-assembly.yml` + `env/hybrid-assembly-linux-64.lock`：fastp 1.3.6、
  SPAdes 4.3.0、BWA 0.7.19、Polypolish 0.6.1、MetaQUAST 5.3.0 与发表图绘制栈
- `env/gene-catalog.yml` + `env/gene-catalog-linux-64.lock`：Prodigal 2.6.3、
  MMseqs2 9.d36de、CD-HIT 4.8.1、SeqKit 2.11.0 与 Python 3.10.14
- `env/gene-abundance.yml` + `env/gene-abundance-linux-64.lock`：Bowtie2 2.5.5、
  SAMtools 1.23.1、HTSeq 2.1.2、seqtk 1.5、DIAMOND 2.2.4、Python 3.10.20
  及固定的 NumPy/Matplotlib/Pillow 作图栈
- `env/eggnog-annotation.yml` + `env/eggnog-annotation-linux-64.lock`：
  eggNOG-mapper 2.1.15、DIAMOND 2.0.15、Python 3.11.9，以及固定的
  NumPy 2.4.6、Matplotlib 3.10.1 与 Pillow 10.3.0 作图栈
- `env/cazyme.yml` + `env/cazyme-linux-64.lock`：dbCAN 5.2.9、DIAMOND
  2.2.4、pyHMMER、Pyrodigal 与固定 Python 作图栈
- `env/resistome.yml` + `env/resistome-linux-64.lock`：RGI 6.0.8、DIAMOND
  2.2.4、Prodigal 2.6.3 与 CARD 4.0.1 运行依赖
- `env/virulome.yml` + `env/virulome-linux-64.lock`：ABRicate 1.4.0、
  BLAST+ 2.17.0 与 VFDB core/full screening 依赖
- `env/bgc-gecco.yml` + `env/bgc-gecco-linux-64.lock`：GECCO 0.10.3；
  `env/antismash8.yml` + `env/antismash8-linux-64.lock`：antiSMASH 8.0.4、
  Prodigal 2.6.3 及其独立数据库合同
- `env/binning.yml` + `env/binning-linux-64.lock`：SemiBin2 2.3.0、VAMB/
  TaxVAMB 5.0.4、Kraken2 2.17.1 与 Python 3.12；MetaBAT2/DAS Tool 从锁定的
  assembly 环境调用
- `env/mag-qc.yml` + `env/mag-qc-linux-64.lock`：Binette 1.2.1、CheckM2
  1.1.0、barrnap executable 1.10.5、tRNAscan-SE 2.0.13 与 Prodigal 2.6.3
- `env/gunc.yml` + `env/gunc-linux-64.lock`：GUNC 1.1.0、DIAMOND 2.1.24
  与 ProGenomes 2.1 数据库运行依赖
- `env/checkm1.yml` + `env/checkm1-linux-64.lock`：CheckM 1.2.5 与
  2015-01-16 reference data；pplacer 固定单线程以控制内存
- `env/relink-biobakery-entrypoints.sh`：重建 bioBakery 后恢复并验收被 HUMAnN
  同路径覆盖的 Bowtie2/DIAMOND 精确入口
- `env/read-qc.yml` + `env/read-qc-linux-64.lock`：FastQC 0.12.1、fastp
  1.3.6、MultiQC 1.35、Python 3.14.6 与 171-package Linux-64 transaction
- `env/host-removal.yml` + `env/host-removal-linux-64.lock`：Hostile 2.0.2、
  Bowtie2 2.5.4、Samtools 1.21、fastp 1.3.6、SeqKit 2.10.0 与
  142-package Linux-64 transaction
- `env/kraken.yml` + `env/kraken-linux-64.lock`：Kraken2 2.17.1、
  Bracken package 3.1p1、Python 3.12.13 与 187-package Linux-64 transaction
- `env/motus.yml` + `env/motus-linux-64.lock`：mOTUs 4.1.0、BWA
  0.7.19-r1273、VSEARCH 2.31.0、Python 3.12.13 与 97-package Linux-64 transaction
- `env/differential-abundance.yml` + `env/differential-abundance-linux-64.lock`：
  第 24 篇 MaAsLin3、ANCOM-BC2 与 ALDEx2 独立方法环境
- `env/machine-learning-renv.lock`：第 27 篇 ranger、xgboost、pROC 与绘图依赖
- `env/cross-cohort-renv.lock`：第 28 篇 glmnet、metafor、pROC 与绘图依赖
- `env/network-renv.lock`：第 29 篇 huge、igraph、ggplot2 与网络审计依赖

恢复 R 环境：

```bash
Rscript -e 'install.packages("renv", repos = "https://cloud.r-project.org")'
Rscript -e 'renv::restore(lockfile = "env/renv.lock", prompt = FALSE)'
```

## 验证入口

检查 manifest：

```bash
python ../skills/best-practice-tutorial-style/scripts/render_tutorial.py \
  tutorial.yaml --variant github --check-only
```

检查 77 篇目录与第 01–44 篇结构：

```bash
python scripts/validate_series.py \
  --project-root . \
  --manifest tutorial.yaml \
  --output scaffold_validation.json
```

执行 Pilot QA：

```bash
python ../skills/tutorial-execution-qa/scripts/run_tutorial_qa.py \
  tutorial.yaml \
  --workspace .tutorial_runs/articles-01-44-local-v1 \
  --output qa_report.json
```

直接渲染网站：

```bash
quarto render
```

从已通过 QA 的网页生成与 16S 系列同款的本地公众号预览：

```bash
python scripts/build_wechat_review_bundle.py \
  --project-root . \
  --manifest tutorial.yaml \
  --qa-report qa_report.json \
  --article-number 1 \
  --article-number 23
```

不传 `--article-number` 时生成完整 77 篇审阅包。每篇目录包含
`article.html`、`draft.json`、`cover.jpg` 和优化后的正文图片；生成器只做本地派生，
不会上传图片、创建草稿、发布或群发。正文沿用 16S 系列的绿灰/米色排版、
figure-only 封面和代码保留规则：删除通用安装与主题函数块，保留会改变分析决策的
代码和全部有意义的结果图。

公众号草稿标题由生成器统一派生为
`宏基因组最佳实践｜N. 主题`。`N` 使用不补零的章节编号；QMD 与网页继续保留
聚焦科学问题的主题标题。标题上限固定为 64 个字符，只有超长主题才在
`tutorial.yaml` 中用 `wechat_title` 缩写，系列前缀和编号不缩写。

审计 77 篇完成度：

```bash
python scripts/audit_chapter_completion.py \
  --project-root . \
  --manifest tutorial.yaml \
  --qa-report qa_report.json \
  --output completion_report.json
```

## 主要入口

- `tutorial.yaml`：内容、执行与发布契约
- `_quarto.yml`：77 篇 Quarto Book 目录
- `index.qmd`：第 01 篇正式文章
- `chapters/02-three-analysis-layers.qmd`：第 02 篇正式文章
- `chapters/03-study-design-power.qmd`：第 03 篇正式文章
- `chapters/04-sequencing-depth.qmd`：第 04 篇正式文章
- `chapters/05-library-prep-absolute-quantification.qmd`：第 05 篇正式文章
- `chapters/06-host-depletion-low-biomass.qmd`：第 06 篇正式文章
- `chapters/07-contamination-controls.qmd`：第 07 篇正式文章
- `chapters/08-fastq-short-long-reads.qmd`：第 08 篇正式文章
- `chapters/09-wsl2-conda.qmd`：第 09 篇正式文章
- `chapters/10-computing-hpc-cloud.qmd`：第 10 篇正式文章
- `chapters/11-install-biobakery-assembly.qmd`：第 11 篇正式文章
- `chapters/12-install-r-cmd.qmd`：第 12 篇正式文章
- `chapters/13-read-qc-fastp.qmd`：第 13 篇正式文章
- `chapters/14-host-removal-complexity-duplicates.qmd`：第 14 篇正式文章
- `chapters/15-metaphlan4.qmd`：第 15 篇正式文章
- `chapters/16-kraken2-bracken.qmd`：第 16 篇正式文章
- `chapters/17-kraken2-database-confidence.qmd`：第 17 篇正式文章
- `chapters/18-profiler-benchmark.qmd`：第 18 篇正式文章
- `chapters/19-humann3.qmd`：第 19 篇正式文章
- `chapters/20-functional-profile-normalization.qmd`：第 20 篇正式文章
- `chapters/21-metagenomic-table-semantics.qmd`：第 21 篇正式文章
- `chapters/22-alpha-beta-diversity.qmd`：第 22 篇正式文章
- `chapters/23-pcoa-cap-permanova.qmd`：第 23 篇正式文章
- `chapters/24-differential-abundance.qmd`：第 24 篇正式文章
- `chapters/25-composition-core-microbiome.qmd`：第 25 篇正式文章
- `chapters/26-curated-metagenomic-data-lineage.qmd`：第 26 篇正式文章
- `chapters/27-machine-learning-roc.qmd`：第 27 篇正式文章
- `chapters/28-cross-cohort-validation.qmd`：第 28 篇真正 LODO 正式文章
- `chapters/29-cooccurrence-network.qmd`：第 29 篇条件关联网络与拓扑审计正式文章
- `chapters/30-short-read-assembly.qmd`：第 30 篇短读长 single/co 六分支组装正式文章
- `chapters/31-long-read-assembly.qmd`：第 31 篇 ONT/HiFi 四分支长读组装正式文章
- `chapters/32-hybrid-assembly-polishing.qmd`：第 32 篇七分支 hybrid assembly、Polypolish 与 truth-aware 审计正式文章
- `chapters/33-assembly-qc.qmd`：第 33 篇 QUAST/MetaQUAST、N50 正控与 task-specific usability 正式文章
- `chapters/34-nonredundant-gene-catalog.qmd`：第 34 篇 Prodigal、MMseqs2/CD-HIT、两阶段 membership 与 truth-aware 基因目录正式文章
- `chapters/35-gene-abundance.qmd`：第 35 篇 read-to-gene ledger、CPM/RPKM/TPM 与 UniRef90-to-MetaCyc 功能汇总正式文章
- `chapters/36-eggnog-functional-annotation.qmd`：第 36 篇 eggNOG 字段覆盖、四态功能暗物质、GO evidence 与多标签守恒正式文章
- `chapters/37-cazymes-dbcan.qmd`：第 37 篇 dbCAN 三工具、CAZy family/class、底物与 CGC 正式文章
- `chapters/38-resistome.qmd`：第 38 篇 CARD-RGI 证据分级、正控与 resistome 表型边界正式文章
- `chapters/39-virulome.qmd`：第 39 篇 VFDB core/full、阈值敏感性与 virulome 表型边界正式文章
- `chapters/40-bgc-natural-products.qmd`：第 40 篇 antiSMASH/GECCO、fragmentation、MIBiG 与 BGC abundance 正式文章
- `chapters/41-read-mapping-depth.qmd`：第 41 篇 Bowtie2/Samtools/JGI 完整坐标回比与深度矩阵正式文章
- `chapters/42-binning.qmd`：第 42 篇 MetaBAT2/SemiBin2/VAMB/TaxVAMB 多分支盲法分箱正式文章
- `chapters/43-bin-refinement.qmd`：第 43 篇 DAS Tool/Binette、独立 QC 与 truth-blind 选集正式文章
- `chapters/44-mag-qc-checkm2-gunc-mimag.qmd`：第 44 篇 CheckM2/GUNC、完整 MIMAG 与 assembly graph 正式文章
- `chapters/62-element-cycling.qmd`：第 62 篇 C/N/S 元素循环、环境关联与 MAG carrier 审计正式文章
- `R/theme_pub.R`：整仓库共享作图函数
- `data/small/README.md`：小数据来源与再生说明
- `scripts/audit_chapter_completion.py`：77 篇完成状态与 QA 新鲜度审计
- `scripts/build_wechat_review_bundle.py`：与 16S 系列同款的本地公众号审阅包派生器
- `scripts/render_wechat_preview.py`：兼容旧命令的入口，转交给上述生成器
- `scripts/validate_article11_installation.py`：双环境、入口证据与数据库门禁验证器
- `scripts/retrieve_article12_cmd.R`：第 12 篇一次性 ExperimentHub 真实取回与离线 cache replay
- `scripts/validate_article12_r_cmd.R`：第 12 篇不联网的包、lock、资源、对象与图片验收器
- `scripts/build_article13_fastq_subset.py`：第 13 篇同步 ENA FASTQ prefix builder
- `scripts/run_article13_read_qc.sh`：第 13 篇一次性 FastQC/fastp/MultiQC 运行器
- `scripts/validate_article13_read_qc.py`：第 13 篇离线冻结产物、账本与图片验收器
- `scripts/build_article14_fastq_controls.py`：第 14 篇双控制同步 FASTQ prefix builder
- `scripts/run_article14_host_removal.sh`：第 14 篇一次性 Hostile 与敏感性分支运行器
- `scripts/validate_article14_host_removal.py`：第 14 篇离线冻结产物、账本与图片验收器
- `scripts/run_article15_metaphlan4.sh`：第 15 篇一次性完整数据库 profiling 与深度分支运行器
- `scripts/validate_article15_metaphlan4.py`：第 15 篇离线冻结证据、marker 分账与图片验收器
- `scripts/run_article16_kraken2_bracken.sh`：第 16 篇一次性 Kraken2 分类与四组 Bracken 分支运行器
- `scripts/validate_article16_kraken2_bracken.py`：第 16 篇离线冻结证据、paired-fragment 守恒、重分配与图片验收器
- `scripts/run_article17_kraken2_database_confidence.sh`：第 17 篇 restart-safe 三数据库、confidence 与 hit-group 一次性运行器
- `scripts/validate_article17_kraken2_database_confidence.py`：第 17 篇 reference-aware 真值、human-lineage 分区、空 Bracken、资源与图片离线验收器
- `scripts/run_article18_profiler_benchmark.sh`：第 18 篇 mOTUs 一次映射、多 `g` 分支、crosswalk 与资源账本一次性运行器
- `scripts/validate_article18_profiler_benchmark.py`：第 18 篇三套 feature space、共同域、组成一致性、资源与图片离线验收器
- `scripts/bootstrap_article19_humann_databases.sh`：第 19 篇四归档分层下载、校验、解包与 inventory 入口
- `scripts/download_ranged_archive.sh`：支持 exact range coverage 与重试的显式并行归档下载器
- `scripts/run_article19_humann3.sh`：第 19 篇兼容 prescreen、HUMAnN 两级检索、reaction regroup、归一化分支与冻结证据一次性运行器
- `scripts/validate_article19_humann3.py`：第 19 篇 read ledger、单位、regroup/strata、数据库、资源与四张图片的离线验收器
- `scripts/prepare_article20_cmd_pathways.R`：第 20 篇真实队列 pathway abundance/coverage 资源取回、对齐与冻结器
- `scripts/run_article20_humann_normalization.sh`：第 20 篇四个真实 HUMAnN pathway 归一化分支运行器
- `scripts/freeze_article20_humann_normalization.py`：第 20 篇运行证据、日志、版本与 checksum 冻结器
- `scripts/validate_article20_functional_normalization.py`：第 20 篇 24 分支、special rows、coverage/prevalence、零值与四张图片离线验收器
- `scripts/run_article21_microbecensus.sh`：第 21 篇固定源码、同一 clean sequence universe 与双 read-length MicrobeCensus 一次性运行器
- `scripts/run_article21_microbecensus.py`：第 21 篇 RAPsearch2 preflight 兼容层与 AGS/GE 运行入口
- `scripts/freeze_article21_microbecensus.py`：第 21 篇版本、资源、日志、命令与 checksum 冻结器
- `scripts/validate_article21_table_semantics.py`：第 21 篇单位、分母、closure、zero、strata、RPKG 与四张图片离线验收器
- `data/small/21-table-semantics-frozen/`：第 21 篇九项 checksum-locked MicrobeCensus 小产物
- `results/21-table-semantics/`：第 21 篇语义合同、变换裁决与验证结果
- `scripts/prepare_article22_diversity_data.R`：第 22 篇 cMD 与 Figshare 真实表冻结器
- `scripts/validate_article22_diversity.R`：第 22 篇多分辨率 alpha/beta、PCoA、敏感性与四张图片离线验收器
- `data/small/22-diversity-inputs/`：第 22 篇九项 checksum-locked 物种、基因与 MAG 冻结输入
- `results/22-alpha-beta-diversity/`：第 22 篇多样性、排序、重复测量、解释边界与验证结果
- `scripts/prepare_article23_ordination_data.R`：第 23 篇 sample-to-spring 聚合、合同与 checksum 冻结器
- `scripts/validate_article23_ordination_permanova.R`：第 23 篇 PCoA/CAP、受限 PERMANOVA、PERMDISP、pairwise 与敏感性验收器
- `data/small/23-ordination-permanova/`：第 23 篇七项 checksum-locked spring-level 冻结输入
- `results/23-ordination-permanova/`：第 23 篇排序、置换空间、location/dispersion、pairwise 与验证结果
- `scripts/prepare_article28_cross_cohort.R`：第 28 篇八队列 cMD 原始资源重建器
- `scripts/validate_article28_cross_cohort.R`：第 28 篇 meta-analysis、LODO、bootstrap、迁移与图片验收器
- `data/small/28-cross-cohort/`：第 28 篇六项 checksum-locked 冻结输入
- `results/28-cross-cohort/`：第 28 篇逐队列效应、预测、调参、异质性、稳定性与验证结果
- `scripts/prepare_article29_network_data.R`：第 29 篇 spring-level MAG 输入、筛选与分析合同冻结器
- `scripts/validate_article29_network.R`：第 29 篇 graphical lasso、StARS、bootstrap、Zi-Pi、零模型与图片验收器
- `data/small/29-network/`：第 29 篇五项 checksum-locked 网络输入
- `results/29-network/`：第 29 篇边、节点角色、稳定性、敏感性、零模型与删除鲁棒性结果
- `scripts/select_article30_read_pairs.py`：第 30 篇一次遍历、无放回、同步 read-pair 精确抽样器
- `scripts/download_article30_assembly_reads.sh`：第 30 篇 ENA bytes/MD5 门禁与 restart-safe 下载入口
- `scripts/run_article30_short_read_assembly.sh`：第 30 篇 QC、六分支组装、索引、八次 mapping 与资源记录运行器
- `scripts/freeze_article30_short_read_assembly.py`：第 30 篇可移植日志、版本、表格、contigs 与 checksum 冻结器
- `scripts/validate_article30_short_read_assembly.py`：第 30 篇离线重算、解释边界与四张图片验收器
- `data/small/30-short-read-assembly-frozen/`：第 30 篇 66 项 checksum-locked 固化证据
- `results/30-short-read-assembly/`：第 30 篇来源、抽样、QC、组装、回帖、资源与验证审计结果
- `scripts/download_article31_long_reads.sh`：第 31 篇完整 ONT/HiFi ENA archives 的 restart-safe bytes/MD5 下载入口
- `scripts/bootstrap_article31_hifiasm_meta.sh`：第 31 篇 hifiasm-meta tag/commit/source checksum 锁定编译入口
- `scripts/run_article31_long_read_assembly.sh`：第 31 篇四分支组装、base-level read-back、junction 与资源记录运行器
- `scripts/freeze_article31_long_read_assembly.py`：第 31 篇可移植日志、表格、四套 contigs 与 checksum 冻结器
- `scripts/validate_article31_long_read_assembly.py`：第 31 篇离线重算、解释边界与四张图片验收器
- `data/small/31-long-read-assembly-frozen/`：第 31 篇 68 项 checksum-locked 固化证据
- `results/31-long-read-assembly/`：第 31 篇来源、工具、组装、read-back、junction、资源与验证审计结果
- `scripts/download_article32_hybrid_sources.sh`：第 32 篇四个 ENA archives、补充表与作者仓库的 restart-safe 下载和身份门禁
- `scripts/select_article32_read_pairs.py`：第 32 篇一次遍历、无放回、同步 10,000,000 read-pair 精确抽样器
- `scripts/run_article32_hybrid_assembly.sh`：第 32 篇七分支组装、BWA/Polypolish、MetaQUAST 与资源账本运行器
- `scripts/freeze_article32_hybrid_assembly.py`：第 32 篇参考、assemblies、日志、表格与 checksum 固化器
- `scripts/validate_article32_hybrid_assembly.py`：第 32 篇离线重算、解释边界与四张图片验收器
- `data/small/32-hybrid-assembly-polishing-frozen/`：第 32 篇 209 项 checksum-locked 固化证据
- `results/32-hybrid-assembly-polishing/`：第 32 篇来源、truth、七分支、polishing、MetaQUAST、资源与验证审计结果
- `scripts/download_article33_qc_sources.sh`：第 33 篇作者 benchmark repository 的 commit 与 truth inventory 下载门禁
- `scripts/prepare_article33_qc_inputs.py`：第 33 篇 15 套 assembly 身份、71/87 truth 和两个确定性正控准备器
- `scripts/run_article33_assembly_qc.sh`：第 33 篇一次 reference-free QUAST、三次 MetaQUAST 与资源账本运行器
- `scripts/summarize_article33_assembly_qc.py`：第 33 篇 17 分支、1,271 行逐 reference 与正控效应汇总器
- `scripts/freeze_article33_assembly_qc.py`：第 33 篇小表、报告、日志、资源和脚本固化器
- `scripts/validate_article33_assembly_qc.py`：第 33 篇离线身份、指标、正控、正文与四张图片验收器
- `data/small/33-assembly-qc-frozen/`：第 33 篇 45 项 checksum-locked 固化证据
- `results/33-assembly-qc/`：第 33 篇来源、指标、正控、正文、图片与验证审计结果
- `scripts/download_article34_gene_catalog_sources.sh`：第 34 篇论文 XML、benchmark commit 与 truth genomes 下载门禁
- `scripts/prepare_article34_gene_catalog_inputs.py`：第 34 篇六套 assembly 身份、87-genome truth 与 ORF 输入准备器
- `scripts/run_article34_gene_catalog.py`：第 34 篇 Prodigal、MMseqs2/CD-HIT、多策略与阈值分支运行器
- `scripts/summarize_article34_gene_catalog.py`：第 34 篇两阶段 membership、真值与资源汇总器
- `scripts/freeze_article34_gene_catalog.py`：第 34 篇主目录、表格、日志、环境、脚本与 checksum 固化器
- `scripts/validate_article34_gene_catalog.py`：第 34 篇离线身份、守恒、真值、边界、正文与四张图片验收器
- `data/small/34-nonredundant-gene-catalog-frozen/`：第 34 篇 44 项 checksum-locked 固化证据
- `results/34-nonredundant-gene-catalog/`：第 34 篇来源、目录、真值、正文、图片与验证审计结果
- `scripts/download_article35_gene_abundance_sources.sh`：第 35 篇 ENA reads、paired fastp 与 catalog source 门禁
- `scripts/prepare_article35_gene_abundance_inputs.py`：第 35 篇 reads、catalog、UniRef90、reaction mapping 与论文 XML 身份准备器
- `scripts/parse_article35_sam.py`：第 35 篇 streaming SAM、CIGAR identity/query coverage 与四策略计数器
- `scripts/run_article35_gene_abundance.py`：第 35 篇论文兼容、全量 mapping 与 iterative-sensitive DIAMOND 运行器
- `scripts/summarize_article35_gene_abundance.py`：第 35 篇 units、functional crosswalk、mass ledger 与资源汇总器
- `scripts/freeze_article35_gene_abundance.py`：第 35 篇小表、日志、环境、脚本与 checksum 固化器
- `scripts/validate_article35_gene_abundance.py`：第 35 篇离线身份、守恒、单位、功能、正文与四张图片验收器
- `data/small/35-gene-abundance-frozen/`：第 35 篇 checksum-locked gene/reaction abundance 证据
- `results/35-gene-abundance/`：第 35 篇 mapping、normalization、function、正文、图片与验证审计结果
- `scripts/bootstrap_article36_eggnog_database.sh`：第 36 篇 eggNOG 5.0.2 archives 下载、installed files 与 SQLite release 门禁
- `scripts/prepare_article36_eggnog_inputs.py`：第 36 篇 catalog、abundance、environment 与四项 database asset 身份准备器
- `scripts/run_article36_eggnog_annotation.py`：第 36 篇 non-electronic 主注释与同 seed all-evidence GO 敏感性一次性运行器
- `scripts/summarize_article36_eggnog_annotation.py`：第 36 篇字段覆盖、四态、ORF strata、COG/KO/GO fractional ledger 与资源汇总器
- `scripts/freeze_article36_eggnog_annotation.py`：第 36 篇 annotations、seed orthologs、日志、环境、脚本与 checksum 固化器
- `scripts/validate_article36_eggnog_annotation.py`：第 36 篇离线身份、守恒、GO policy、正文与四张图片验收器
- `data/small/36-eggnog-functional-annotation-frozen/`：第 36 篇 45 项 checksum-locked annotation 证据
- `results/36-eggnog-functional-annotation/`：第 36 篇 database、annotation、正文、图片与验证审计结果
- `scripts/prepare_article37_cazymes.py`、`run_article37_cazymes.py`、`summarize_article37_cazymes.py`、`freeze_article37_cazymes.py`：第 37 篇输入、dbCAN 实跑、汇总与固化链
- `scripts/prepare_article38_resistome.py`、`run_article38_resistome.py`、`summarize_article38_resistome.py`、`freeze_article38_resistome.py`：第 38 篇 CARD-RGI 实跑与固化链
- `scripts/prepare_article39_virulome.py`、`run_article39_virulome.py`、`summarize_article39_virulome.py`、`freeze_article39_virulome.py`：第 39 篇 VFDB-ABRicate 实跑与固化链
- `scripts/prepare_article40_bgc.py`、`run_article40_bgc.py`、`summarize_article40_bgc.py`、`freeze_article40_bgc.py`：第 40 篇 antiSMASH/GECCO 实跑与固化链
- `scripts/validate_article37_40.py`：第 37–40 篇共享 checksum、科学合同、正文与三格式图片验收器；四个 article-specific wrappers 是正式入口
- `data/small/37-cazymes-dbcan-frozen/` 至 `data/small/40-bgc-natural-products-frozen/`：第 37–40 篇 checksum-locked 小表、日志、资源、环境与脚本证据
- `results/37-cazymes-dbcan/` 至 `results/40-bgc-natural-products/`：第 37–40 篇离线校验与发表图审计结果
- `scripts/prepare_article41_mapping_depth.py`、`run_article41_mapping_depth.py`、`summarize_article41_mapping_depth.py`、`freeze_article41_mapping_depth.py`：第 41 篇完整坐标 mapping/depth 实跑与固化链
- `scripts/prepare_article42_binning.py`、`run_article42_binning.py`、`summarize_article42_binning.py`、`freeze_article42_binning.py`：第 42 篇五分支分箱、独立 QC 与 truth audit 固化链
- `scripts/prepare_article43_refinement.py`、`run_article43_refinement.py`、`summarize_article43_refinement.py`、`freeze_article43_refinement.py`：第 43 篇 DAS Tool/Binette 精炼与 truth-blind 选集固化链
- `scripts/prepare_article44_mag_qc.py`、`run_article44_mag_qc.py`、`summarize_article44_mag_qc.py`、`freeze_article44_mag_qc.py`：第 44 篇完整 MIMAG、domain mode 与 assembly-graph 审计固化链
- `scripts/validate_article41_mapping_depth.py` 至 `scripts/validate_article44_mag_qc.py`：第 41–44 篇 checksum、科学合同、正文与三格式图片的离线验收入口
- `data/small/41-read-mapping-depth-frozen/` 至 `data/small/44-mag-qc-mimag-graph-frozen/`：第 41–44 篇 checksum-locked 深度、binning、refinement 与 MAG-QC 证据
- `results/41-read-mapping-depth/` 至 `results/44-mag-qc-mimag-graph/`：第 41–44 篇离线验证与发表图审计结果
- `scripts/download_article62_element_data.py`、`run_article62_element_cycles.py`、`plot_article62_element_cycles.py`、`freeze_article62_element_cycles.py`：第 62 篇真实 hot-spring community/MAG 数据下载、分析、作图与固化链
- `scripts/validate_article62_element_cycles.py`：第 62 篇 checksum、DNF marker、置换/bootstrap、解释边界、正文与三格式图片验收器
- `data/small/62-element-cycling-frozen/`：第 62 篇 37 项 checksum-locked 输入、结果、环境与再生脚本证据
- `qa/article62/validation-summary.json`：第 62 篇 251/251 项离线验证摘要
- `qa_report.json`：本地发布门禁报告
- `rendered/wechat_review_01_77_16s_style/01/` 至
  `rendered/wechat_review_01_77_16s_style/77/`：
  与 16S 系列同款的 77 篇手机端本地审阅包

## 许可证

教程文字与原创图使用 CC BY 4.0；原创代码使用 MIT。第三方数据、软件和论文材料
保持各自许可证。
