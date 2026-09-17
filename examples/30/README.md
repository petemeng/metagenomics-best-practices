# 短读长组装：MEGAHIT vs metaSPAdes、单样本 vs 混合组装

数据与代码清单：[files.tsv](files.tsv)。约 115.1 MiB；不含原始 FASTQ 或大型参考数据库。

清单中 `input` 是用于本例的公开输入或有来源记录的处理后表格；`saved_result` 是已经计算的结果，便于核对图表，不代表重新拟合模型；`implementation` 是分析脚本；`environment` 是环境记录。

下载后按[本章完整说明](../../chapters/30-short-read-assembly.qmd)运行。原始数据的研究出处、样本筛选、统计单位和局限见正文。重计算流程还需要正文指定的软件、原始输入和数据库，下载这份清单不等于完成上游重跑。

在空文件夹中下载并运行 `download.R`，文件会保留正文所用的相对目录。已存在且与清单不一致的文件会报错，不会覆盖。
