# Bioinformatics NGS Pipeline

## 项目概述
这个仓库包含了一系列用于下一代测序(NGS)数据分析的生物信息学流程脚本，是我在研究生期间自行总结并在课题组的数据分析中经过长期实战验证的pipeline。这些脚本涵盖了从原始数据处理到下游分析的多个方面，包括染色质免疫沉淀-测序(ChIP-seq)、双亚硫酸盐测序(Bisulfite-seq)、RNA测序(RNA-seq)和转录起始位点(TSS)富集分析等。

## 包含的流程

### TSS 富集热图分析
`TSS_enrich_heatMap` 目录包含用于生成转录起始位点(TSS)和转录终止位点(TES)周围富集水平热图的脚本。这对于评估ChIP-seq或ATAC-seq数据质量和可视化启动子区域的蛋白质结合模式非常有用。

### 双亚硫酸盐测序分析
`bisulfite-seq` 目录提供了处理和分析双亚硫酸盐测序数据的脚本，用于DNA甲基化研究。

### ChIP-seq组蛋白分析
`chip-seq` 目录中的脚本专门用于处理ChIP-seq数据，特别关注组蛋白修饰的分析。

### tDNA插入重测序
`resequencing_tDNA_insert` 目录包含用于分析转基因DNA(tDNA)插入位点的脚本，利用双末端(PE)测序数据。

### RNA-seq分析
`rna-seq` 目录提供了RNA测序数据分析的完整流程，用于基因表达研究。

## 使用方法

### 环境要求
- Bash shell环境(我使用Linux系统，其余系统未尝试)
- 生物信息学工具：流程中提到的软件皆可用bioconda下载
- R 统计环境(用于数据可视化)
- Python 3

### 如何使用
仓库中的脚本分为两种类型：
1. **可执行程序**：这些脚本可以直接运行，通过自定义并调用函数进行快速大规模raw_data处理
2. **代码笔记**：类似于工作记录的代码片段集合，方便复制粘贴使用

请查看各个脚本的注释和文件头部说明以了解具体用法。大部分脚本都有详细的参数说明和使用示例。

## 文件说明

### 核心功能

各目录下的主要脚本功能简述：

- **TSS_enrich_heatMap/**
  - 生成基因TSS/TES区域±x kb的信号热图

- **bisulfite-seq/**
  - 从原始FASTQ文件到甲基化位点分析的完整流程

- **chip-seq/**
  - 组蛋白修饰ChIP-seq分析专用脚本

- **resequencing_tDNA_insert/**
  - 转基因插入位点精确定位和鉴定

- **rna-seq/**
  - 从测序数据到差异表达分析的完整流程

## 常见问题

- **Q: 脚本需要修改哪些参数？**  
  A: 大多数脚本需要修改输入/输出路径和参考基因组路径，请查看脚本开头的变量设置部分。

- **Q: 如何安装依赖工具？**  
  A: 推荐使用conda/bioconda安装，例如：`conda install -c bioconda bwa samtools bedtools`

## 贡献与反馈

欢迎通过Issues提出问题或建议。如果您对改进这些脚本有兴趣，也欢迎提交Pull Request。

## 许可证

本项目采用MIT许可证 - 查看 [LICENSE](LICENSE) 文件了解详情

## 引用
如果您在研究中使用了这些脚本，请考虑引用本仓库：FROST-web3. (2025). Bioinformatics_NGS_pipeline: A collection of validated NGS analysis pipelines. GitHub.
