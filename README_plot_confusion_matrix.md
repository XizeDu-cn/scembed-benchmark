# 混淆矩阵绘图脚本使用说明

## 概述

`plot_confusion_matrix.R` 是一个命令行R脚本，用于生成混淆矩阵热图，比较真实标签与聚类结果。

## 依赖要求

- R环境（推荐使用scib-benchmark conda环境）
- 必需的R包：
  - argparse
  - ggplot2
  - dplyr
  - ggpubr

## 使用方法

```bash
Rscript plot_confusion_matrix.R --clusters <聚类文件> --truth <真实标签文件> --output <输出文件>
```

### 参数说明

- `--clusters`: 输入聚类结果TSV文件路径
- `--truth`: 输入真实标签TSV文件路径  
- `--output`: 输出PNG图片文件路径

### 示例

```bash
# 使用scib-benchmark环境运行
/opt/anaconda3/envs/scib-benchmark/bin/Rscript plot_confusion_matrix.R \
  --clusters results/Cell_line_mixing/clustering/scbasset/seed0/r0.5.tsv \
  --truth truth_labels.tsv \
  --output confusion_matrix_plot.png
```

## 输入文件格式

### 聚类结果文件
TSV格式，包含以下列：
- `cell`: 细胞ID
- `cluster`: 聚类标签

### 真实标签文件  
TSV格式，包含以下列：
- `cell`: 细胞ID
- `truth_label`: 真实标签

## 输出

脚本会生成一个PNG格式的混淆矩阵热图，包含：
- 主热图：显示真实标签与聚类结果的对应关系
- 侧边条形图：显示同质性和完整性指标
- ARI和ARI2指标显示

## 注意事项

1. 脚本需要用户提供 `adjusted_wallance_indices` 函数的定义（在脚本中的注释占位符位置）
2. 确保输入文件格式正确，且都包含 `cell` 列用于数据合并
3. 输出图片尺寸为 10x12 英寸，DPI为300

## 集成到Snakemake工作流

可以在Snakemake规则中使用此脚本：

```python
rule plot_confusion_matrix:
    input:
        clusters = "results/{scenario}/clustering/{method}/seed{seed}/r{res}.tsv",
        truth = "truth_labels.tsv"
    output:
        plot = "results/{scenario}/plots/{method}_seed{seed}_r{res}_confusion.png"
    shell:
        "/opt/anaconda3/envs/scib-benchmark/bin/Rscript plot_confusion_matrix.R "
        "--clusters {input.clusters} "
        "--truth {input.truth} "
        "--output {output.plot}"
```


