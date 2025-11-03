# scATAC-seq嵌入基准测试

基于Snakemake的scATAC-seq嵌入方法评估工作流，用于系统性比较不同降维方法的聚类性能。

## 项目概述

本项目使用Snakemake工作流管理系统，自动化执行以下任务：
- 对多种嵌入方法进行多分辨率Leiden聚类
- 使用多个随机种子确保结果稳定性
- 计算多种聚类评估指标（ARI、NMI、Silhouette等）
- 生成混淆矩阵可视化图表
- 汇总并比较不同方法的性能

## 环境配置

### 使用Conda创建环境（推荐）

```bash
conda env create -f env.yaml
conda activate scib-benchmark
```

或者手动创建环境：

```bash
conda create -n scib-benchmark -c conda-forge -c bioconda \
  python=3.12 \
  snakemake \
  scanpy pandas numpy scipy \
  r-base r-ggplot2 r-dplyr r-ggpubr r-patchwork r-argparse
```

### 使用pip安装Python依赖

```bash
pip install -r requirements.txt
```

## 快速开始

### 1. 配置参数

编辑 `config.yaml` 文件：

```yaml
adata_path: "/path/to/your/adata.h5ad"
label_key: "cell_type"
batch_key: "batch"
scenario: "Cell_line_mixing"
```

### 2. 运行工作流

```bash
# 查看工作流DAG
snakemake --dag | dot -Tpng > dag.png

# 干运行，查看将要执行的任务
snakemake -n

# 运行完整工作流
snakemake --cores 8
```

## 工作流说明

### 主要步骤

1. **数据准备** - 提取真实标签
2. **聚类分析** - 对每个嵌入方法使用不同分辨率和随机种子进行Leiden聚类
3. **评估计算** - 计算聚类质量指标（ARI、NMI、Silhouette等）
4. **可视化** - 生成混淆矩阵热图
5. **结果汇总** - 聚合所有方法的最佳性能

### 输出结构

```
results/
└── {scenario}/
    ├── clustering/
    │   └── {method}/
    │       └── seed{seed}/
    │           └── r{resolution}.tsv
    ├── evaluation/
    │   └── {method}/
    │       └── seed{seed}/
    │           └── r{resolution}/
    │               ├── metrics.tsv
    │               └── confusion_matrix.png
    ├── ground_truth.tsv
    ├── index.tsv
    ├── summary.tsv
    └── best_per_method.tsv
```

## 直接运行脚本（不使用Snakemake）

如果你想单独运行基准测试：

```bash
python run_benchmark.py \
  --adata data/adata.h5ad \
  --label-key cell_type \
  --batch-key batch \
  --obsm-keys X_cellspace X_simba X_scbasset X_gfetm \
  --resolutions 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 \
  --seeds 0 2 5 \
  --n-neighbors 15 \
  --scenario Cell_line_mixing \
  --output-dir results
```

## 混淆矩阵可视化

单独生成混淆矩阵图：

```bash
Rscript plot_confusion_matrix.R \
  --clusters results/Cell_line_mixing/clustering/simba/seed0/r0.5.tsv \
  --truth results/Cell_line_mixing/ground_truth.tsv \
  --output confusion_matrix.png
```

详细说明见 [README_plot_confusion_matrix.md](README_plot_confusion_matrix.md)

## 评估指标

本工作流计算以下聚类评估指标：

### 监督指标（需要真实标签）
- **ARI** (Adjusted Rand Index) - 聚类准确性
- **NMI** (Normalized Mutual Information) - 信息一致性
- **AMI** (Adjusted Mutual Information) - 调整后互信息

### 无监督指标
- **Silhouette Score** - 簇内紧密度
- **Calinski-Harabasz Index** - 簇间分离度
- **Davies-Bouldin Index** - 簇质量

### 批次效应指标（需要batch_key）
- **Graph Connectivity** - 批次混合程度
- **Batch Silhouette** - 批次校正效果

## 项目结构

```
.
├── Snakefile                           # Snakemake工作流定义
├── config.yaml                         # 配置文件
├── env.yaml                            # Conda环境配置
├── requirements.txt                    # Python依赖
├── run_benchmark.py                    # 主基准测试脚本
├── aggregate_results.py                # 结果聚合脚本
├── extract_truth_labels.py             # 标签提取工具
├── plot_confusion_matrix.R             # 混淆矩阵绘图脚本
├── README.md                           # 本文件
├── README_plot_confusion_matrix.md     # R脚本说明
└── results/                            # 输出目录（被gitignore）
```

## 参数说明

### 聚类参数
- `--resolutions`: 聚类分辨率列表，值越大簇越多
- `--seeds`: 随机种子列表，用于确保结果可重复性
- `--n-neighbors`: KNN图的邻居数，默认15

### 数据参数
- `--label-key`: AnnData对象中真实标签的列名（`adata.obs`）
- `--batch-key`: 批次信息列名（可选）
- `--obsm-keys`: 要评估的嵌入方法列表

### 输出参数
- `--scenario`: 场景名称，用于组织输出目录
- `--output-dir`: 结果输出根目录
- `--select-metric`: 选择最优分辨率的指标，默认ARI

## 最佳实践

1. **使用多个随机种子**（如0, 2, 5）确保结果稳定性
2. **设置合理的分辨率范围**，覆盖真实类别数附近的值
3. **查看DAG图**了解工作流执行计划
4. **使用Snakemake的并行功能**加速计算：`snakemake --cores 8`

## 引用

如果本工具对你的研究有帮助，请引用相关的方法论文。

## 许可证

MIT License

## 联系方式

如有问题或建议，请在GitHub仓库提issue。
