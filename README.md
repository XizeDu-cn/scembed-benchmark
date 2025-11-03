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

### 1. 准备数据

#### 数据要求

你需要准备一个AnnData对象（`.h5ad`格式），必须包含：

1. **嵌入矩阵**（存储在 `adata.obsm`）
   - 每个嵌入方法对应一个矩阵，形状为 `(n_cells, n_dims)`
   - 键名示例：`cellspace`, `simba`, `scbasset`, `gfetm` 等
   - 也支持带前缀的键名：`X_cellspace`, `X_SIMBA` 等

2. **细胞注释**（存储在 `adata.obs`）
   - `cell_type` - 真实细胞类型标签（**必需**）
   - `batch` - 批次信息（可选，用于批次效应评估）

#### 数据准备示例

```python
import scanpy as sc
import numpy as np

# 1. 读取你的数据
adata = sc.read_h5ad("your_raw_data.h5ad")

# 2. 运行不同的嵌入方法，并保存到obsm
# 示例：添加LSI嵌入（如果你有scATAC-seq数据）
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.tl.pca(adata, n_comps=50)
adata.obsm['lsi'] = adata.obsm['X_pca']  # 保存为lsi

# 示例：添加其他方法的嵌入
# adata.obsm['cellspace'] = your_cellspace_embedding
# adata.obsm['simba'] = your_simba_embedding

# 3. 检查数据结构
print("嵌入方法:", list(adata.obsm.keys()))
print("细胞数:", adata.n_obs)
print("标签列:", adata.obs.columns.tolist())

# 4. 确保标签正确
print("\n真实标签统计:")
print(adata.obs['cell_type'].value_counts())

# 5. 检查是否有缺失值
assert adata.obs['cell_type'].notna().all(), "cell_type列包含缺失值！"

# 6. 保存准备好的数据
adata.write_h5ad("adata_ready_for_benchmark.h5ad")
print("\n数据已保存，可以运行基准测试！")
```

#### 快速检查数据

运行以下命令检查你的数据是否符合要求：

```python
import scanpy as sc

adata = sc.read_h5ad("your_data.h5ad")

# 检查1: 嵌入是否存在
assert len(adata.obsm.keys()) > 0, "错误：obsm中没有嵌入！"
print(f"✓ 找到 {len(adata.obsm.keys())} 个嵌入方法")

# 检查2: 标签是否存在
assert 'cell_type' in adata.obs.columns, "错误：缺少cell_type列！"
print(f"✓ cell_type列存在，包含 {adata.obs['cell_type'].nunique()} 个类别")

# 检查3: 无缺失值
assert adata.obs['cell_type'].notna().all(), "错误：cell_type包含缺失值！"
print("✓ cell_type无缺失值")

# 检查4: 嵌入维度
for key in adata.obsm.keys():
    print(f"  - {key}: {adata.obsm[key].shape}")

print("\n✓ 数据检查通过！可以开始基准测试。")
```

### 2. 配置参数

编辑 `config.yaml` 文件：

```yaml
# 基本输入
scenario: Cell_line_mixing
adata: /path/to/your/adata.h5ad
label_key: cell_type
batch_key: batch  # 可选

# 要评估的嵌入方法（对应adata.obsm中的键名）
obsm_keys: [cellspace, simba, scbasset, gfetm]

# 运行参数
seeds: [0, 2, 5]           # 随机种子，建议3个
n_neighbors: 15            # KNN图的邻居数
select_metric: ARI         # 选择最优分辨率的指标
output_dir: ./results

# 分辨率配置（二选一）
# 方式1: 网格表达式（推荐）
res_grid: "0.02:2.00:0.08"  # 从0.02到2.00，步长0.08

# 方式2: 显式列表
# resolutions: [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
```

### 3. 运行工作流

```bash
# 激活环境
conda activate scib-benchmark

# 查看工作流DAG图
snakemake --dag | dot -Tpng > workflow_dag.png

# 干运行，查看将要执行的任务
snakemake -n

# 运行完整工作流（推荐添加 --scheduler greedy 优化任务调度）
snakemake --cores 8 --scheduler greedy

# 如果遇到调度问题，可以使用 ilp 调度器
snakemake --cores 8 --scheduler ilp
```

**重要提示**：
- 使用 `--scheduler greedy` 可以提高大规模任务的调度效率
- `--cores` 指定并行任务数，根据你的CPU核心数调整
- 首次运行建议先 `snakemake -n` 确认任务列表

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

如果你想单独运行基准测试脚本，不通过Snakemake：

### 基本用法

```bash
python run_benchmark.py \
  --adata data/adata.h5ad \
  --label-key cell_type \
  --batch-key batch \
  --obsm-keys cellspace simba scbasset gfetm \
  --res-grid "0.1:2.0:0.1" \
  --seeds 0 2 5 \
  --n-neighbors 15 \
  --scenario Cell_line_mixing \
  --output-dir results
```

**注意**：`obsm-keys` 参数不需要 `X_` 前缀，脚本会自动处理。

### 使用早停优化（推荐）

对于大规模参数扫描，使用早停策略可以显著加速：

```bash
python run_benchmark.py \
  --adata data/adata.h5ad \
  --label-key cell_type \
  --obsm-keys cellspace simba \
  --res-grid "0.02:2.0:0.02" \
  --seeds 0 2 5 \
  --infer-target-k-from-label \
  --early-stop \
  --k-margin 3 \
  --scenario Cell_line_mixing \
  --output-dir results
```

**早停参数说明**：
- `--infer-target-k-from-label`: 从真实标签自动推断目标簇数
- `--early-stop`: 启用早停，当簇数达到目标后停止扫描更高分辨率
- `--k-margin 3`: 在目标簇数基础上多扫描3个簇的冗余

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

## 高级用法

### Snakemake常用选项

```bash
# 只运行特定方法
snakemake --cores 4 results/Cell_line_mixing/best_per_method.tsv

# 重新运行失败的任务
snakemake --cores 8 --rerun-incomplete

# 生成工作流报告
snakemake --report report.html

# 清理所有输出（谨慎使用！）
snakemake --delete-all-output

# 使用Conda环境管理（如果env.yaml配置了环境）
snakemake --cores 8 --use-conda
```

### 任务调度优化

对于大规模任务（数千个任务），推荐使用优化的调度器：

```bash
# greedy调度器（快速，适合大多数情况）
snakemake --cores 16 --scheduler greedy

# ilp调度器（最优化，但较慢）
snakemake --cores 16 --scheduler ilp
```

### 部分运行特定方法

如果只想测试某几个方法，临时修改 `config.yaml` 中的 `obsm_keys`：

```yaml
obsm_keys: [cellspace, simba]  # 只运行这两个方法
```

## 最佳实践

1. **数据准备**
   - 确保 `adata.obsm` 中的键名简洁清晰（如 `cellspace` 而非 `X_cellspace_embedding_final`）
   - 真实标签列 `cell_type` 应该是分类型变量，无缺失值
   - 如果有批次效应，确保 `batch` 列正确标注

2. **参数设置**
   - 使用至少**3个随机种子**（0, 2, 5）确保结果稳定性
   - 分辨率范围应**覆盖真实类别数**附近（建议 `0.1` 到 `2.0`）
   - KNN邻居数建议 **10-30** 之间，默认15

3. **运行策略**
   - 首次运行先用 `snakemake -n` 检查任务数量
   - 启用 `--scheduler greedy` 优化大规模任务调度
   - 使用 `--cores` 根据CPU核心数合理设置并行度
   - 对于探索性分析，先用少量分辨率和种子测试

4. **结果解读**
   - 查看 `best_per_method.tsv` 了解每个方法的最优性能
   - 检查 `summary.tsv` 分析不同参数对结果的影响
   - 使用混淆矩阵图直观评估聚类质量

## 常见问题

### Q1: 为什么要使用Snakemake而不是直接运行Python脚本？

Snakemake提供了：
- **自动并行化**：同时运行多个独立任务
- **断点续传**：失败后只重新运行未完成的任务
- **依赖管理**：自动处理文件依赖关系
- **可重复性**：完整记录工作流程

### Q2: 如何选择合适的分辨率范围？

建议使用网格表达式 `res_grid: "0.02:2.00:0.08"`，这样：
- 低分辨率（0.02-0.5）：适合粗粒度聚类
- 中分辨率（0.5-1.0）：通常接近真实类别数
- 高分辨率（1.0-2.0）：用于细粒度分群

### Q3: obsm中的嵌入键名格式要求？

Snakefile会自动处理嵌入键名。在 `config.yaml` 中使用简化名称：
- 配置文件: `obsm_keys: [cellspace, simba]`
- 实际键名可以是: `X_cellspace`, `cellspace`, `X_CELLSPACE` 等

### Q4: 如何处理大型数据集？

- 使用 `--early-stop` 减少不必要的分辨率扫描
- 增加 `--cores` 提高并行度
- 考虑先在数据子集上测试参数

## 输出文件说明

### metrics.tsv
包含所有评估指标的长表格式，字段包括：
- `method`: 嵌入方法名
- `resolution`: Leiden分辨率
- `seed`: 随机种子
- `metric`: 指标名称（ARI, NMI, Silhouette等）
- `value`: 指标值
- `n_clusters`: 聚类数量

### best_per_method.tsv
每个方法的最优配置（按select_metric选择）：
- 最优分辨率
- 对应的簇数
- 所有指标的最优值

### confusion_matrix.png
混淆矩阵可视化，包含：
- 主热图：真实标签 vs 预测簇
- 同质性条形图：每个簇的纯度
- 完整性条形图：每个类别的集中度
- ARI指标显示

## 引用

如果本工具对你的研究有帮助，请引用相关的方法论文。

## 许可证

MIT License

## 联系方式

- GitHub: https://github.com/XizeDu-cn/scembed-benchmark
- 问题反馈：请在GitHub仓库提交Issue
