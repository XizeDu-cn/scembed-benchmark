# Python 基于 scIB 的嵌入聚类基准评估（scATAC/scRNA 通用）

本工具针对 `adata.h5ad` 中的多个嵌入（存于 `obsm`）执行多分辨率、多随机种子聚类评估，并输出与示例类似的分层目录和汇总表格（TSV）。指标以 Python scIB 为主，必要时补充 sklearn。

## 安装依赖

```bash
pip install -r requirements.txt
```

## 运行示例

```bash
python run_benchmark.py \
  --adata /path/to/adata.h5ad \
  --label-key cell_type \
  --batch-key batch \
  --obsm-keys X_lsi X_pca X_umap X_tsne \
  --resolutions 0.05 0.1 0.2 0.5 1.0 \
  --seeds 0 2 5 \
  --n-neighbors 15 \
  --scenario Cell_line_mixing \
  --output-dir /path/to/output
```

- 若未指定 `--obsm-keys`，默认使用 `adata.obsm.keys()` 中的全部嵌入（建议恰好四个）。
- 若无 `--label-key`（真值标签），将跳过 ARI/AMI/NMI，用 Silhouette/CH/DB 等无监督指标选择最优分辨率。

## 输出结构

```
{output-dir}/{scenario}/
  ├─ evaluation/{method}/seed{seed}/r{resolution}/metrics.tsv
  ├─ clustering/{method}/seed{seed}/r{resolution}.tsv
  ├─ index.tsv                   # 逐次运行的文件与元数据映射
  ├─ summary.tsv                 # 所有 metrics.tsv 的长表汇总
  └─ best_per_method.tsv         # 每个 method 的最优分辨率/簇数（按选择指标）
```

- metrics.tsv（长表）包含：`n_clusters, metric, value, method, long_method, feature_type, ndim, resolution, seed, scenario, dataset, dataset2, clustering_file, snn_file, k_optimal` 等字段。
- index.tsv 映射每次运行的 `file → 元数据`，便于快速检索。

## 常见参数说明
- **--label-key**: 真值标签列名，位于 `adata.obs`。
- **--batch-key**: 批次列名（若提供，将额外评估 `graph_connectivity` 和 `silhouette_batch` 相关指标）。
- **--select-metric**: 选择最优分辨率的指标，默认 `ARI`；若无标签自动回退到 `SILHOUETTE_CLUSTER`。

## 备注
- 聚类使用 Scanpy Leiden（按每个嵌入和分辨率/种子分别建立图）。
- 主指标来自 scIB（ARI/NMI/Silhouette/Graph Connectivity 等），AMI/CH/DB 由 sklearn 计算。
