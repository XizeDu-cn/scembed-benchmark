#!/usr/bin/env python3

import argparse
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

import scanpy as sc
from sklearn.metrics import (
    adjusted_mutual_info_score,
    adjusted_rand_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    normalized_mutual_info_score,
    silhouette_score,
)

try:
    # 尝试导入scIB指标计算模块
    from scib.metrics import ari as scib_ari
    from scib.metrics import nmi as scib_nmi
    from scib.metrics import silhouette as scib_silhouette
    from scib.metrics import graph_connectivity as scib_graph_connectivity
except Exception:
    # 导入失败时设为None
    scib_ari = None
    scib_nmi = None
    scib_silhouette = None
    scib_graph_connectivity = None


def parse_args() -> argparse.Namespace:
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="基于 scIB 对 AnnData 的 obsm 嵌入进行多分辨率/多种子聚类评估，并输出汇总表格"
    )
    parser.add_argument("--adata", required=True, help="输入 AnnData .h5ad 文件路径")
    parser.add_argument(
        "--obsm-keys",
        nargs="+",
        default=None,
        help="要评估的 obsm 键名（默认使用 adata.obsm.keys() 全部）",
    )
    parser.add_argument(
        "--label-key",
        default=None,
        help="真值标签列（位于 adata.obs），用于计算 ARI/NMI/AMI 等有监督指标",
    )
    parser.add_argument(
        "--batch-key",
        default=None,
        help="批次列（位于 adata.obs），若提供，将评估图连通性等批次相关指标",
    )
    parser.add_argument(
        "--resolutions",
        nargs="+",
        type=float,
        default=[0.05, 0.1, 0.2, 0.5, 1.0],
        help="Leiden 分辨率列表",
    )
    parser.add_argument(
        "--res-grid",
        type=str,
        default=None,
        help="以 start:end:step 形式给出分辨率网格（提供则覆盖 --resolutions）",
    )
    parser.add_argument(
        "--seeds",
        nargs="+",
        type=int,
        default=[0, 2, 5],
        help="随机种子列表",
    )
    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=15,
        help="构图时的 k（近邻数）",
    )
    parser.add_argument(
        "--scenario",
        type=str,
        default=None,
        help="场景（用于输出目录名）；默认取文件名不含扩展名",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./py_benchmark_results",
        help="输出根目录",
    )
    parser.add_argument(
        "--select-metric",
        type=str,
        default="ARI",
        choices=["ARI", "NMI", "AMI", "SILHOUETTE_CLUSTER", "CH", "DB"],
        help="选择最优分辨率的指标（若无标签，将自动回退到无监督指标）",
    )
    parser.add_argument(
        "--single-combo",
        action="store_true",
        help="仅运行单一 (method, seed, resolution) 组合并写出对应 metrics/clustering；不生成全局 summary/index/best",
    )
    # 目标簇数与早停策略
    parser.add_argument(
        "--target-k",
        type=int,
        default=None,
        help="目标簇数上界（期望覆盖到 1..target_k+k_margin）；与 --early-stop 配合可提前停止扫描",
    )
    parser.add_argument(
        "--k-margin",
        type=int,
        default=3,
        help="在目标簇数基础上的冗余覆盖量",
    )
    parser.add_argument(
        "--early-stop",
        action="store_true",
        help="开启早停：当观察到的最大簇数 ≥ target_k + k_margin 时，提前结束该 method/seed 的分辨率扫描",
    )
    parser.add_argument(
        "--infer-target-k-from-label",
        action="store_true",
        help="若提供 label_key，则从标签唯一值个数推断 target_k（当 --target-k 未显式给出时）",
    )
    return parser.parse_args()


def ensure_dir(path: str) -> None:
    """确保目录存在，如不存在则创建"""
    os.makedirs(path, exist_ok=True)


def load_adata(path: str):
    """加载AnnData文件"""
    adata = sc.read_h5ad(path)
    return adata


def list_obsm_keys(adata, override: Optional[List[str]]) -> List[str]:
    """列出要处理的obsm键名列表"""
    if override is not None and len(override) > 0:
        for key in override:
            if key not in adata.obsm.keys():
                raise ValueError(f"指定的 obsm 键不存在: {key}")
        return override
    return list(adata.obsm.keys())


def set_neighbors(adata, rep_key: str, n_neighbors: int, random_state: int, metric: str = "euclidean") -> None:
    """使用指定嵌入构建KNN图"""
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        use_rep=rep_key,
        metric=metric,
    )


def leiden_cluster(adata, resolution: float, seed: int, key_added: str) -> None:
    """执行Leiden聚类"""
    sc.tl.leiden(
        adata,
        resolution=resolution,
        random_state=seed,
        key_added=key_added,
    )


def compute_metrics(
    adata,
    embed_key: str,
    cluster_key: str,
    label_key: Optional[str],
    batch_key: Optional[str],
) -> Dict[str, float]:
    """计算各种评估指标
    
    包括:
    - 基础指标: 聚类数
    - 无监督指标: 轮廓系数、CH分数、DB分数
    - 有监督指标(需要label_key): ARI、NMI、AMI等
    - 批次相关指标(需要batch_key): 图连通性、批次ASW
    """
    X = adata.obsm[embed_key]
    labels_cluster = adata.obs[cluster_key].astype(str).values

    results: Dict[str, float] = {}

    # 基础聚类个数
    n_clusters = len(pd.unique(labels_cluster))
    results["N_CLUSTERS"] = float(n_clusters)

    # 无监督指标：轮廓系数/CH/DB
    try:
        results["SILHOUETTE_CLUSTER"] = float(
            silhouette_score(X, labels_cluster, metric="euclidean")
        )
    except Exception:
        results["SILHOUETTE_CLUSTER"] = np.nan

    try:
        results["CH"] = float(calinski_harabasz_score(X, labels_cluster))
    except Exception:
        results["CH"] = np.nan

    try:
        results["DB"] = float(davies_bouldin_score(X, labels_cluster))
    except Exception:
        results["DB"] = np.nan

    # 若存在真值标签，计算有监督指标
    if label_key is not None and label_key in adata.obs.columns:
        labels_true = adata.obs[label_key].astype(str).values
        # scIB 的 ARI/NMI
        if scib_ari is not None:
            try:
                results["ARI"] = float(scib_ari(adata, label_key, cluster_key))
            except Exception:
                results["ARI"] = float(adjusted_rand_score(labels_true, labels_cluster))
        else:
            results["ARI"] = float(adjusted_rand_score(labels_true, labels_cluster))

        if scib_nmi is not None:
            try:
                results["NMI"] = float(scib_nmi(adata, label_key, cluster_key))
            except Exception:
                results["NMI"] = float(
                    normalized_mutual_info_score(labels_true, labels_cluster)
                )
        else:
            results["NMI"] = float(
                normalized_mutual_info_score(labels_true, labels_cluster)
            )

        # AMI（sklearn）
        results["AMI"] = float(
            adjusted_mutual_info_score(labels_true, labels_cluster)
        )

        # Silhouette（按真值标签）
        if scib_silhouette is not None:
            try:
                results["SILHOUETTE_LABEL"] = float(
                    scib_silhouette(adata, label_key=label_key, embed=embed_key)
                )
            except Exception:
                results["SILHOUETTE_LABEL"] = float(
                    silhouette_score(X, labels_true, metric="euclidean")
                )
        else:
            results["SILHOUETTE_LABEL"] = float(
                silhouette_score(X, labels_true, metric="euclidean")
            )

    # 若提供批次，评估图连通性和批次ASW
    if batch_key is not None and batch_key in adata.obs.columns:
        # 图连通性
        if scib_graph_connectivity is not None:
            try:
                results["GRAPH_CONNECTIVITY"] = float(
                    scib_graph_connectivity(adata, label_key=batch_key)
                )
            except Exception:
                results["GRAPH_CONNECTIVITY"] = np.nan
                
        # 批次ASW
        batch_labels = adata.obs[batch_key].astype(str).values
        try:
            results["SILHOUETTE_BATCH"] = float(
                silhouette_score(X, batch_labels, metric="euclidean")
            )
        except Exception:
            results["SILHOUETTE_BATCH"] = np.nan

    return results


def to_long_records(
    metrics: Dict[str, float],
    *,
    method: str,
    resolution: float,
    seed: int,
    scenario: str,
    ndim: int,
    clustering_file: str,
    snn_file: str,
    dataset: str,
    dataset2: str,
) -> List[Dict[str, object]]:
    """将指标字典转换为长格式记录列表"""
    long_method = f"Leiden_on_{method}"
    records: List[Dict[str, object]] = []
    n_clusters = int(metrics.get("N_CLUSTERS", np.nan) if not np.isnan(metrics.get("N_CLUSTERS", np.nan)) else -1)
    for metric_name, value in metrics.items():
        if metric_name == "N_CLUSTERS":
            # 也记录为一条度量
            pass
        record = {
            "n_clusters": n_clusters,
            "metric": metric_name,
            "value": value,
            "method": method,
            "long_method": long_method,
            "feature_type": "embedding",
            "ndim": ndim,
            "resolution": resolution,
            "seed": seed,
            "scenario": scenario,
            "dataset": dataset,
            "dataset2": dataset2,
            "clustering_file": clustering_file,
            "snn_file": snn_file,
            # k_optimal 在后续 summary 聚合中更新；此处先与 n_clusters 对齐
            "k_optimal": n_clusters,
        }
        records.append(record)
    return records


def choose_best_resolution(
    df: pd.DataFrame,
    method: str,
    seed: int,
    prefer_metric: str,
) -> Tuple[float, int, float]:
    """为给定method/seed选择最优分辨率
    
    返回: (最优分辨率, 对应簇数, 对应指标分数)
    """
    sub = df[(df["method"] == method) & (df["seed"] == seed)]

    if prefer_metric not in sub["metric"].unique():
        # 回退策略
        fallback_order = ["SILHOUETTE_CLUSTER", "CH", "DB"]
        for fb in fallback_order:
            if fb in sub["metric"].unique():
                prefer_metric = fb
                break

    metric_df = sub[sub["metric"] == prefer_metric].copy()
    if metric_df.empty:
        return (np.nan, -1, np.nan)

    if prefer_metric in ["CH", "SILHOUETTE_CLUSTER", "ARI", "NMI", "AMI", "SILHOUETTE_LABEL", "GRAPH_CONNECTIVITY"]:
        best_row = metric_df.loc[metric_df["value"].idxmax()]
    else:  # DB 越小越好
        best_row = metric_df.loc[metric_df["value"].idxmin()]

    best_resolution = float(best_row["resolution"])
    best_k = int(best_row["n_clusters"]) if not pd.isna(best_row["n_clusters"]) else -1
    best_score = float(best_row["value"]) if not pd.isna(best_row["value"]) else np.nan
    return (best_resolution, best_k, best_score)


def main() -> None:
    """主函数:
    1. 解析参数
    2. 加载数据
    3. 对每个嵌入/种子/分辨率组合:
       - 构建KNN图
       - 执行聚类
       - 计算评估指标
    4. 汇总结果并输出
    """
    args = parse_args()

    adata = load_adata(args.adata)
    scenario = args.scenario or os.path.splitext(os.path.basename(args.adata))[0]

    # 解析分辨率
    resolutions: List[float]
    if args.res_grid:
        try:
            start_s, end_s, step_s = args.res_grid.split(":")
            start_v = float(start_s)
            end_v = float(end_s)
            step_v = float(step_s)
            vals = np.arange(start_v, end_v + 1e-12, step_v)
            # 规避浮点小误差
            resolutions = [float(f"{v:.6f}") for v in vals]
        except Exception as e:
            raise ValueError(f"--res-grid 解析失败，应为 start:end:step，例如 0.01:1.0:0.01；错误: {e}")
    else:
        resolutions = args.resolutions

    def fmt_res(v: float) -> str:
        s = f"{v:.6f}".rstrip("0").rstrip(".")
        return s if s != "" else "0"

    obsm_keys = list_obsm_keys(adata, args.obsm_keys)
    if len(obsm_keys) == 0:
        raise ValueError("未检测到任何 obsm 嵌入，请通过 --obsm-keys 指定")

    root_out = os.path.join(args.output_dir, scenario)
    ensure_dir(root_out)

    index_rows: List[Dict[str, object]] = []
    all_records: List[Dict[str, object]] = []

    for method in obsm_keys:
        if method not in adata.obsm:
            print(f"[警告] 跳过不存在的嵌入: {method}")
            continue
        ndim = int(adata.obsm[method].shape[1])

        for seed in args.seeds:
            # 每个 method/seed 先构图
            metric = "cosine" if method == "cellspace" else "euclidean"
            set_neighbors(adata, method, args.n_neighbors, seed, metric=metric)

            # 推断目标簇数（若未提供且可从标签推断）
            target_k_val = args.target_k
            if target_k_val is None and args.infer_target_k_from_label and args.label_key and args.label_key in adata.obs.columns:
                try:
                    target_k_val = int(adata.obs[args.label_key].nunique())
                except Exception:
                    target_k_val = None

            observed_max_k = -1
            for res in resolutions:
                # 聚类
                res_tag = fmt_res(res)
                clust_key = f"leiden_r{res_tag}_s{seed}_{method}"
                leiden_cluster(adata, res, seed, clust_key)

                # 计算指标
                metrics = compute_metrics(
                    adata,
                    embed_key=method,
                    cluster_key=clust_key,
                    label_key=args.label_key,
                    batch_key=args.batch_key,
                )

                # 输出路径
                eval_dir = os.path.join(
                    root_out, "evaluation", method, f"seed{seed}", f"r{res_tag}"
                )
                clust_dir = os.path.join(
                    root_out, "clustering", method, f"seed{seed}"
                )
                ensure_dir(eval_dir)
                ensure_dir(clust_dir)

                metrics_file = os.path.join(eval_dir, "metrics.tsv")
                clustering_file = os.path.join(clust_dir, f"r{res_tag}.tsv")

                # 保存聚类标签
                pd.DataFrame(
                    {
                        "cell": adata.obs_names,
                        "cluster": adata.obs[clust_key].astype(str).values,
                    }
                ).to_csv(clustering_file, sep="\t", index=False)

                # 记录长表
                records = to_long_records(
                    metrics,
                    method=method,
                    resolution=res,
                    seed=seed,
                    scenario=scenario,
                    ndim=ndim,
                    clustering_file=clustering_file,
                    snn_file="NA",
                    dataset=scenario.replace("_", " "),
                    dataset2=scenario,
                )
                df_metrics = pd.DataFrame.from_records(records)
                df_metrics.to_csv(metrics_file, sep="\t", index=False)

                all_records.extend(records)

                # 写入 index 行（模拟示例中的文件索引）
                index_rows.append(
                    {
                        "file": os.path.relpath(metrics_file, start=root_out),
                        "scenario": scenario,
                        "workflow": "evaluation",
                        "method": method,
                        "feature_type": "embedding",
                        "tile_size": 0,
                        "distance": "default",
                        "ndim": ndim,
                        "seed": seed,
                        "filename": os.path.basename(metrics_file),
                        "resolution": res,
                        "rds_file": "NA",
                        "clustering_file": os.path.relpath(
                            clustering_file, start=root_out
                        ),
                        "long_method": f"Leiden_on_{method}",
                        "snn_file": "NA",
                        "n_clusters": int(metrics.get("N_CLUSTERS", np.nan))
                        if not np.isnan(metrics.get("N_CLUSTERS", np.nan))
                        else -1,
                        "dataset": scenario.replace("_", " "),
                        "dataset2": scenario,
                        "k_optimal": int(metrics.get("N_CLUSTERS", np.nan))
                        if not np.isnan(metrics.get("N_CLUSTERS", np.nan))
                        else -1,
                    }
                )

                # 早停：根据簇数上界覆盖判断
                try:
                    current_k = int(metrics.get("N_CLUSTERS", np.nan)) if not np.isnan(metrics.get("N_CLUSTERS", np.nan)) else -1
                except Exception:
                    current_k = -1
                observed_max_k = max(observed_max_k, current_k)
                if args.early_stop and target_k_val is not None and observed_max_k >= (target_k_val + args.k_margin):
                    # 已达到期望覆盖，提前结束该 method/seed 的扫描
                    break

    # 汇总输出（单组合模式关闭）
    if not args.single_combo:
        index_df = pd.DataFrame.from_records(index_rows)
        index_path = os.path.join(root_out, "index.tsv")
        if not index_df.empty:
            index_df.to_csv(index_path, sep="\t", index=False)

        summary_df = pd.DataFrame.from_records(all_records)
        summary_path = os.path.join(root_out, "summary.tsv")
        if not summary_df.empty:
            summary_df.to_csv(summary_path, sep="\t", index=False)

        # 选择每个 method/seed 的最优分辨率
        best_rows: List[Dict[str, object]] = []
        if not summary_df.empty:
            for method in summary_df["method"].unique():
                for seed in summary_df["seed"].unique():
                    best_res, best_k, best_score = choose_best_resolution(
                        summary_df, method, seed, args.select_metric
                    )
                    if not np.isnan(best_res):
                        best_rows.append(
                            {
                                "method": method,
                                "seed": int(seed),
                                "best_resolution": best_res,
                                "k_optimal": best_k,
                                "select_metric": args.select_metric,
                                "best_score": best_score,
                            }
                        )
        best_df = pd.DataFrame.from_records(best_rows)
        best_path = os.path.join(root_out, "best_per_method.tsv")
        if not best_df.empty:
            best_df.to_csv(best_path, sep="\t", index=False)

    print(f"完成。输出目录: {root_out}")


if __name__ == "__main__":
    main()
