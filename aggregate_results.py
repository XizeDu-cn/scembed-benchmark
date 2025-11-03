#!/usr/bin/env python3
import argparse
import os
import glob
import numpy as np
import pandas as pd
from typing import List, Tuple


def parse_args():
    p = argparse.ArgumentParser("聚合 run_benchmark 单组合输出，生成 summary/index/best")
    p.add_argument("--root", required=True, help="场景根目录，如 py_benchmark/results/Cell_line_mixing")
    p.add_argument("--scenario", required=True)
    p.add_argument("--select-metric", default="ARI")
    return p.parse_args()


def load_all_metrics(root: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    metrics_files = glob.glob(os.path.join(root, "evaluation", "*", "seed*", "r*", "metrics.tsv"))
    records: List[pd.DataFrame] = []
    index_rows: List[dict] = []
    for mf in metrics_files:
        df = pd.read_csv(mf, sep="\t")
        records.append(df)
        parts = mf.split(os.sep)
        # 路径结构: .../evaluation/{method}/seed{seed}/r{res_tag}/metrics.tsv
        method = parts[-4]
        seed = int(parts[-3].replace("seed", ""))
        res_tag = parts[-2].replace("r", "")
        ndim = int(df["ndim"].iloc[0]) if "ndim" in df.columns else -1
        n_clusters = int(df["n_clusters"].iloc[0]) if "n_clusters" in df.columns else -1
        index_rows.append(
            {
                "file": os.path.relpath(mf, start=root),
                "scenario": df["scenario"].iloc[0] if "scenario" in df.columns else "",
                "workflow": "evaluation",
                "method": method,
                "feature_type": "embedding",
                "tile_size": 0,
                "distance": "default",
                "ndim": ndim,
                "seed": seed,
                "filename": os.path.basename(mf),
                "resolution": res_tag,
                "rds_file": "NA",
                "clustering_file": os.path.relpath(
                    os.path.join(root, "clustering", method, f"seed{seed}", f"r{res_tag}.tsv"), start=root
                ),
                "long_method": f"Leiden_on_{method}",
                "snn_file": "NA",
                "n_clusters": n_clusters,
                "dataset": df["dataset"].iloc[0] if "dataset" in df.columns else "",
                "dataset2": df["dataset2"].iloc[0] if "dataset2" in df.columns else "",
                "k_optimal": n_clusters,
            }
        )
    summary = pd.concat(records, ignore_index=True) if records else pd.DataFrame()
    index_df = pd.DataFrame(index_rows)
    return summary, index_df


def choose_best_resolution(df: pd.DataFrame, prefer_metric: str) -> pd.DataFrame:
    best_rows = []
    if df.empty:
        return pd.DataFrame()
    for method in df["method"].unique():
        for seed in df["seed"].unique():
            sub = df[(df["method"] == method) & (df["seed"] == seed)]
            metric = prefer_metric
            if metric not in sub["metric"].unique():
                for fb in ["SILHOUETTE_CLUSTER", "CH", "DB"]:
                    if fb in sub["metric"].unique():
                        metric = fb
                        break
            mdf = sub[sub["metric"] == metric]
            if mdf.empty:
                continue
            if metric in ["CH", "SILHOUETTE_CLUSTER", "ARI", "NMI", "AMI", "SILHOUETTE_LABEL", "GRAPH_CONNECTIVITY"]:
                row = mdf.loc[mdf["value"].idxmax()]
            else:
                row = mdf.loc[mdf["value"].idxmin()]
            best_rows.append(
                {
                    "method": method,
                    "seed": int(seed),
                    "best_resolution": row["resolution"],
                    "k_optimal": int(row["n_clusters"]) if not pd.isna(row["n_clusters"]) else -1,
                    "select_metric": metric,
                    "best_score": float(row["value"]) if not pd.isna(row["value"]) else np.nan,
                }
            )
    return pd.DataFrame(best_rows)


def main():
    args = parse_args()
    summary, index_df = load_all_metrics(args.root)
    if not summary.empty:
        summary.to_csv(os.path.join(args.root, "summary.tsv"), sep="\t", index=False)
    if not index_df.empty:
        index_df.to_csv(os.path.join(args.root, "index.tsv"), sep="\t", index=False)
    best_df = choose_best_resolution(summary, args.select_metric)
    if not best_df.empty:
        best_df.to_csv(os.path.join(args.root, "best_per_method.tsv"), sep="\t", index=False)
    print("聚合完成", args.root)


if __name__ == "__main__":
    main()
