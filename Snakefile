import os

# 配置
configfile: "config.yaml"

SCENARIO = config.get("scenario")
ADATA = config.get("adata")
LABEL_KEY = config.get("label_key", None)
BATCH_KEY = config.get("batch_key", None)
OBSM_KEYS = config.get("obsm_keys", [])
SEEDS = config.get("seeds", [0, 2, 5])
N_NEIGHBORS = int(config.get("n_neighbors", 15))
SELECT_METRIC = config.get("select_metric", "ARI")
OUTPUT_DIR = config.get("output_dir", "./py_benchmark_results")

# 分辨率网格（扩展到2.0），可用两种方式：显式列表 或 start:end:step
RESOLUTIONS = config.get("resolutions", [])
RES_GRID = config.get("res_grid", None)

def gen_resolutions():
    if RES_GRID:
        start_s, end_s, step_s = RES_GRID.split(":")
        start_v = float(start_s)
        end_v = float(end_s)
        step_v = float(step_s)
        vals = []
        v = start_v
        while v <= end_v + 1e-12:
            vals.append(float(f"{v:.6f}"))
            v += step_v
        return vals
    if RESOLUTIONS:
        return [float(x) for x in RESOLUTIONS]
    # 默认：0.02 步，直到 2.00
    return [float(f"{i/50:.2f}") for i in range(1, 101)]

RES_LIST = gen_resolutions()

# 提取真实标签规则
rule extract_true_labels:
    input:
        adata=ADATA
    output:
        truth_tsv=os.path.join(OUTPUT_DIR, SCENARIO, "ground_truth.tsv")
    params:
        label_key=LABEL_KEY
    conda:
        "env.yaml"
    shell:
        "/opt/anaconda3/envs/scib-benchmark/bin/python extract_truth_labels.py {input.adata} {params.label_key} {output.truth_tsv}"

# 单组合运行，生成 metrics 与 clustering
rule run_single_combo:
    conda:
        "env.yaml"
    output:
        metrics=os.path.join(
            OUTPUT_DIR, SCENARIO, "evaluation", "{method}", "seed{seed}", "r{res_tag}", "metrics.tsv"
        ),
        clustering=os.path.join(
            OUTPUT_DIR, SCENARIO, "clustering", "{method}", "seed{seed}", "r{res_tag}.tsv"
        ),
    params:
        res=lambda wildcards: float(wildcards.res_tag)
    wildcard_constraints:
        res_tag=r"[0-9]+\.?[0-9]*",
    shell:
        (
            f"/opt/anaconda3/envs/scib-benchmark/bin/python run_benchmark.py "
            f"--adata {ADATA} "
            + (f"--label-key {LABEL_KEY} " if LABEL_KEY else "")
            + (f"--batch-key {BATCH_KEY} " if BATCH_KEY else "")
            + "--obsm-keys {wildcards.method} --resolutions {params.res} "
            + "--seeds {wildcards.seed} "
            + f"--n-neighbors {N_NEIGHBORS} --scenario {SCENARIO} --output-dir {OUTPUT_DIR} "
            + ("--infer-target-k-from-label " if LABEL_KEY else "")
            + "--early-stop "
            + "--single-combo"
        )

# 生成单个混淆矩阵图规则
rule plot_single_confusion_matrix:
    input:
        clustering=os.path.join(
            OUTPUT_DIR, SCENARIO, "clustering", "{method}", "seed{seed}", "r{res_tag}.tsv"
        ),
        truth=rules.extract_true_labels.output.truth_tsv
    output:
        plot=os.path.join(
            OUTPUT_DIR, SCENARIO, "evaluation", "{method}", "seed{seed}", "r{res_tag}", "confusion_matrix.png"
        )
    conda:
        "env.yaml"
    params:
        r_script="plot_confusion_matrix.R" # R脚本路径
    shell:
        "if [ -s {input.truth} ]; then /opt/anaconda3/envs/scib-benchmark/bin/Rscript {params.r_script} --clusters {input.clustering} --truth {input.truth} --output {output.plot}; else touch {output.plot}; fi"

# 聚合所有 metrics 生成 summary/index/best
rule aggregate:
    conda:
        "env.yaml"
    input:
        expand(
            os.path.join(OUTPUT_DIR, SCENARIO, "evaluation", "{method}", "seed{seed}", "r{res_tag}", "metrics.tsv"),
            method=OBSM_KEYS,
            seed=SEEDS,
            res_tag=[str(r).rstrip("0").rstrip(".") if isinstance(r, float) else str(r) for r in RES_LIST],
        ),
    output:
        summary=os.path.join(OUTPUT_DIR, SCENARIO, "summary.tsv"),
        index=os.path.join(OUTPUT_DIR, SCENARIO, "index.tsv"),
        best=os.path.join(OUTPUT_DIR, SCENARIO, "best_per_method.tsv"),
    params:
        methods=OBSM_KEYS,
        seeds=SEEDS,
        select_metric=SELECT_METRIC,
    shell:
        (
            "python aggregate_results.py "
            "--root {root} --scenario {sc} --select-metric {sm}"
        ).format(
            root=os.path.join(OUTPUT_DIR, SCENARIO),
            sc=SCENARIO,
            sm=SELECT_METRIC,
        )

# 动态目标：所有组合
rule cluster_combo:
    input:
        rules.run_single_combo.output

# 便捷目标：展开所有组合
rule expand_all_combos:
    input:
        expand(
            os.path.join(OUTPUT_DIR, SCENARIO, "evaluation", "{method}", "seed{seed}", "r{res_tag}", "metrics.tsv"),
            method=OBSM_KEYS,
            seed=SEEDS,
            res_tag=[str(r).rstrip("0").rstrip(".") if isinstance(r, float) else str(r) for r in RES_LIST],
        )

# 主目标规则
rule all:
    input:
        os.path.join(OUTPUT_DIR, SCENARIO, "summary.tsv"),
        os.path.join(OUTPUT_DIR, SCENARIO, "index.tsv"),
        os.path.join(OUTPUT_DIR, SCENARIO, "best_per_method.tsv"),
        # 添加混淆矩阵图的展开目标
        expand(
            rules.plot_single_confusion_matrix.output.plot,
            method=OBSM_KEYS,
            seed=SEEDS,
            res_tag=[str(r).rstrip("0").rstrip(".") if isinstance(r, float) else str(r) for r in RES_LIST],
        ),
