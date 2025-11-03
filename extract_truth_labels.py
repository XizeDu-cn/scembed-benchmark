#!/usr/bin/env python3

import sys
import scanpy as sc
import pandas as pd

def main():
    if len(sys.argv) != 4:
        print("用法: python extract_truth_labels.py <adata_file> <label_key> <output_file>")
        sys.exit(1)
    
    adata_file = sys.argv[1]
    label_key = sys.argv[2]
    output_file = sys.argv[3]
    
    # 读取AnnData文件
    adata = sc.read_h5ad(adata_file)
    
    if label_key in adata.obs.columns:
        # 创建真实标签数据框
        df = pd.DataFrame({
            'cell': adata.obs.index,
            'truth_label': adata.obs[label_key]
        })
        df.to_csv(output_file, sep='\t', index=False)
        print(f"真实标签已保存到: {output_file}")
    else:
        # 如果没有标签键，创建空文件
        with open(output_file, 'w') as f:
            f.write("cell\ttruth_label\n")
        print(f"未找到标签键 '{label_key}'，创建空文件: {output_file}")

if __name__ == "__main__":
    main()


