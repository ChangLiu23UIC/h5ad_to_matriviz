import scanpy as sc
import os
import glob
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from tqdm import tqdm

# --- 配置路径 ---
ref_path = r"E:\UIC_PHD\Python_projects\parquet_read\src\map_h5ad\refs\lymph_node_ref.h5ad"
input_dir = r"F:\12_30_2025_test_azimuth_annotation"
output_dir = r"E:\UIC_PHD\PythonProject\12-4-2025-Matriviz\raw\Mapped_lymph_node_h5ad"

# 设置全局随机种子
SEED = 42
np.random.seed(SEED)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# --- 1. 准备 Reference 数据 (在循环外执行一次) ---
print("Loading and preprocessing reference data...")
ref_raw = sc.read_h5ad(ref_path)
# 清理基因名版本号
ref_raw.var_names = ref_raw.var_names.str.split('.').str[0]

# 预处理 Reference：确保其具有 PCA 和 Neighbors 布局
# ingest 需要参考集已经完成了这些步骤
sc.pp.pca(ref_raw, random_state=SEED)
sc.pp.neighbors(ref_raw, n_neighbors=14, random_state=SEED)
sc.tl.umap(ref_raw, random_state=SEED)

# --- 2. 获取文件列表 ---
query_files = glob.glob(os.path.join(input_dir, "*.h5ad"))
pbar = tqdm(query_files, desc="Processing Samples", unit="sample")

for query_path in pbar:
    file_name = os.path.basename(query_path)
    pbar.set_description(f"Processing: {file_name}")
    save_path = os.path.join(output_dir, file_name)

    try:
        # A. 读取 Query 数据
        adata = sc.read_h5ad(query_path)
        adata.var_names = adata.var_names.str.split('.').str[0]

        # B. 基因对齐与补零 (关键步骤)
        # 创建一个全零矩阵，形状为 (Query细胞数, Reference基因数)
        n_obs = adata.n_obs
        n_vars_ref = ref_raw.n_vars
        new_X = csr_matrix((n_obs, n_vars_ref))

        # 创建用于 Mapping 的临时对象，其 var 必须与 ref_raw 完全一致
        adata_map = sc.AnnData(X=new_X, obs=adata.obs.copy(), var=ref_raw.var.copy())

        # 找出 Query 中存在的基因
        common_genes = adata.var_names.intersection(ref_raw.var_names)

        # 将 unscaled 层中存在的基因数据填充到全零矩阵中
        if "unscaled" in adata.layers:
            # 提取 unscaled 层中共同基因的数据
            source_data = adata[:, common_genes].layers["unscaled"]
            # 映射到 adata_map 中的对应位置
            # 使用 pandas 指数匹配来确保位置准确
            var_index_map = {gene: i for i, gene in enumerate(ref_raw.var_names)}
            target_indices = [var_index_map[gene] for gene in common_genes]

            # 将稀疏矩阵的数据填入
            # 注意：AnnData 的 .X 赋值需要处理
            temp_X = adata_map.X.tolil()
            temp_X[:, target_indices] = source_data
            adata_map.X = temp_X.tocsr()

            # C. 归一化对齐
            sc.pp.normalize_total(adata_map, target_sum=1e4)
            sc.pp.log1p(adata_map)
        else:
            tqdm.write(f"Warning: 'unscaled' layer not found in {file_name}. Skipped padding.")
            continue

        # D. 执行映射 (使用循环外已经算好的 ref_raw)
        # 此时 adata_map.var_names 和 ref_raw.var_names 完全一致
        sc.tl.ingest(adata_map, ref_raw, obs='cell_type')

        # E. 将结果同步回原始 adata
        # 我们保留原始的 adata（包含原本的 .X 和所有层），仅更新预测的元数据
        adata.obs['azimuth_label'] = adata_map.obs['cell_type']
        adata.obsm['X_umap_proj'] = adata_map.obsm['X_umap']

        # F. 保存 (保留了原始 .X 的内容)
        adata.write(save_path)

    except Exception as e:
        tqdm.write(f"Error processing {file_name}: {e}")

print("\nAll tasks completed! Random state fixed at 42.")