import scanpy as sc
import os
import glob
import numpy as np  # 用于设置全局种子
from tqdm import tqdm

# --- 配置路径 ---
ref_path = r"E:\UIC_PHD\Python_projects\parquet_read\src\map_h5ad\refs\bladder_ref.h5ad"
input_dir = r"F:\12_30_2025_test_azimuth_annotation"
output_dir = r"E:\UIC_PHD\PythonProject\12-4-2025-Matriviz\raw\Mapped_bladder_h5ad"

# 设置全局随机种子（双重保险）
SEED = 42
np.random.seed(SEED)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# --- 1. 准备 Reference 数据 ---
print("Loading reference data...")
ref_raw = sc.read_h5ad(ref_path)
ref_raw.var_names = ref_raw.var_names.str.split('.').str[0]

# --- 2. 获取文件列表 ---
query_files = glob.glob(os.path.join(input_dir, "*.h5ad"))
pbar = tqdm(query_files, desc="Processing Samples", unit="sample")

for query_path in pbar:
    file_name = os.path.basename(query_path)
    pbar.set_description(f"Processing: {file_name}")
    save_path = os.path.join(output_dir, file_name)

    try:
        adata = sc.read_h5ad(query_path)
        adata.var_names = adata.var_names.str.split('.').str[0]

        # 1. 找出共同基因
        common_genes = ref_raw.var_names.intersection(adata.var_names)

        # 2. 提取子集
        ref_sub = ref_raw[:, common_genes].copy()
        adata_sub = adata[:, common_genes].copy()

        # 3. 对针对该样本的 ref_sub 进行快速预处理
        # 在 PCA 和 UMAP 中显式设置 random_state
        sc.pp.pca(ref_sub, random_state=SEED)
        sc.pp.neighbors(ref_sub, n_neighbors=14, random_state=SEED)
        sc.tl.umap(ref_sub, random_state=SEED)

        # 4. 映射
        # ingest 内部会继承 ref_sub 的布局信息
        sc.tl.ingest(adata_sub, ref_sub, obs='cell_type')

        # 5. 整理元数据
        adata_sub.obs['azimuth_label'] = adata_sub.obs['cell_type']
        adata_sub.obsm['X_umap_proj'] = adata_sub.obsm['X_umap']

        # 6. 保存
        adata_sub.write(save_path)

    except Exception as e:
        tqdm.write(f"Error processing {file_name}: {e}")

print("\nAll tasks completed! Random state fixed at 42.")