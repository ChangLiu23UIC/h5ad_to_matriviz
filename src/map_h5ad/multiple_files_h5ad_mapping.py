import scanpy as sc
import os
import glob
from tqdm import tqdm  # 导入进度条库

# --- 配置路径 ---
ref_path = r"E:\UIC_PHD\Python_projects\parquet_read\src\map_h5ad\refs\bladder_ref.h5ad"
input_dir = r"G:\Secondary_analysis\Organized_by_Tissue\Bladder"
output_dir = r"F:\12_30_2025_test_azimuth_annotation\bladder_mapped"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# --- 1. 准备 Reference 数据 ---
print("Loading reference data...")
ref_data = sc.read_h5ad(ref_path)
ref_data.var_names = ref_data.var_names.str.split('.').str[0]

# --- 2. 获取文件列表 ---
query_files = glob.glob(os.path.join(input_dir, "*.h5ad"))

# --- 3. 使用 tqdm 遍历 ---
# desc 是进度条前的描述文字，unit 是单位
pbar = tqdm(query_files, desc="Processing Samples", unit="sample")

for query_path in pbar:
    file_name = os.path.basename(query_path)

    # 实时更新进度条左侧的文字，显示当前处理的文件
    pbar.set_description(f"Processing: {file_name}")

    save_path = os.path.join(output_dir, file_name)

    try:
        adata = sc.read_h5ad(query_path)
        adata.var_names = adata.var_names.str.split('.').str[0]

        common_genes = adata.var_names.intersection(ref_data.var_names)
        adata_sub = adata[:, common_genes].copy()
        ref_sub = ref_data[:, common_genes].copy()
        adata_sub = adata_sub[:, ref_sub.var_names].copy()

        # 计算
        sc.pp.pca(ref_sub)
        sc.pp.neighbors(ref_sub, n_neighbors=14, use_rep='X_pca')
        sc.tl.umap(ref_sub, random_state=42)

        sc.tl.ingest(adata_sub, ref_sub, obs='cell_type')

        adata_sub.obs['azimuth_label'] = adata_sub.obs['cell_type']
        adata_sub.obsm['X_umap_proj'] = adata_sub.obsm['X_umap']



    except Exception as e:
        # 使用 tqdm.write 可以在不破坏进度条的情况下打印错误信息
        tqdm.write(f"Error processing {file_name}: {e}")

print("\nAll tasks completed!")