from cgitb import small

import numpy as np
import scanpy as sc
import pandas as pd

bladder_path = r"E:\UIC_PHD\Python_projects\parquet_read\src\map_h5ad\refs\bladder_ref.h5ad"
knee_path = r"E:\UIC_PHD\Python_projects\parquet_read\src\map_h5ad\refs\knee_ref.h5ad"
large_path = r"E:\UIC_PHD\Python_projects\parquet_read\src\map_h5ad\refs\large_intestine_ref.h5ad"
lymph_path = r"E:\UIC_PHD\Python_projects\parquet_read\src\map_h5ad\refs\lymph_node_ref.h5ad"
small_path = r"E:\UIC_PHD\Python_projects\parquet_read\src\map_h5ad\refs\small_intestine_ref.h5ad"
spleen_path = r"E:\UIC_PHD\Python_projects\parquet_read\src\map_h5ad\refs\spleen_ref.h5ad"

df = sc.read_h5ad(spleen_path)

df_cen = df.obs
lst = list(set(df_cen["cell_type"].to_list()))

df_final = pd.DataFrame(lst, columns=["Type"])
df_final.to_parquet("centroids/spleen_centroid.parquet")