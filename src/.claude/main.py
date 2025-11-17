import pandas as pd

import pyarrow.parquet as pq


import pandas as pd

df = pd.read_parquet(r"E:\UIC_PHD\Matriviz_project\RNA_h5_test\pnsdf_rna_centroid_v1.parquet")
print("Cortex columns:", df.columns.tolist())
print("Cortex shape:", df.shape)
print("First few rows:")
print(df.head())


df = pd.read_parquet(r"E:\UIC_PHD\Matriviz_project\MatriViz-Data\cortex_centroid_v1.parquet")
print("PNSDF columns:", df.columns.tolist())
print("PNSDF shape:", df.shape)
print("First few rows:")
print(df.head())

# # Replace with your file path
file = pq.ParquetFile(r"E:\UIC_PHD\Matriviz_project\RNA_h5_test\asdfasdf_rna_centroid_v1.parquet")
file1 = pq.ParquetFile(r"E:\UIC_PHD\Matriviz_project\MatriViz-Data\cortex_centroid_v1.parquet")
#
#
# # Show metadata summary
print("pnsdf",file.metadata)
print("cortex",file1.metadata)
#
# # Or access Parquet format version directly
# print("Parquet format pnsdf:", file.metadata.format_version)
# print("Parquet format cortex:", file1.metadata.format_version)
#
#
#
# df = pd.read_parquet(r"E:\UIC_PHD\Matriviz_project\RNA_h5_test\pnsdf_rna_centroid_v1.parquet")
# df1 = pd.read_parquet(r"E:\UIC_PHD\Matriviz_project\MatriViz-Data\cortex_centroid_v1.parquet")