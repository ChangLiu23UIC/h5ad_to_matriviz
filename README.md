# H5AD to MatriViz Converter

This tool converts h5ad files (AnnData format) to MatriViz-compatible parquet format for visualization in the MatriViz application.

## Features

- Convert h5ad files to parquet format compatible with MatriViz
- Extract UMAP coordinates and gene expression data
- Generate category files and overview JSON
- GUI and command-line interfaces
- Support for h5seurat conversion (requires R)

## Installation

1. Install Python dependencies:
```bash
pip install -r requirements.txt
```

2. For h5seurat support, install R and SeuratDisk:
```R
install.packages("SeuratDisk")
```

## Usage

### GUI Interface

Run the GUI application:
```bash
python src/file_to_parquet.py
```

1. Select your input .h5ad file
2. Enter tissue name (used for file naming)
3. Select output directory
4. (Optional) Select category JSON template
5. Click "Convert"

### Command Line Interface

Use the test script for batch processing:
```bash
python src/test_conversion.py input.h5ad "Tissue Name" ./output
```

Optional: Specify category JSON template:
```bash
python src/test_conversion.py input.h5ad "Lung" ./output template.json
```

## Output Files

After conversion, you'll get these files:

- `{tissue}_v1.parquet` - Main expression data with UMAP coordinates
- `{tissue}_centroid_v1.parquet` - Cell type centroids (if available)
- `{tissue}_v1_category.json` - Gene category definitions
- `{tissue}.json` - Overview configuration for MatriViz

## MatriViz Compatibility

The converter produces files with this structure:

### Parquet File Schema
- `index` (string): Cell identifier
- `umap_1` (numeric): X coordinate
- `umap_2` (numeric): Y coordinate
- `gene1`, `gene2`, ... (numeric): Gene expression values

### Overview JSON
```json
{
  "fileType": "matriviz",
  "version": "0.0.1",
  "category_name": "tissue_name",
  "category_description": "Tissue Name",
  "parquet_file": "tissue_name_v1.parquet",
  "category_file": "tissue_name_v1_category.json",
  "centroid_file": "tissue_name_centroid_v1.parquet"
}
```

## Requirements

- Python 3.8+
- h5ad file with UMAP coordinates in `adata.obsm['X_umap']`
- Gene expression matrix in `adata.X`

## Troubleshooting

### Common Issues

1. **"No UMAP coordinates found"**
   - Ensure your h5ad file has UMAP coordinates in `adata.obsm['X_umap']`

2. **Memory errors with large files**
   - The converter loads the full expression matrix into memory
   - Consider subsetting large datasets before conversion

3. **Column name conflicts**
   - The converter automatically removes columns containing "celltype"
   - Check your column names if genes are missing

### Performance Tips

- Use zstd compression for smaller file sizes
- Consider filtering lowly expressed genes before conversion
- For very large datasets, use the command-line interface for batch processing

## License

This tool is provided for use with the MatriViz application.