import sys, os, json
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import anndata as ad
from pathlib import Path
from PyQt6.QtWidgets import (
    QApplication, QWidget, QLabel, QPushButton, QLineEdit, QFileDialog,
    QVBoxLayout, QHBoxLayout, QMessageBox, QProgressBar
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal


# =========================
# Worker Thread (Background Conversion)
# =========================
class ConverterThread(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(bool, str)

    def __init__(self, input_file, tissue, output_dir):
        super().__init__()
        self.input_file = str(input_file)
        self.tissue = str(tissue).strip()
        self.output_dir = str(output_dir)

    def run(self):
        try:
            self.progress.emit("Loading data...")
            os.makedirs(self.output_dir, exist_ok=True)

            # --- 支持 h5ad 和 h5/h5seurat ---
            if self.input_file.endswith(".h5ad"):
                adata = ad.read_h5ad(self.input_file)
            elif self.input_file.endswith(".h5") or self.input_file.endswith(".h5seurat"):
                self.progress.emit("Converting .h5seurat → .h5ad ...")
                temp_h5ad = os.path.join(self.output_dir, "temp_converted.h5ad")

                import subprocess
                r_script = f"""
                suppressMessages(library(SeuratDisk));
                Convert("{self.input_file}", dest = "h5ad", overwrite = TRUE);
                """
                result = subprocess.run(["Rscript", "-e", r_script], capture_output=True, text=True)
                if result.returncode != 0:
                    raise RuntimeError(f"R conversion failed:\\n{result.stderr}")

                adata = ad.read_h5ad(temp_h5ad)
            else:
                raise ValueError("Supported input: .h5ad or .h5/.h5seurat")

            print("=== AnnData loaded ===")
            print(f"obs: {adata.obs_names.shape}, vars: {adata.var_names.shape}")

            # --- Expression Matrix ---
            self.progress.emit("Extracting expression matrix...")
            expr = adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X
            expr_df = pd.DataFrame(expr, index=adata.obs_names, columns=adata.var_names)
            expr_df.columns = expr_df.columns.astype(str)

            # --- UMAP Coordinates ---
            if "X_umap" not in adata.obsm:
                raise ValueError("No UMAP coordinates found in object.")
            umap_df = pd.DataFrame(
                adata.obsm["X_umap"], columns=["UMAP_1", "UMAP_2"], index=adata.obs_names
            )
            umap_df["index"] = umap_df.index

            # --- Metadata ---
            meta_df = adata.obs.copy()
            meta_df["index"] = meta_df.index
            celltype_cols = [c for c in meta_df.columns if "celltype" in c.lower()]
            label_col = celltype_cols[-1] if celltype_cols else None

            # --- Merge Expression + UMAP ---
            self.progress.emit("Merging data...")
            merged = umap_df.join(expr_df, how="inner")
            numeric_part = merged.iloc[:, 2:].apply(pd.to_numeric, errors="coerce")
            merged["All_Genes"] = numeric_part.mean(axis=1, skipna=True)
            merged = merged.loc[:, ~merged.columns.str.contains("celltype", case=False)]
            merged.columns = merged.columns.astype(str)

            # =============================
            # 🔹 File Naming Convention
            # =============================
            base_name = str(self.tissue).replace(" ", "_")  # preserve capitalization
            expr_filename = f"{base_name}_v1_ss.parquet"
            centroid_filename = f"{base_name}_centroid_v1.parquet"
            category_filename = f"{base_name}_v1_category.json"
            overview_filename = f"{base_name}.json"

            expr_path = Path(self.output_dir) / expr_filename
            centroid_path = Path(self.output_dir) / centroid_filename
            category_path = Path(self.output_dir) / category_filename
            overview_path = Path(self.output_dir) / overview_filename

            # --- Write Expression Parquet ---
            self.progress.emit("Writing expression parquet...")
            pq.write_table(pa.Table.from_pandas(merged.reset_index(drop=True)), expr_path, compression="zstd")

            # --- Compute Centroids ---
            self.progress.emit("Computing centroids...")
            if label_col:
                meta_sub = meta_df[[label_col, "index"]].merge(umap_df, on="index", how="left")
                centroids = (
                    meta_sub.groupby(label_col)[["UMAP_1", "UMAP_2"]]
                    .median()
                    .reset_index()
                    .rename(columns={label_col: "Type", "UMAP_1": "cen_x", "UMAP_2": "cen_y"})
                )
            else:
                centroids = pd.DataFrame(columns=["Type", "cen_x", "cen_y"])

            pq.write_table(pa.Table.from_pandas(centroids), centroid_path, compression="zstd")

            # =============================
            # 🔹 Matrisome Category JSON (matches cortex_v1_category.json)
            # =============================
            self.progress.emit("Generating matrisome category JSON...")
            matrisome_path = Path(__file__).parent / "matrisome_list.json"
            if not matrisome_path.exists():
                raise FileNotFoundError(f"matrisome_list.json not found at {matrisome_path}")

            with open(matrisome_path, "r", encoding="utf-8") as f:
                matrisome_data = json.load(f)

            categories = {
                "ECM Glycoproteins": [],
                "Secreted Factors": [],
                "Collagens": [],
                "ECM Regulators": [],
                "Proteoglycans": [],
                "ECM-affiliated Proteins": []
            }

            # Flatten gene mapping from matrisome list
            def collect_genes(node, parent_key=None):
                if isinstance(node, dict):
                    for k, v in node.items():
                        collect_genes(v, k)
                elif isinstance(node, list) and parent_key in categories:
                    categories[parent_key].extend(node)

            collect_genes(matrisome_data)

            # Filter by expressed genes
            all_genes = set(expr_df.columns.str.upper())
            for key in categories:
                categories[key] = sorted([g for g in categories[key] if g.upper() in all_genes])

            # Write the category JSON
            with open(category_path, "w", encoding="utf-8") as f:
                json.dump(categories, f, indent=4)

            # =============================
            # 🔹 Overview JSON (matches cortex.json)
            # =============================
            self.progress.emit("Writing overview JSON...")
            overview_json = {
                "fileType": "matriviz",
                "version": "0.0.1",
                "category_name": base_name.lower(),
                "category_description": str(self.tissue),
                "parquet_file": expr_filename,
                "category_file": category_filename,
                "centroid_file": centroid_filename
            }

            with open(overview_path, "w", encoding="utf-8") as f:
                json.dump(overview_json, f, indent=4)

            self.finished.emit(True, f"✅ Conversion completed!\n\nOutput: {self.output_dir}")

        except Exception as e:
            import traceback
            tb = traceback.format_exc()
            self.finished.emit(False, f"❌ Error: {str(e)}\n\nTraceback:\n{tb}")


# =========================
# Main GUI Window
# =========================
class AzimuthApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Azimuth → MatriViz Converter")
        self.setFixedSize(480, 300)

        layout = QVBoxLayout()

        self.input_label = QLabel("Input .h5ad file:")
        self.input_line = QLineEdit()
        self.input_btn = QPushButton("Browse")
        self.input_btn.clicked.connect(self.select_input)
        hl1 = QHBoxLayout()
        hl1.addWidget(self.input_line)
        hl1.addWidget(self.input_btn)

        self.tissue_label = QLabel("Tissue name:")
        self.tissue_line = QLineEdit()

        self.output_label = QLabel("Output directory:")
        self.output_line = QLineEdit()
        self.output_btn = QPushButton("Browse")
        self.output_btn.clicked.connect(self.select_output)
        hl2 = QHBoxLayout()
        hl2.addWidget(self.output_line)
        hl2.addWidget(self.output_btn)

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(False)

        self.convert_btn = QPushButton("Convert")
        self.convert_btn.clicked.connect(self.start_conversion)

        layout.addWidget(self.input_label)
        layout.addLayout(hl1)
        layout.addWidget(self.tissue_label)
        layout.addWidget(self.tissue_line)
        layout.addWidget(self.output_label)
        layout.addLayout(hl2)
        layout.addWidget(self.convert_btn)
        layout.addWidget(self.progress_bar)
        self.setLayout(layout)

    def select_input(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select .h5ad File", "", "H5AD files (*.h5ad)")
        if file:
            self.input_line.setText(file)

    def select_output(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if folder:
            self.output_line.setText(folder)

    def start_conversion(self):
        input_file = self.input_line.text().strip()
        tissue = self.tissue_line.text().strip()
        output_dir = self.output_line.text().strip()

        if not all([input_file, tissue, output_dir]):
            QMessageBox.warning(self, "Missing Info", "Please fill all fields.")
            return

        self.convert_btn.setEnabled(False)
        self.progress_bar.setVisible(True)

        self.worker = ConverterThread(input_file, tissue, output_dir)
        self.worker.progress.connect(self.show_progress)
        self.worker.finished.connect(self.finish_conversion)
        self.worker.start()

    def show_progress(self, msg):
        self.progress_bar.setFormat(msg)

    def finish_conversion(self, success, msg):
        self.convert_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        QMessageBox.information(self, "Result", msg)


# =========================
# Run App
# =========================
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = AzimuthApp()
    window.show()
    sys.exit(app.exec())
