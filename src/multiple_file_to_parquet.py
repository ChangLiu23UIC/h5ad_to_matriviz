# --- keep all imports the same as your original code ---
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


# =========================================================
# Multi-file Converter Thread
# =========================================================
class ConverterThread(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(bool, str)

    def __init__(self, input_files, tissue, output_dir, json_template=None):
        super().__init__()
        self.input_files = input_files  # now a list
        self.tissue = tissue
        self.output_dir = output_dir
        self.json_template = json_template

    def run(self):

        results = []
        try:
            for idx, input_file in enumerate(self.input_files):
                basename = os.path.splitext(os.path.basename(input_file))[0]
                file_folder = os.path.join(self.output_dir, basename)

                os.makedirs(file_folder, exist_ok=True)
                self.progress.emit(f"[{idx+1}/{len(self.input_files)}] Processing {basename} ...")

                msg = self.process_single_file(input_file, file_folder)
                results.append(msg)

            final_msg = "\n\n".join(results)
            self.finished.emit(True, "All files processed!\n\n" + final_msg)

        except Exception as e:
            import traceback
            tb = traceback.format_exc()
            self.finished.emit(False, f"❌ Error: {str(e)}\n\nTraceback:\n{tb}")

    def process_single_file(self, input_file, output_dir):
        """
        Exactly the same as your original logic,
        just moved into this function.
        """

        self.progress.emit("Loading data...")

        # --- Load DATA ---
        if input_file.endswith(".h5ad"):
            adata = ad.read_h5ad(input_file)
        else:
            raise ValueError("Only .h5ad supported in multi-file mode")

        # --- Expression Matrix ---
        self.progress.emit("Extracting expression matrix...")

        expr = adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X
        expr_df = pd.DataFrame(expr, index=adata.obs_names, columns=adata.var_names)
        expr_df.columns = expr_df.columns.astype(str)

        # Add All_Genes
        expr_df["All_Genes"] = expr_df.mean(axis=1)

        # UMAP
        if "X_umap" not in adata.obsm:
            raise ValueError("No UMAP found")

        umap_df = pd.DataFrame(
            adata.obsm["X_umap"],
            columns=["umap_1", "umap_2"],
            index=adata.obs_names
        )
        umap_df["index"] = umap_df.index.astype(str)

        # Metadata
        meta_df = adata.obs.copy()
        meta_df["index"] = meta_df.index

        celltype_cols = [c for c in meta_df.columns if "celltype" in c.lower()]
        label_col = celltype_cols[-1] if celltype_cols else None

        # Merge
        self.progress.emit("Merging data...")

        merged = umap_df.join(expr_df, how="inner")
        merged = merged.reset_index().rename(columns={"index": "index"})
        merged['index'] = merged['index'].astype(str)

        merged = merged[["index", "umap_1", "umap_2"] +
                        [c for c in merged.columns if c not in ["index","umap_1","umap_2"]]]

        # --- File naming ---
        base_name = os.path.splitext(os.path.basename(input_file))[0]

        expr_path = Path(output_dir) / f"{base_name}_v1.parquet"
        centroid_path = Path(output_dir) / f"{base_name}_centroid_v1.parquet"
        category_path = Path(output_dir) / f"{base_name}_v1_category.json"
        overview_path = Path(output_dir) / f"{base_name}.json"

        # Write parquet
        self.progress.emit("Writing parquet...")
        pq.write_table(pa.Table.from_pandas(merged), expr_path, compression="zstd")

        # Centroids
        if label_col:
            meta_sub = meta_df[[label_col, "index"]].merge(umap_df, on="index", how="left")
            centroids = (
                meta_sub.groupby(label_col)[["umap_1", "umap_2"]]
                .median()
                .reset_index()
                .rename(columns={label_col: "Type", "umap_1": "cen_x", "umap_2": "cen_y"})
            )
        else:
            centroids = pd.DataFrame(columns=["Type","cen_x","cen_y"])

        pq.write_table(pa.Table.from_pandas(centroids), centroid_path, compression="zstd")

        # Category JSON
        self.progress.emit("Writing category JSON...")
        genes = list(adata.var_names.astype(str))
        category_json = { "All_Genes": genes[:200] }
        with open(category_path, "w") as f:
            json.dump(category_json, f, indent=4)

        # Overview
        self.progress.emit("Writing overview JSON...")
        overview_json = {
            "fileType": "matriviz",
            "version": "0.0.1",
            "category_name": base_name,
            "category_description": self.tissue,
            "parquet_file": os.path.basename(expr_path),
            "category_file": os.path.basename(category_path),
            "centroid_file": os.path.basename(centroid_path)
        }

        with open(overview_path, "w") as f:
            json.dump(overview_json, f, indent=4)

        return f"✔ Finished: {base_name}"


# =========================================================
# GUI
# =========================================================
class AzimuthApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Azimuth → MatriViz Multi-file Converter")
        self.setFixedSize(520, 400)

        layout = QVBoxLayout()

        # Input files
        self.input_label = QLabel("Input .h5ad files (multi-select supported):")
        self.input_line = QLineEdit()
        self.input_btn = QPushButton("Browse")
        self.input_btn.clicked.connect(self.select_input)
        hl1 = QHBoxLayout()
        hl1.addWidget(self.input_line)
        hl1.addWidget(self.input_btn)

        # Tissue
        self.tissue_label = QLabel("Tissue name:")
        self.tissue_line = QLineEdit()

        # Output
        self.output_label = QLabel("Output directory:")
        self.output_line = QLineEdit()
        self.output_btn = QPushButton("Browse")
        self.output_btn.clicked.connect(self.select_output)
        hl2 = QHBoxLayout()
        hl2.addWidget(self.output_line)
        hl2.addWidget(self.output_btn)

        # JSON template (optional)
        self.json_label = QLabel("Category template JSON:")
        self.json_line = QLineEdit()
        self.json_btn = QPushButton("Browse")
        self.json_btn.clicked.connect(self.select_json)
        hl3 = QHBoxLayout()
        hl3.addWidget(self.json_line)
        hl3.addWidget(self.json_btn)

        # Progress
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(False)

        self.convert_btn = QPushButton("Convert")
        self.convert_btn.clicked.connect(self.start_conversion)

        # Layout
        layout.addWidget(self.input_label)
        layout.addLayout(hl1)
        layout.addWidget(self.tissue_label)
        layout.addWidget(self.tissue_line)
        layout.addWidget(self.output_label)
        layout.addLayout(hl2)
        layout.addWidget(self.json_label)
        layout.addLayout(hl3)
        layout.addWidget(self.convert_btn)
        layout.addWidget(self.progress_bar)

        self.setLayout(layout)

    def select_input(self):
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select Input Files", "", "H5AD Files (*.h5ad)"
        )
        if files:
            self.input_line.setText(";".join(files))

    def select_output(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if folder:
            self.output_line.setText(folder)

    def select_json(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select JSON Template", "", "JSON (*.json)")
        if file:
            self.json_line.setText(file)

    def start_conversion(self):
        file_list = self.input_line.text().split(";")
        file_list = [f.strip() for f in file_list if f.strip()]

        tissue = self.tissue_line.text().strip()
        outdir = self.output_line.text().strip()
        json_template = self.json_line.text().strip() or None

        if not file_list or not tissue or not outdir:
            QMessageBox.warning(self, "Missing Info", "Please fill all fields.")
            return

        self.convert_btn.setEnabled(False)
        self.progress_bar.setVisible(True)

        self.worker = ConverterThread(file_list, tissue, outdir, json_template)
        self.worker.progress.connect(self.show_progress)
        self.worker.finished.connect(self.finish_conversion)
        self.worker.start()

    def show_progress(self, msg):
        self.progress_bar.setFormat(msg)

    def finish_conversion(self, success, msg):
        self.convert_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        QMessageBox.information(self, "Result", msg)


# Run
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = AzimuthApp()
    window.show()
    sys.exit(app.exec())
