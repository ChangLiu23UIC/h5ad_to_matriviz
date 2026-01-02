import sys, os, json
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import anndata as ad
from pathlib import Path
from PyQt6.QtWidgets import (
    QApplication, QWidget, QLabel, QPushButton, QLineEdit, QFileDialog,
    QVBoxLayout, QHBoxLayout, QMessageBox, QProgressBar, QCheckBox
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal


# =========================================================
# Multi-file Converter Thread
# =========================================================
class ConverterThread(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(bool, str)

    def __init__(self, input_files, tissue, output_dir, json_template=None, downsample=True):
        super().__init__()
        self.input_files = input_files
        self.tissue = tissue
        self.output_dir = output_dir
        self.json_template = json_template
        self.downsample = downsample

    def run(self):
        results = []
        skipped_files = []
        try:
            os.makedirs(self.output_dir, exist_ok=True)

            for idx, input_file in enumerate(self.input_files):
                basename = os.path.splitext(os.path.basename(input_file))[0]
                file_folder = os.path.join(self.output_dir, basename)

                self.progress.emit(f"[{idx + 1}/{len(self.input_files)}] Checking {basename} ...")

                try:
                    # Logic to check and process
                    msg = self.process_single_file(input_file, file_folder)
                    results.append(msg)
                except ValueError as ve:
                    # Catch the "Incompatible" files
                    skipped_files.append(f"{input_file} : {str(ve)}")
                    continue

            # Write skipped files to text file
            if skipped_files:
                skip_log_path = os.path.join(self.output_dir, "skipped_files.txt")
                with open(skip_log_path, "w") as f:
                    f.write("\n".join(skipped_files))
                results.append(f"⚠️ {len(skipped_files)} files skipped. See skipped_files.txt")

            final_msg = "\n\n".join(results)
            self.finished.emit(True, "Process Completed!\n\n" + final_msg)

        except Exception as e:
            import traceback
            tb = traceback.format_exc()
            self.finished.emit(False, f"❌ Fatal Error: {str(e)}\n\nTraceback:\n{tb}")

    def process_single_file(self, input_file, output_dir):
        # --- Load DATA ---
        adata = ad.read_h5ad(input_file)

        # CHECK REQUIREMENTS: Must have X_umap_proj and specific metadata
        required_obs = ['azimuth_label', 'predicted_label']
        has_umap_proj = "X_umap_proj" in adata.obsm
        has_metadata = all(col in adata.obs.columns for col in required_obs)

        if not (has_umap_proj and has_metadata):
            raise ValueError("Missing 'X_umap_proj' or Azimuth metadata (azimuth_label/predicted_label)")

        os.makedirs(output_dir, exist_ok=True)
        self.progress.emit(f"Extracting data for {os.path.basename(input_file)}...")

        # --- Expression Matrix ---
        expr = adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X
        expr_df = pd.DataFrame(expr, index=adata.obs_names, columns=adata.var_names)
        expr_df.columns = expr_df.columns.astype(str)
        expr_df["All_Genes"] = expr_df.mean(axis=1)

        # --- UMAP (Using X_umap_proj as requested) ---
        umap_df = pd.DataFrame(
            adata.obsm["X_umap_proj"],
            columns=["umap_1", "umap_2"],
            index=adata.obs_names
        )

        # --- Metadata (Extract azimuth_label and predicted_label) ---
        meta_df = adata.obs[['azimuth_label', 'predicted_label']].copy()
        meta_df["index"] = meta_df.index.astype(str)

        # --- Merge ---
        merged = umap_df.join(expr_df, how="inner")
        merged = merged.join(meta_df.drop(columns="index"), how="inner")
        merged = merged.reset_index().rename(columns={"index": "index"})

        # Ensure 'index', 'umap_1', 'umap_2' are at the front
        cols = ["index", "umap_1", "umap_2", "azimuth_label", "predicted_label"]
        other_cols = [c for c in merged.columns if c not in cols]
        merged = merged[cols + other_cols]

        # --- DOWNSAMPLING LOGIC ---
        if self.downsample:
            SCREEN_W, SCREEN_H = 600, 500
            x_min, x_max = merged["umap_1"].min(), merged["umap_1"].max()
            y_min, y_max = merged["umap_2"].min(), merged["umap_2"].max()

            if x_max == x_min: x_max += 1
            if y_max == y_min: y_max += 1

            merged["px"] = ((merged["umap_1"] - x_min) / (x_max - x_min) * (SCREEN_W - 1)).astype(int)
            merged["py"] = ((merged["umap_2"] - y_min) / (y_max - y_min) * (SCREEN_H - 1)).astype(int)
            merged = merged.groupby(["px", "py"]).sample(1, random_state=0).drop(columns=["px", "py"])

        # --- Save Files ---
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        expr_path = Path(output_dir) / f"{base_name}_v1.parquet"
        centroid_path = Path(output_dir) / f"{base_name}_centroid_v1.parquet"
        category_path = Path(output_dir) / f"{base_name}_v1_category.json"
        overview_path = Path(output_dir) / f"{base_name}.json"

        # Write Parquet
        pq.write_table(pa.Table.from_pandas(merged), expr_path, compression="zstd")

        # Centroids (Based on predicted_label)
        centroids = (
            merged.groupby("predicted_label")[["umap_1", "umap_2"]]
            .median()
            .reset_index()
            .rename(columns={"predicted_label": "Type", "umap_1": "cen_x", "umap_2": "cen_y"})
        )
        pq.write_table(pa.Table.from_pandas(centroids), centroid_path, compression="zstd")

        # Category JSON (Top 200 genes)
        genes = list(adata.var_names.astype(str))
        category_json = {"All_Genes": genes[:200]}
        with open(category_path, "w") as f:
            json.dump(category_json, f, indent=4)

        # Overview JSON
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
# GUI (Same structure as original)
# =========================================================
class AzimuthApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Azimuth → MatriViz (UMAP_Proj Mode)")
        self.setFixedSize(520, 450)

        layout = QVBoxLayout()
        self.input_line = QLineEdit()
        self.input_btn = QPushButton("Browse Files")
        self.input_btn.clicked.connect(self.select_input)
        hl1 = QHBoxLayout();
        hl1.addWidget(self.input_line);
        hl1.addWidget(self.input_btn)

        self.tissue_line = QLineEdit()
        self.output_line = QLineEdit()
        self.output_btn = QPushButton("Browse Folder")
        self.output_btn.clicked.connect(self.select_output)
        hl2 = QHBoxLayout();
        hl2.addWidget(self.output_line);
        hl2.addWidget(self.output_btn)

        self.json_line = QLineEdit()
        self.json_btn = QPushButton("Browse JSON")
        self.json_btn.clicked.connect(self.select_json)
        hl3 = QHBoxLayout();
        hl3.addWidget(self.json_line);
        hl3.addWidget(self.json_btn)

        self.downsample_cb = QCheckBox("Enable Downsampling (reduces file size)")
        self.downsample_cb.setChecked(True)

        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        self.convert_btn = QPushButton("Convert")
        self.convert_btn.clicked.connect(self.start_conversion)

        layout.addWidget(QLabel("Input .h5ad files:"))
        layout.addLayout(hl1)
        layout.addWidget(QLabel("Tissue name:"))
        layout.addWidget(self.tissue_line)
        layout.addWidget(QLabel("Output directory:"))
        layout.addLayout(hl2)
        layout.addWidget(QLabel("Category template (Optional):"))
        layout.addLayout(hl3)
        layout.addWidget(self.downsample_cb)
        layout.addWidget(self.convert_btn)
        layout.addWidget(self.progress_bar)
        self.setLayout(layout)

    def select_input(self):
        files, _ = QFileDialog.getOpenFileNames(self, "Select Files", "", "H5AD (*.h5ad)")
        if files: self.input_line.setText(";".join(files))

    def select_output(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Folder")
        if folder: self.output_line.setText(folder)

    def select_json(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select JSON", "", "JSON (*.json)")
        if file: self.json_line.setText(file)

    def start_conversion(self):
        files = [f.strip() for f in self.input_line.text().split(";") if f.strip()]
        tissue = self.tissue_line.text().strip()
        outdir = self.output_line.text().strip()
        if not files or not tissue or not outdir:
            QMessageBox.warning(self, "Error", "Fill all required fields.")
            return

        self.convert_btn.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)

        self.worker = ConverterThread(files, tissue, outdir, self.json_line.text(), self.downsample_cb.isChecked())
        self.worker.progress.connect(lambda m: self.progress_bar.setFormat(m))
        self.worker.finished.connect(self.finish_conversion)
        self.worker.start()

    def finish_conversion(self, success, msg):
        self.convert_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        QMessageBox.information(self, "Done", msg)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = AzimuthApp()
    win.show()
    sys.exit(app.exec())