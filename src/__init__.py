"""
H5AD to MatriViz Converter

A Python tool for converting h5ad files (AnnData format) to MatriViz-compatible parquet format.
"""

__version__ = "1.0.0"
__author__ = "UIC PHD"

from .file_to_parquet import AzimuthApp, ConverterThread

__all__ = ["AzimuthApp", "ConverterThread"]