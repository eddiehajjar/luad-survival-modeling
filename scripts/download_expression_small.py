from pathlib import Path
from src.data.gdc_expression import list_luad_expression_files, download_expression_files

files = list_luad_expression_files(workflow_type="STAR - Counts", sample_type="Primary Tumor")
files = files[:25]

out_dir = Path("data/raw/expression/HTSeq-FPKM")
download_expression_files(files, out_dir=out_dir)

print("Downloaded:", len(list(out_dir.glob("*.tsv"))), "files into", out_dir)