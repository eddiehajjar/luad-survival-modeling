from pathlib import Path
from src.data.gdc_expression import list_luad_expression_files

files = list_luad_expression_files(workflow_type="STAR - Counts", sample_type="Primary Tumor")
print("Total expression files returned:", len(files))
print("First 5:")
for f in files[:5]:
    print(f)