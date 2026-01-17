from collections import Counter
from src.data.gdc_expression import list_luad_expression_files

files = list_luad_expression_files(workflow_type=None, sample_type="Primary Tumor", max_files=5000)
print("Total files:", len(files))
print("Top workflow_type:")
from collections import Counter
c = Counter([f.workflow_type for f in files])
for k, v in c.most_common(20):
    print(v, k)