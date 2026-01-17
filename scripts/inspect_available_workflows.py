from collections import Counter
from src.data.gdc_expression import _post_json, GDC_FILES_ENDPOINT

def main():
    # Same query, but DO NOT filter on analysis.workflow_type
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-LUAD"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "cases.samples.sample_type", "value": ["Primary Tumor"]}},
        ],
    }

    fields = [
        "file_id",
        "analysis.workflow_type",
        "data_format",
        "experimental_strategy",
    ]

    payload = {"filters": filters, "fields": ",".join(fields), "format": "JSON", "size": 5000}
    out = _post_json(GDC_FILES_ENDPOINT, payload)

    hits = out.get("data", {}).get("hits", [])
    print("Total hits:", len(hits))

    wf = []
    fmts = []
    for h in hits:
        wf.append((h.get("analysis") or {}).get("workflow_type"))
        fmts.append(h.get("data_format"))

    print("\nTop workflow_type values:")
    for k, v in Counter([x for x in wf if x]).most_common(20):
        print(f"{v:5d}  {k}")

    print("\nTop data_format values:")
    for k, v in Counter([x for x in fmts if x]).most_common(20):
        print(f"{v:5d}  {k}")

if __name__ == "__main__":
    main()