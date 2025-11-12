import json
from mygene import MyGeneInfo

# ========= Load your JSON =========
with open("matrisome_list.json") as f:
    data = json.load(f)

# ========= Extract all gene symbols =========
def extract_genes(d):
    genes = []
    if isinstance(d, dict):
        for v in d.values():
            genes.extend(extract_genes(v))
    elif isinstance(d, list):
        genes.extend(d)
    return genes

genes = list(set(extract_genes(data)))

# ========= Query MyGeneInfo =========
mg = MyGeneInfo()
results = mg.querymany(genes, scopes="symbol", fields="ensembl.gene", species="human")

# ========= Build mapping dict =========
mapping = {}
for r in results:
    gene = r.get("query")
    ens = None
    if "ensembl" in r:
        val = r["ensembl"]
        if isinstance(val, list):
            ens = val[0].get("gene")
        elif isinstance(val, dict):
            ens = val.get("gene")
    mapping[gene] = ens

# ========= Recursively replace genes with ENSG =========
def replace_genes(d):
    if isinstance(d, dict):
        return {k: replace_genes(v) for k, v in d.items()}
    elif isinstance(d, list):
        return [mapping.get(g, g) if isinstance(g, str) else replace_genes(g) for g in d]
    else:
        return d

updated_data = replace_genes(data)

# ========= Save to JSON =========
with open("matrisome_list_ENSG.json", "w") as f:
    json.dump(updated_data, f, indent=2)

print("✅ Saved: matrisome_list_ENSG.json")
