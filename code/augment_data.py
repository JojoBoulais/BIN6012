import scanpy
import pandas as pd
import re
import os
from pathlib import Path
import jax
import jax.numpy as jnp

data_fp = os.path.join(Path(__file__).parent, "data")

fp_anndata = os.path.join(data_fp, "myocardial_infarction.h5ad")
anndata = scanpy.read_h5ad(fp_anndata)

anndata_csm = anndata[anndata.obs["cell_type"] == "cardiac muscle myoblast", :].copy()

max_cells = 10000
anndata_csm = anndata_csm[:max_cells, :].copy()

# print(anndata_csm.uns["organism_ontology_term_id"])
# exit()

#Shuffle propre (garde obs aligné)
# key = jax.random.key(0)
# perm = jax.random.permutation(key, anndata_csm.n_obs)
# anndata_csm = anndata_csm.X[perm, :]

# Subsample


# anndata_csm.obs["species"] =  "Homo sapiens"
# anndata_csm.obs["organism"] =  "Homo sapiens"

anndata_csm.obs["organism_ontology_term_id"] = "NCBITaxon:9606"

print(anndata_csm)

# Save
anndata_csm.write_h5ad(os.path.join(data_fp, "myocardial_infarction_subset.h5ad"))
