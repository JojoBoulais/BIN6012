from utils_filepath import *
import scanpy

fp = get_data_filepath("myocardial_infarction_preprocessed_2.h5ad")

print(fp)

anndata = scanpy.read_h5ad(fp)

print(anndata)

print("*" * 20)

print(anndata.var["organism"].unique())
print(anndata.var["organism_id"])


print(anndata.var_names[:20].tolist())
print(anndata.var.columns.tolist())
print(anndata.var_names.is_unique)