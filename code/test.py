from utils_filepath import *
import scanpy

fp = get_data_filepath("myocardial_infarction.h5ad")

print(fp)


anndata = scanpy.read_h5ad(fp)

print(anndata)