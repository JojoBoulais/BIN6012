from utils_filepath import *
from _utils import *

import scanpy
import torch

dataset_fp = get_filepath("data", "myocardial_infarction_preprocessed_4.h5ad")
model_fp = get_filepath("model", "small-v2.ckpt")

model = load_model_with_cuda_if_avail(model_fp)
adata = scanpy.read_h5ad(dataset_fp)

print(model.organisms)

topk = list(set(adata.var.index[adata.var.highly_variable]))

print("topk[:20]")
print(topk[:20])



for gene in ["AL359922", "AL359922.1", "ENST00000582560",
             "ENST00000582560.1", "ENSG00000266732", "ENSG00000266732.1",
             "LPA", "ENSG00000198670", "ENST00000316300.10", "ENSG00000112137", "ENST00000676159.1"]:

    print(gene in model.genes)

# model.organisms = adata.var["organism"].unique()
# torch.save(model, get_model_filepath("small-v2_organisms_subset.ckpt"))