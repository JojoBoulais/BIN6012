from pathlib import Path
import os

import torch
from scprint2 import scPRINT2
import scanpy

# ================ Load data ================
fp_anndata = os.path.join(Path(__file__).parent, "data", "alzheimer.h5ad")
anndata = scanpy.read_h5ad(fp_anndata)

# ================ Load model ================

model_fp = 'small-v2.ckpt'

model = scPRINT2.load_from_checkpoint(
    model_fp,
    precpt_gene_emb=None,
    gene_pos_file=None,
)
if not torch.cuda.is_available():
    model = model.to(torch.float32)

model = model.to("cuda" if torch.cuda.is_available() else "cpu")

# ================ finetune model ================

