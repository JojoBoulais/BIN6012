import os

from lightning.pytorch import Trainer
from scprint2 import scPRINT2
from scdataloader import utils, Preprocessor, DataModule
import torch
from pathlib import Path
import numpy as np
from utils_filepath import *

curr_dir = Path(__file__).parent

# setup a datamodule to train scprint2 from scratch
# preprocess datasets

def custom_pre_process(adata, 
                       max_cells=10000,
                       ontology_term_id="NCBITaxon:9606", 
                       cell_type="cardiac muscle myoblast"):

    if max_cells:
        adata = adata[:max_cells, :].copy()

    if cell_type:
        adata = adata[adata.obs["cell_type"] == cell_type, :].copy()

    if ontology_term_id:
        adata.obs["organism_ontology_term_id"] = ontology_term_id
        adata = adata.copy()

    return adata

preprocessor = Preprocessor(
    do_postp=False,
    force_preprocess=True,
    additional_preprocess=custom_pre_process
)

adata = preprocessor(adata)

# art = ln.Artifact(adata, description="test")
# art.save()
# ln.Collection(art, key="test", description="test").save()

datamodule = DataModule(
    collection_name="test",
    organisms=["NCBITaxon:9606"], #organism that we will work on
    how="most expr", # for the collator (most expr genes only will be selected)
    max_len=1000, # only the 1000 most expressed
    batch_size=64,
    num_workers=1,
    validation_split=0.1,
)

# setup a model parameter
model_fp = os.path.join(curr_dir, "small-v2.ckpt")

model = scPRINT2.load_from_checkpoint(
    model_fp,
    precpt_gene_emb=None,
    gene_pos_file=None,
)
if not torch.cuda.is_available():
    model = model.to(torch.float32)

model = model.to("cuda" if torch.cuda.is_available() else "cpu")


# to train / fit / test the model setup a trainer
trainer = Trainer(...)

# call the fit function
trainer.fit(model, datamodule=datamodule)


# to do predictions Denoiser, Embedder, GNInfer
# denoiser = Denoiser(...)
# adata = sc.read_h5ad(...)
# denoiser(model, adata=adata)

grn = np.random.rand(10,10)
adata = anndata.AnnData(X=np.random.rand(10,10))

grn = GRNAnnData(adata, grn=grn)

print(grn) #shows the number of elements
grn.varp['GRN'] or grn.grn #shows the GRN
subgrn = grn.get(['gene1', 'gene2']) #only gets some elements from the GRN
subgrn.targets #shows the target connections
subgrn.plot() # displays the network

subgrn.write_h5ad('grn.h5ad') #writes it
read_h5ad('grn.h5ad') #reads it