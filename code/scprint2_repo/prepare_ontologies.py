from scdataloader.utils import _adding_scbasecamp_genes, populate_my_ontology, load_genes
import os.path
import urllib.request
import torch

# %load_ext autoreload
# %autoreload 2


populate_my_ontology(
    organisms_clade=["vertebrates"],
    sex=["PATO:0000384", "PATO:0000383"],
    organisms=["NCBITaxon:10090", "NCBITaxon:9606"],
    # celltypes=None,
    # ethnicities=None,
    # assays=None,
    # tissues=None,
    # diseases=None,
    # dev_stages=None,
)
_adding_scbasecamp_genes()