from argparse import ArgumentParser
from scprint2.tasks import GNInfer
from anndata.utils import make_index_unique

from lightning.pytorch import Trainer
from scprint2 import scPRINT2
from pathlib import Path
import numpy as np
from _utils import load_model_with_cuda_if_avail
from utils_filepath import *
import scanpy
from grnndata import read_h5ad

def main(argv):

    model = load_model_with_cuda_if_avail(get_model_filepath(argv.model_file))

    adata = scanpy.read_h5ad(get_data_filepath(argv.h5ad_file))

    # if argv.cell_type:
    #     sub = adata[adata.obs.celltype == argv.cell_type].copy()
    #     scanpy.pp.log1p(sub)
    #     scanpy.pp.highly_variable_genes(sub, n_top_genes=6000)
    #     topk = set(sub.var.index[sub.var.highly_variable])
    #     is_expr = set(
    #         adata.var.index[
    #             np.array(adata[adata.obs.celltype == argv.cell_type].X.sum(0) > 0)[0]
    #         ]
    #     )
    # else:
    #     scanpy.pp.log1p(adata)
    #     scanpy.pp.highly_variable_genes(adata, n_top_genes=6000)
    #     topk = set(adata.var.index[adata.var.highly_variable])
    #     is_expr = set(
    #         adata.var.index[
    #             np.array(adata[adata.obs.celltype == argv.cell_type].X.sum(0) > 0)[0]
    #         ]
    #     )

    # genes_h = list(is_expr & set(model.genes) & topk)

    # grn_inferer = GNInfer(
    #         how="some",
    #         preprocess="softmax",
    #         head_agg="mean",
    #         filtration="none",
    #         genelist=genes_h,
    #         max_cells=1024,
    #         num_workers=4,
    #         batch_size=8,
    #         precomp_attn=False,
    #         cell_type_col="celltype",
    #     )
    
    # - "most var across": select the most variable genes across all cell types
    # - "most var within": select the most variable genes within a cell type
    # - "random expr": select random expressed genes
    # - "some": select a subset of genes defined in genelist
    # - "most expr": select the most expressed genes in the cell type

    grn_inferer = GNInfer(
        how="most var across",
        preprocess="softmax",
        head_agg="mean",
        filtration="none",
        num_genes=4000,
        max_cells=1024,
        num_workers=8,
        batch_size=16,
        precomp_attn=True
    )

    grn = grn_inferer(model, adata, cell_type=argv.cell_type)
    grn.var.index = make_index_unique(grn.var["symbol"].astype(str))
    grn.var.index.name = "index"
    grn.varp["GRN"][np.isnan(grn.varp["GRN"])] = 0

    grn.plot_subgraph(
        grn.genelist[0], only=55, interactive=False, max_genes=15,
    )

if __name__ == "__main__":

    parser = ArgumentParser(
        prog='Preprocess scPrint2 data',
        description="Preprocessing h5ad file to make it compatible with scPrint2."
        )
    
    parser.add_argument("model_file")
    parser.add_argument("h5ad_file")
    parser.add_argument("--cell_type", default="cardiac muscle myoblast")


    main(parser.parse_args())