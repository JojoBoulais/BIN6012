from argparse import ArgumentParser
from pathlib import Path
import shutil
import uuid
from pathlib import Path
import os

from anndata.utils import make_index_unique
from scprint2.tasks import GNInfer
from lightning.pytorch import Trainer
from scprint2 import scPRINT2
import numpy as np
import scanpy
from grnndata import read_h5ad
from grnndata import utils as grnutils

from utils_filepath import *
from _utils import load_model_with_cuda_if_avail


def main(argv):

    if os.path.exists(argv.model_file):
        model_fp = argv.model_file
    else:
        model_fp = get_filepath("model", argv.model_file)
        if not os.path.exists(model_fp):
            print(f"{model_fp} does not exists.")
            exit(1)

    if os.path.exists(argv.h5ad_file):
        adata_fp = argv.h5ad_file
    else:
        adata_fp = get_filepath("data", argv.h5ad_file)
        if not os.path.exists(adata_fp):
            print(f"{adata_fp} does not exists.")
            exit(1)

    print("loading model...")
    model = load_model_with_cuda_if_avail(model_fp)

    print("loading dataset.")
    adata = scanpy.read_h5ad(adata_fp)

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


    
    # - "most var across": select the most variable genes across all cell types
    # - "most var within": select the most variable genes within a cell type
    # - "random expr": select random expressed genes
    # - "some": select a subset of genes defined in genelist
    # - "most expr": select the most expressed genes in the cell type

    # This is generating file 'grn_fromscprint.h5ad' to loc
    # sub = adata[adata.obs.celltype == argv.cell_type].copy()
    topk = list(set(adata.var.index[adata.var.highly_variable])) + ["ENSG00000198670",  "ENSG00000112137"]

    #['ENSG00000165325', 'ENSG00000259124', 'ENSG00000234840', 'ENSG00000263708', 'ENSG00000188177', 'ENSG00000134909', 'ENSG00000213316', 'ENSG00000064692', 'ENSG00000229599', 'ENSG00000198626', 'ENSG00000204103', 'ENSG00000099725', 'ENSG00000145675', 'ENSG00000053524', 'ENSG00000137486', 'ENSG00000227531', 'ENSG00000165244', 'ENSG00000225913', 'ENSG00000286861', 'ENSG00000163535']

    # is_expr = set(
    #     adata.var.index[
    #         np.array(adata[adata.obs.celltype == argv.cell_type].X.sum(0) > 0)[0]
    #     ]
    # )

    # genes = list(is_expr & set(model.genes) & topk)
    genes = topk

    # we now compute eigen centrality creating a sparse network by only keeping the top 20 neighbors for each gene in the network
    # TOP = 20
    # grnutils.get_centrality(grn_h, TOP, top_k_to_disp=0)

    model.organisms = adata.var["organism"].unique()

    print("Creating GNInfer")
    grn_inferer = GNInfer(
        how="most var across",
        preprocess="softmax",
        head_agg="mean",
        filtration="none",
        # genelist=genes,
        num_genes=4000,
        max_cells=1024,
        num_workers=8,
        batch_size=16,
        precomp_attn=True
    )

    print("processing grn_inferer...")
    grn = grn_inferer(model, adata, cell_type=argv.cell_type) 

    if argv.out_filename:
        grn_out_filename = argv.out_filename
    else:
        grn_out_filename = Path(adata_fp).stem + Path(model_fp) + ".h5ad"

    grn_out_filename = get_filepath("grn", grn_out_filename)

    # tmp path
    tmp_loc = os.path.join(get_filepath("tmp"), uuid.uuid4())
    tmp_filename = os.path.join(tmp_loc, "grn_fromscprint.h5ad")
    grn.loc = tmp_loc

    shutil.move(tmp_filename, grn_out_filename)
    os.rmdir(tmp_loc)
    print(f"Generated GRN file : {grn_out_filename}.")

    grn.var.index = make_index_unique(grn.var["symbol"].astype(str))
    grn.var.index.name = "index"
    grn.varp["GRN"][np.isnan(grn.varp["GRN"])] = 0

    gene_of_interest = "ENSG00000198670"

    print(f"plotting using {gene_of_interest} gene.")
    grn.plot_subgraph(
        seed=gene_of_interest, only=55, interactive=False, max_genes=15,
    )

if __name__ == "__main__":

    parser = ArgumentParser(
        prog='scPrint2 GRN on a preprocessed dataset',
        description="Generating a h5ad file from a preprocessed dataset."
        )
    
    parser.add_argument("model_file")
    parser.add_argument("h5ad_file")
    parser.add_argument("--out_filename", default="", required=False)
    parser.add_argument("--cell_type", default="cardiac muscle myoblast")

    main(parser.parse_args())