from argparse import ArgumentParser
from ast import arg
from pathlib import Path
import shutil
import uuid
from pathlib import Path
import os

from customize_libs import GNInfer
from lightning.pytorch import Trainer

import scanpy
from grnndata import read_h5ad
from grnndata import utils as grnutils

from utils_filepath import *
from utils import load_model_with_cuda_if_avail
from Pathways import Pathway


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

    print(adata.obs["disease"].unique())

    print(len((adata.obs["disease"] == "myocardial infarction")
                & (adata.obs["cell_type"] == "cardiac muscle myoblast")))

    # Processing only genes of interest
    pathways = Pathway.load_all_pathways()
    gene_list = []
    for pathway in pathways:
        gene_list.extend(pathway.get_genes_ensembl())

    gene_list = list(set(gene_list))

    print(f"{len(gene_list)} genes included in GNinfer.")

    # Forcing human organism
    model.organisms = adata.var["organism"].unique()

    print("Creating GNInfer")
    grn_inferer = GNInfer(
        how="most var across", # most var across
        preprocess="softmax",
        head_agg="mean",
        filtration="none",
        genelist=gene_list, #gene_list, # processing only these gene with how="some"
        num_genes=2000,
        max_cells=0, # all cells
        num_workers=8,
        batch_size=16,
        precomp_attn=False # <============ was True before
    )

    cell_type = str(argv.cell_type).replace(" ", "_")

    cell_diseases = ["normal", "myocardial infarction"]
    for disease in cell_diseases:

        # if disease != "all":
        #     adata = adata[adata.obs["disease"] == disease, :].copy()

        disease = str(disease).replace(" ", "_")

        grn_out_filename = "_".join([str(Path(adata_fp).stem),
                                    str(Path(model_fp).stem),
                                    cell_type,
                                    disease])
        
        # ====================================> test
        grn_out_filename = grn_out_filename + "_test.h5ad"

        grn_out_filename = get_filepath("grn", grn_out_filename)

        # tmp path
        tmp_loc = os.path.join(get_filepath("tmp"), str(uuid.uuid4())) + "/"
        tmp_filename = os.path.join(tmp_loc, "grn_fromscprint.h5ad")

        os.makedirs(tmp_loc, exist_ok=True)

        grn_inferer.loc = tmp_loc

        print("processing grn_inferer...")
        grn_inferer(model, 
                    adata,
                    cell_type=argv.cell_type)

        shutil.move(tmp_filename, grn_out_filename)
        os.rmdir(tmp_loc)
        print(f"Generated GRN file : {grn_out_filename}.")


if __name__ == "__main__":

    parser = ArgumentParser(
        prog='scPrint2 GRN on a preprocessed dataset',
        description="Generating a h5ad file from a preprocessed dataset."
        )
    
    parser.add_argument("model_file")
    parser.add_argument("h5ad_file")
    parser.add_argument("--cell_type", default="cardiac muscle myoblast")

    main(parser.parse_args())
