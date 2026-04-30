from argparse import ArgumentParser
from ast import arg
from pathlib import Path
import shutil
import uuid
from pathlib import Path
import os

import pandas as pd
import scanpy
from grnndata import read_h5ad
from grnndata import utils as grnutils
import matplotlib.pyplot as plt
import seaborn as sns
from anndata.utils import make_index_unique
import numpy as np

from customize_libs import customize_GRNAnnData, BenGRN, compute_pr, GNInfer
from utils_filepath import *
from utils import load_model_with_cuda_if_avail
from Pathways import Pathway, PathwayLoader


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

    # Processing only genes of interest
    pathways = PathwayLoader(argv.pathways_version).load_all_pathways()

    # Loading ground truths
    fp_collectri = get_filepath("data", "collectri.parquet")
    df_collectri = pd.read_parquet(fp_collectri)

    fp_omnipath = get_filepath("data", "omnipath.parquet")
    df_omnipath = pd.read_parquet(fp_omnipath)

    # Forcing human organism
    model.organisms = adata.var["organism"].unique()

    for p, pathway in enumerate(pathways):

        print(
            f"Processing pathway {p+1}/{len(pathways)+1}. ({pathway.get_name_as_path()})")

        genes = pathway.get_genes_ensembl()

        cell_type = str(argv.cell_type).replace(" ", "_")

        # CIs = {}

        cell_diseases = ["normal", "myocardial infarction"]
        for disease in cell_diseases:

            print("Creating GNInfer")
            grn_inferer = GNInfer(
                how="some",
                preprocess="softmax",
                head_agg="mean",
                filtration="none",
                genelist=genes,
                num_genes=0,
                max_cells=0,  # all cells
                num_workers=8,
                batch_size=16,
                precomp_attn=False,
                drop_unexpressed=False
            )

            adata_d = adata[adata.obs["disease"] == disease, :].copy()

            disease = str(disease).replace(" ", "_")

            # tmp path
            tmp_loc = os.path.join(get_filepath(
                "tmp"), str(uuid.uuid4())) + "/"
            tmp_filename = os.path.join(tmp_loc, "grn_fromscprint.h5ad")

            os.makedirs(tmp_loc, exist_ok=True)

            grn_inferer.loc = tmp_loc

            print("processing grn_inferer...")
            grn_for_plot = grn_inferer(model,
                                       adata_d,
                                       cell_type=argv.cell_type)

            # shutil.move(tmp_filename, grn_out_filename)
            os.remove(tmp_filename)
            os.rmdir(tmp_loc)

            grn_for_plot = customize_GRNAnnData(grn_for_plot)

            # Results filepaths
            results_fp = get_filepath("results")
            results_fp = os.path.join(results_fp, pathway.get_name_as_path())
            out_filepath = "_".join([str(Path(adata_fp).stem),
                                     cell_type,
                                     disease])

            out_filepath = out_filepath.replace("/", "_")

            if not os.path.exists(results_fp):
                os.makedirs(results_fp)

            out_filepath_graph = out_filepath + "_graph.png"

            out_filepath_graph = os.path.join(results_fp, out_filepath_graph)

            # ================  Genes network ================
            G, mat = grn_for_plot.plot_subgraph(
                pathway_name=argv.pathway,
                genes=pathway.get_genes_ensembl(),
                out_filepath=out_filepath_graph
            )

            rn = grn_for_plot.rn

            genes_present = [g for g in genes if g in grn_for_plot.grn.index]

            genes_present = [gene for gene in genes_present
                             if gene in rn.keys()]

            mat = grn_for_plot.grn.loc[genes_present,
                                       genes_present].rename(index=rn, columns=rn)

            # ================ heatmap ================
            out_filepath_heatmap = out_filepath + "_heatmap.png"
            out_filepath_heatmap = os.path.join(
                results_fp, out_filepath_heatmap)
            plt.figure(figsize=(15, 15))
            sns.heatmap(mat,
                        annot=True,
                        cmap='viridis',
                        vmin=0,
                        vmax=1.0)
            plt.title(pathway.get_name())
            plt.tight_layout()
            plt.savefig(out_filepath_heatmap)

            # ================ GTs ================
            grn_for_plot.varp["GRN"][np.isnan(grn_for_plot.varp["GRN"])] = 0
            grn_for_plot.var.index = make_index_unique(
                grn_for_plot.var["symbol"].astype(str))

            print("Generating metrics")
            metrics = {}

            # Calculating AURPC et OR
            metrics["collectri"] = BenGRN(
                grn_for_plot).compare_to(to=df_collectri)
            metrics["omnipath"] = BenGRN(
                grn_for_plot).compare_to(to=df_omnipath)

            # Rounding
            for gt in metrics:
                for metric in metrics[gt].keys():
                    metrics[gt][metric] = round(metrics[gt][metric], 5)

            out_filepath_metrics = out_filepath + "_metrics.tsv"
            out_filepath_metrics = os.path.join(
                results_fp, out_filepath_metrics)

            df_metrics = pd.DataFrame(metrics)
            df_metrics.to_csv(out_filepath_metrics, sep="\t")


if __name__ == "__main__":

    parser = ArgumentParser(
        prog='scPrint2 GRN on a preprocessed dataset',
        description="Generating a h5ad file from a preprocessed dataset."
    )

    parser.add_argument("model_file")
    parser.add_argument("h5ad_file")
    parser.add_argument("--cell_type",
                        default="cardiac muscle myoblast",
                        choices=["cardiac muscle myoblast",
                                 "cardiac endothelial cell",
                                 "fibroblast of cardiac tissue"])
    parser.add_argument("--pathways_version",
                        default=3)

    main(parser.parse_args())
