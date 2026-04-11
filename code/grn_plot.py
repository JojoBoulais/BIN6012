from argparse import ArgumentParser
import os
from pathlib import Path
from bengrn import get_GT_db, get_perturb_gt, get_sroy_gt

from anndata.utils import make_index_unique
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from grnndata import read_h5ad
from scdataloader.utils import load_genes
import bionty as bt
import numpy as np
from scipy.stats import ttest_1samp
from grnndata import from_adata_and_longform

from utils_filepath import *
from customize_libs import customize_GRNAnnData, BenGRN
from Pathways import Pathway

def main(argv):

    diseases = ["normal", "myocardial infarction"]

    for d, in_filepath in enumerate([argv.normal_h5ad, argv.disease_h5ad]):
    
        disease = diseases[d]
        
        print(disease)

        if os.path.exists(in_filepath):
            adata_filepath = in_filepath
        else:
            adata_filepath = get_filepath("grn", in_filepath)
            if not os.path.exists(adata_filepath):
                print(f"{adata_filepath} does not exists.")
                exit(1)

        pathways = Pathway.pathways if argv.pathway == "all" else argv.pathway

        # Loading ground truths
        fp_collectri = get_filepath("data", "collectri.parquet")
        df_collectri = pd.read_parquet(fp_collectri)

        fp_omnipath = get_filepath("data", "omnipath.parquet")
        df_omnipath = pd.read_parquet(fp_omnipath)

        for p in list(pathways):

            print(f"Loading {adata_filepath}...")
            grn_for_plot = read_h5ad(adata_filepath)

            # Adding modified plot_subgraph method here
            grn_for_plot = customize_GRNAnnData(grn_for_plot)

            grn_for_plot.var.index.name = "index"

            # Importing Pathway
            pathway = Pathway(p)

            results_fp = get_filepath("results")
            results_fp = os.path.join(results_fp, pathway.get_name_as_path())
            out_filepath = "_".join([str(Path(adata_filepath).stem), disease, p])

            out_filepath = out_filepath.replace("/", "-")

            if not os.path.exists(results_fp):
                os.makedirs(results_fp)

            out_filepath_graph = out_filepath + "_graph" +".png"

            out_filepath_graph = os.path.join(results_fp, out_filepath_graph)

            # ================  Genes network ================
            # G, mat = grn_for_plot.plot_subgraph(
            #         pathway_name=argv.pathway,
            #         genes=pathway.get_genes_ensembl(),
            #         out_filepath=""
            #     )
            
            rn = grn_for_plot.rn

            genes = pathway.get_genes_ensembl()

            grn_for_plot.available_genes = [gene for gene in genes \
                                            if gene in rn.keys()]

            diff = list(set(genes).difference(set(grn_for_plot.available_genes)))

            if diff:
                print(f"Leaving behind {len(diff)} genes:")
                print(diff)
                genes = grn_for_plot.available_genes

            # Sous-matrice du GRN pour les gènes du pathway
            mask = grn_for_plot.var_names.isin(genes)

            grn_for_plot = grn_for_plot[:, mask].copy()

            # Rename with gene symbol instead of ensembl IDs
            # grn_for_plot.varp["GRN"] = grn_for_plot.varp["GRN"].rename(index=rn, columns=rn)

            mat = grn_for_plot.grn.loc[genes, genes].rename(index=rn, columns=rn)

            # ================ heatmap ================
            out_filepath_heatmap = out_filepath + "_heatmap" +".png"
            out_filepath_heatmap = os.path.join(results_fp, out_filepath_heatmap)
            plt.figure(figsize=(15, 15))
            sns.heatmap(mat, annot=True, cmap='viridis')
            plt.title(argv.pathway)
            plt.tight_layout()
            plt.savefig(out_filepath_heatmap)

            # ================ GTs ================

            grn_for_plot.varp["GRN"][np.isnan(grn_for_plot.varp["GRN"])] = 0
            grn_for_plot.var.index = make_index_unique(grn_for_plot.var["symbol"].astype(str))

            print("Generating metrics")
            metrics = {}

            # Calculating AURPC et OR
            metrics["collectri"] = BenGRN(grn_for_plot).compare_to(to=df_collectri)
            metrics["omnipath"] = BenGRN(grn_for_plot).compare_to(to=df_omnipath)

            df_metrics = pd.DataFrame(metrics)

            # Rounding
            for gt in df_metrics.keys():
                for metric in metrics[gt].keys():
                    metrics[gt][metric] = round(metrics[gt][metric], 4)

            out_filepath_metrics = out_filepath + "_metrics.tsv"
            out_filepath_metrics = os.path.join(results_fp, out_filepath_metrics)

            df_metrics.to_csv(out_filepath_metrics, sep="\t")


if __name__ == "__main__":

    parser = ArgumentParser(
        prog='Plotting GRN results',
        description="Plotting GRN results."
        )
    
    parser.add_argument("--normal_h5ad", type=str, required=True)
    parser.add_argument("--disease_h5ad", type=str, required=True)
    parser.add_argument("--pathway", type=str, required=True, choices=Pathway.pathways + ["all"])

    main(parser.parse_args())
