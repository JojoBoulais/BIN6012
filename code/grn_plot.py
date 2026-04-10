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


        print(f"Loading {adata_filepath}...")
        grn_for_plot = read_h5ad(adata_filepath)

        # Adding modified plot_subgraph method here
        grn_for_plot = customize_GRNAnnData(grn_for_plot)

        grn_for_plot.var.index.name = "index"

        # Importing Pathway
        pathway = Pathway(argv.pathway)

        results_fp = get_filepath("results")

        out_filepath = "_".join([str(Path(adata_filepath).stem), disease, argv.pathway])

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
        grn_for_plot.grn = grn_for_plot.grn.loc[genes, genes]

        # Rename with gene symbol instead of ensembl IDs
        grn_for_plot.grn = grn_for_plot.grn.rename(index=rn, columns=rn)


        # ================ heatmap ================
        out_filepath_heatmap = out_filepath + "_heatmap" +".png"
        out_filepath_heatmap = os.path.join(results_fp, out_filepath_heatmap)
        plt.figure(figsize=(15, 15))
        sns.heatmap(grn_for_plot, annot=True, cmap='viridis')
        plt.title(argv.pathway)
        plt.tight_layout()
        plt.savefig(out_filepath_heatmap)

        # ================ GTs ================

        # “GT ≈ interactions biologiques générales (majoritairement basales/canoniques/sans diseases)”

        grn_for_plot.varp["GRN"][np.isnan(grn_for_plot.varp["GRN"])] = 0
        grn_for_plot.var.index = make_index_unique(grn_for_plot.var["symbol"].astype(str))



        # ==== GWPS ====
        # gwps = get_perturb_gt()
        # gwps.var.index = gwps.var.gene_name
        # gwps = gwps.extract_links()

        # gt_inter_fp = get_filepath("data", "RF2-PPI_scores")

        # df_gt_inter = pd.read_csv(
        #     gt_inter_fp,
        #     sep="\t",
        #     skiprows=7,
        #     header=None,
        #     names=["Pair", "value", "source"],
        # )

        # df_gt_inter["geneA"] = df_gt_inter["Pair"].astype(str).str.split("_").str[0]
        # df_gt_inter["geneB"] = df_gt_inter["Pair"].astype(str).str.split("_").str[1]

        # df_gt_inter = df_gt_inter[["geneA", "geneB", "value"]]

        # Keeping only genes from pathway
        # df_RF2 = df_RF2[
        #     df_RF2["geneA"].isin(grn_for_plot.available_genes)
        #     & df_RF2["geneB"].isin(grn_for_plot.available_genes)
        # ]

        fp_collectri = get_filepath("data", "collectri.parquet")
        df_collectri = pd.read_parquet(fp_collectri)

        fp_omnipath = get_filepath("data", "omnipath.parquet")
        df_omnipath = pd.read_parquet(fp_omnipath)

        print("Generating metrics")
        metrics = {}

        metrics["collectri"] = BenGRN(grn_for_plot).compare_to(to=df_collectri)
        metrics["omnipath"] = BenGRN(grn_for_plot).compare_to(to=df_omnipath)
        # metrics["cellmap"] = BenGRN(grn_for_plot).compare_to(gt_cm)
        # grn_for_plot.varp["GRN"] = grn_for_plot.varp["GRN"].T
        # metrics["gwps"] = BenGRN(grn_for_plot).compare_to(gwps)
        # metrics["interact"] = BenGRN(grn_for_plot).compare_to(df_gt_inter)

        df_metrics = pd.DataFrame(metrics)

        out_filepath_metrics = out_filepath + "_metrics.tsv"
        out_filepath_metrics = os.path.join(results_fp, out_filepath_metrics)

        df_metrics.to_csv(out_filepath_metrics, sep="\t")

        # show = {}
        # for gt, v in metrics.items():
        #     if gt + "_auprc" in show:
        #         show[gt + "_auprc"].append(v["auprc"] / v["rand_precision"])
        #     else:
        #         show[gt + "_auprc"] = [v["auprc"] / v["rand_precision"]]
        #     if gt + "_odd_ratio" in show:
        #         show[gt + "_odd_ratio"].append(v["odd_ratio"])
        #     else:
        #         show[gt + "_odd_ratio"] = [v["odd_ratio"]]


        # df = pd.DataFrame(
        #     {i: v for i, v in show.items() if "_auprc" in i},
        #     index=[
        #         "cardiac muscle myoblast"
        #     ],
        # )

        # df.columns = [col.replace("_auprc", "") for col in df.columns]
        # df_long = df.reset_index().melt(
        #     id_vars="index", var_name="Database", value_name="AUPRC"
        # )

        # df_long = df_long.rename(columns={"index": "Cell Type"})

        # plt.figure(figsize=(10, 6))
        # ax = sns.boxplot(
        #     data=df_long,
        #     x="Database",
        #     y="AUPRC",
        #     # hue="method",
        #     fliersize=0,
        #     palette="Set2",
        #     showmeans=True,
        #     meanprops={"marker": "d", "markerfacecolor": "white", "markeredgecolor": "black"},
        # )
        # sns.stripplot(
        #     data=df_long,
        #     x="Database",
        #     y="AUPRC",
        #     # hue="method",
        #     # palette=["black"],  # Use black for both methods
        #     dodge=True,  # This aligns the strips with the boxplot groups
        #     alpha=0.4,
        #     size=5,
        #     legend=False,  # Avoid duplicate legend entries
        # )
        # # Reference line at 1
        # ax.axhline(1, ls="--", c="gray", lw=1)

        # # Compute and annotate p-values (one-sample test vs 0)
        # max_y = df_long["AUPRC"].max()
        # min_y = df_long["AUPRC"].min()
        # span = max(max_y - min_y, 1e-6)
        # offset = 0.08 * span

        # # We need multiple cell types if we want a p-value

        # # for i, col in enumerate(df.columns):

        # #     x = df[col].dropna().values
        # #     x = x[np.isfinite(x)]
        # #     print()

        # #     print("heeere")
        # #     print(x)

        # #     stat, p = ttest_1samp(x, 1.0, nan_policy="omit")
        # #     y = df[col].max() + offset
        # #     ax.text(i, y, f"p={p:.2e}", ha="left", va="bottom", fontsize=9)

        # ax.set_ylim(min_y - offset, max_y + 3 * offset)
        # plt.title("AUPRC by Database (p-values vs 1)")
        # plt.xlabel("Database")
        # plt.ylabel("AUPRC")
        # plt.tight_layout()

        # out_filepath_gt = out_filepath + "_GTs.png"
        # out_filepath_gt = os.path.join(results_fp, out_filepath_gt)

        # plt.savefig(out_filepath_gt)


if __name__ == "__main__":

    parser = ArgumentParser(
        prog='Plotting GRN results',
        description="Plotting GRN results."
        )
    
    parser.add_argument("--normal_h5ad", type=str, required=True)
    parser.add_argument("--disease_h5ad", type=str, required=True)
    parser.add_argument("--pathway", type=str, required=True, choices=Pathway.pathways)

    main(parser.parse_args())
