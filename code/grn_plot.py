from argparse import ArgumentParser
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt

import os

from anndata.utils import make_index_unique


from grnndata import read_h5ad

from utils_filepath import *
from customize_libs import customize_GRNAnnData
from Pathways import Pathway

def main(argv):

    if os.path.exists(argv.h5ad_file):
        adata_filepath = argv.h5ad_file
    else:
        adata_filepath = get_filepath("grn", argv.h5ad_file)
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

    out_filepath = "_".join([str(Path(adata_filepath).stem), argv.pathway])

    out_filepath = out_filepath.replace("/", "-")

    if not os.path.exists(results_fp):
        os.makedirs(results_fp)

    out_filepath_graph = out_filepath + "_graph" +".png"

    out_filepath_graph = os.path.join(results_fp, out_filepath_graph)

    # Plotting genes network
    adjacency_data = grn_for_plot.plot_subgraph(
            pathway_name=argv.pathway,
            genes=pathway.get_genes_ensembl(),
            out_filepath=out_filepath_graph
        )
    
    # Plotting heatmap
    out_filepath_heatmap = out_filepath + "_heatmap" +".png"
    out_filepath_heatmap = os.path.join(results_fp, out_filepath_heatmap)
    sns.heatmap(adjacency_data, annot=True, cmap='viridis')
    plt.savefig(out_filepath_heatmap)

if __name__ == "__main__":

    parser = ArgumentParser(
        prog='Plotting GRN results',
        description="Plotting GRN results."
        )
    
    parser.add_argument("h5ad_file")
    parser.add_argument("--pathway", type=str, required=True, choices=Pathway.pathways)
    parser.add_argument("--threshold", type=float, default=55.0, required=False)
    parser.add_argument("--max_genes", type=int, default=15, required=False)

    main(parser.parse_args())
