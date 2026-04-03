from argparse import ArgumentParser

import os

from networkx import adjacency_data
import scanpy
from anndata.utils import make_index_unique

from grnndata import utils as grnutils
from grnndata import read_h5ad

from utils_filepath import *


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

    grn_for_plot.var.index = make_index_unique(grn_for_plot.var["symbol"].astype(str))
    grn_for_plot.var.index.name = "index"
    # grn_for_plot.varp["GRN"][np.isnan(grn_for_plot.varp["GRN"])] = 0

    print(f"plotting using {argv.gene_of_interest} gene...")

    print(grn_for_plot.var["symbol"])
    
    for gene in ["PHACTR1", "LPA", "AL359922.1", "AL359922", "ATXN2", "CELSR2", "MIA3", "NECTIN2"]:

        print(gene, gene in grn_for_plot.var["symbol"])


    

    adjacency_data = grn_for_plot.plot_subgraph(
            seed=argv.gene_of_interest, only=55, interactive=True, max_genes=argv.max_genes,
        )



if __name__ == "__main__":

    parser = ArgumentParser(
        prog='Plotting GRN results',
        description="Plotting GRN results."
        )
    
    parser.add_argument("h5ad_file")
    parser.add_argument("--gene_of_interest", type=str, required=True)
    parser.add_argument("--max_genes", type=int, default=15, required=False)
    parser.add_argument("--out_filename", default="", required=False)

    main(parser.parse_args())
