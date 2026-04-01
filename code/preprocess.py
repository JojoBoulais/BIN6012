from argparse import ArgumentParser
from ast import arg
import os
from pathlib import Path

from scdataloader import utils, Preprocessor, DataModule
from utils_filepath import *
import scanpy

def main(argv):

    if os.path.exists(argv.in_file):
        fp = argv.in_file
    else:
        fp = get_filepath("data", argv.in_file)

        if not os.path.exists(fp):
            print(f"{fp} does not exists")
            exit(1)

    adata = scanpy.read_h5ad(fp)

    if argv.max_cells > 0:
        print(f"Keep only first {str(argv.max_cells)} cells")
        adata = adata[:argv.max_cells, :].copy()
        
    if argv.cell_type:
        print(f"Keep only {argv.cell_type} cells type.")
        adata = adata[adata.obs["cell_type"] == argv.cell_type, :].copy()

    if argv.organism:
        print(f"Injecting {argv.organism} as adata.obs['organism_ontology_term_id'] value.")
        adata.obs["organism_ontology_term_id"] = argv.organism
        adata = adata.copy()

    # Running preprocessor
    print("Running Preprocessor on dataset.")
    preprocessor = Preprocessor(
        do_postp=False,
        force_preprocess=True
    )

    preprocessed_adata = preprocessor(adata.copy())

    # These are required for GRN to run
    if not argv.omit_log_transform:
        print("log transforming dataset.")
        scanpy.pp.log1p(preprocessed_adata)

    if argv.n_top_genes > 0:
        print(f"{str(argv.n_top_genes)} highly_variable_genes.")
        scanpy.pp.highly_variable_genes(preprocessed_adata, n_top_genes=argv.n_top_genes)

    # https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.neighbors.html
    if argv.n_neighbors > 0:
        print(f"scanpy.pp.neighbors on dataset ({str(argv.n_neighbors)} neighbors).")
        scanpy.pp.neighbors(preprocessed_adata, n_neighbors=argv.n_neighbors)

    out_filepath = ""
    subset = 1
    while(not out_filepath):
        pfp = Path(fp)
        _out_filepath = str(pfp.with_name(pfp.stem + f"_preprocessed_{subset}.h5ad"))
        if not os.path.exists(_out_filepath):
            out_filepath = _out_filepath

        subset += 1

    preprocessed_adata.write_h5ad(out_filepath)

if __name__ == "__main__":

    parser = ArgumentParser(
        prog='Preprocess scPrint2 data',
        description="Preprocessing h5ad file to make it compatible with scPrint2."
        )
    
    parser.add_argument("in_file")
    parser.add_argument("--max_cells", type=int, default=-1)
    parser.add_argument("--organism", type=str, default="NCBITaxon:9606") # Human
    parser.add_argument("--cell_type", type=str, default="")
    parser.add_argument("--n_neighbors", type=int, default=20)
    parser.add_argument("--n_top_genes", type=int, default=6000)
    parser.add_argument("--omit_log_transform", action="store_true")

    main(parser.parse_args())
    