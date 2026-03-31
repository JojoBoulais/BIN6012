from argparse import ArgumentParser
import os
from pathlib import Path

from scdataloader import utils, Preprocessor, DataModule
from utils_filepath import *
import scanpy

def main(argv):

    fp = get_data_filepath(argv.in_file)

    if not os.path.exists(fp):
        print(f"{fp} does not exists")
        exit(1)

    adata = scanpy.read_h5ad(fp)

    if argv.max_cells > 0:
        adata = adata[:argv.max_cells, :].copy()

    if argv.cell_type:
        adata = adata[adata.obs["cell_type"] == argv.cell_type, :].copy()

    if argv.organism:
        adata.obs["organism_ontology_term_id"] = argv.organism
        adata = adata.copy()

    preprocessor = Preprocessor(
        do_postp=False,
        force_preprocess=True
    )

    preprocessed_adata = preprocessor(adata.copy())

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

    main(parser.parse_args())
    