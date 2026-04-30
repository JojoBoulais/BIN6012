from finetune import *

from argparse import ArgumentParser
from utils_filepath import *
from utils import load_model_with_cuda_if_avail
import scanpy as sc
import torch
import os


def main(argv):
    # Load the training and validation dataset
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
    adata = sc.read_h5ad(adata_fp)

    # Use only the myoblast cells
    adata = adata[adata.obs["cell_type_ontology_term_id"] == "CL:0000513",]

    # Creating the fintuner
    print("Initializing the fintuner...")
    finetuner = FinetuneGN(ft_mode="full")

    # Finetuning the model
    print("Finetuning the model...")
    finetuned_model = finetuner(model=model, adata=adata)

    # Saving the finetuned model
    torch.save(finetuned_model, "fintuned_GRN_model.pth")


if __name__ == "__main__":
    parser = ArgumentParser(
        prog="Finetuning the scprint-2 model.",
        description="Generating a Denoised h5ad file from a preprocessed dataset."
    )

    parser.add_argument("model_file")
    parser.add_argument("h5ad_file")

    main(parser.parse_args())
