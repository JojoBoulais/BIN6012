from argparse import ArgumentParser
import os

from scprint2.tasks import Denoiser
import scanpy

from utils_filepath import *
from utils import *
from pathlib import Path

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

    # Forcing human organism
    model.organisms = adata.var["organism"].unique()

    denoise = Denoiser(
        batch_size=40 if adata.X.sum(1).mean() < 50_000 else 20,
        max_len=8_000,
        max_cells=100_0000,
        doplot=False,
        num_workers=8,
        predict_depth_mult=5,
        downsample_expr=0.7,
    )

    # This add a new layer to anndata: pred_adata.layers["scprint_mu"]
    _, idx, nadata = denoise(model, adata)

    out_filepath = get_filepath("data", "")

    out_filepath = out_filepath + Path(adata_fp).stem + "_denoised.h5ad"
    nadata.write_h5ad(out_filepath)

    
if __name__ == "__main__":

    parser = ArgumentParser(
        prog='scPrint2 Denoising on a preprocessed dataset',
        description="Generating a Denoised h5ad file from a preprocessed dataset."
        )
    
    parser.add_argument("model_file")
    parser.add_argument("h5ad_file")

    main(parser.parse_args())
