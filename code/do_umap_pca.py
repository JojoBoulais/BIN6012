from argparse import ArgumentParser
from utils_filepath import *

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import os
import umap
import matplotlib.pyplot as plt
import scanpy

def main(argv):

    if os.path.exists(argv.h5ad_file):
        adata_fp = argv.h5ad_file
    else:
        adata_fp = get_filepath("data", argv.h5ad_file)
        if not os.path.exists(adata_fp):
            print(f"{adata_fp} does not exist.")
            exit(1)

    print("loading dataset.")
    adata = scanpy.read_h5ad(adata_fp)

    results_fp = get_filepath("results")
    

    cell_types = ["cardiac muscle myoblast", "fibroblast of cardiac tissue"]

    methods = ["pca", "umap"]

    for method in methods:

        out_dir = os.path.join(results_fp, method)
        os.makedirs(out_dir, exist_ok=True)

        for cell_type in cell_types:

            adata_subtype = adata[adata.obs["cell_type"] == cell_type].copy()

            adata_subtype = adata_subtype[::5].copy()

            X_input = adata_subtype.X.toarray() if hasattr(adata_subtype.X, "toarray") else adata_subtype.X

            if method == "umap":
                X = umap.UMAP().fit_transform(X_input)
            else:
                scaler = StandardScaler()
                scaled_data = scaler.fit_transform(X_input)

                # Keep enough components to explain 95% of the variance
                pca = PCA(n_components=2)
                X = pca.fit_transform(scaled_data)

            fig, ax = plt.subplots(figsize=(12, 6))

            for disease in ["normal", "myocardial infarction"]:
                adata_disease = adata_subtype.obs["disease"] == disease

                ax.scatter(
                    X[adata_disease, 0],
                    X[adata_disease, 1],
                    alpha=0.6,
                    label=disease
                )

            ax.legend(
                loc="center left",
                bbox_to_anchor=(1.03, 0.5),
                fancybox=True,
                frameon=False,
                handletextpad=1.3
            )

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            ax.set_title(f"{method.upper()} utilisant les paramètres par défaut")

            plt.tight_layout()

            cell_type_str = cell_type.replace(" ", "_")
            plt.savefig(os.path.join(out_dir, cell_type_str + ".png"))
            plt.close(fig)

if __name__ == "__main__":

    parser = ArgumentParser(
        prog='UMAP on a preprocessed dataset',
        description="UMAP on a preprocessed dataset"
    )

    parser.add_argument("h5ad_file")

    main(parser.parse_args())