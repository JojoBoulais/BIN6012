from bengrn import get_GT_db, get_perturb_gt, get_sroy_gt, BenGRN
from utils_filepath import *
import pandas as pd

DBs = ["collectri", "omnipath"]

for db in DBs:    
    df_db = get_GT_db(db)

    out_filepath = get_filepath("data", db + ".parquet")
    df_db.to_parquet(out_filepath)


gwps = get_perturb_gt()
gwps.var.index = gwps.var.gene_name

out_filepath = get_filepath("data", "gwps.h5ad")
gwps.write_h5ad(out_filepath)



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


# df_gt_inter.to_parquet(gt_inter_fp + ".parquet")
