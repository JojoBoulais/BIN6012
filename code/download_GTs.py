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
