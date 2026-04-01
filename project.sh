#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=slurm_output/grn.%j.out
#SBATCH --error=slurm_output/grn.%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --account=def-gsarah

module load python/3.12

# Activer le virtual environnement qui contient scprint2
source ~/.venv/bin/activate

# ============= CUSTOM SCRIPTS =============
#python3 /home/jordboul/scratch/BIN6012/scPrint-2/code/preprocess.py myocardial_infarction.h5ad --max_cells 10000 #--cell_type "cardiac muscle myoblast"
python3 /home/jordboul/scratch/BIN6012/scPrint-2/code/grn.py small-v2.ckpt myocardial_infarction_preprocessed_4.h5ad

# ============= scPRINT2 scripts =============
#scprint2 gninfer --adata /home/jordboul/scratch/BIN6012/scPrint-2/data/myocardial_infarction_subset.h5ad --ckpt_path /home/jordboul/scratch/BIN6012/scPrint-2/model_ckpt/small-v2.ckpt --output_filename grn.h5ad
#scprint2 gninfer --adata /home/jordboul/scratch/BIN6012/scPrint-2/data/myocardial_infarction_subset.h5ad --ckpt_path /home/jordboul/scratch/BIN6012/scPrint-2/model_ckpt/small-v2.ckpt --species 'NCBITaxon:9606' --cell_type 'cardiac muscle myoblast' --output_filename /home/jordboul/scratch/BIN6012/scPrint-2/grn.h5ad
# lamin init --storage ./testdb --name testdb --modules bionty
# scprint2 easy_setup
#python3 prepare_ontologies.py