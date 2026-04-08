#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=slurm_output/grn_plot.%j.out
#SBATCH --error=slurm_output/grn_plot.%j.err
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

module load python/3.12

# Activer le virtual environnement qui contient scprint2
source ~/scprint2/bin/activate

# ============= CUSTOM SCRIPTS =============

#python3 /home/jordboul/scratch/scPrint2/BIN6012/code/preprocess.py /home/jordboul/scratch/scPrint2/BIN6012/data/myocardial_infarction.h5ad --max_cells -1 #--cell_type "cardiac muscle myoblast"
#python3 /home/jordboul/scratch/scPrint2/BIN6012/code/grn.py /home/jordboul/scratch/scPrint2/BIN6012/model/small-v2.ckpt /home/jordboul/scratch/scPrint2/BIN6012/data/myocardial_infarction_preprocessed_1.h5ad
python3 /home/marco30/BIN6012/code/grn_plot.py "/home/marco30/BIN6012/grn_fromscprint.h5ad" --pathway Wnt/β-catenin

# ============= scPRINT2 scripts =============
#scprint2 gninfer --adata /home/jordboul/scratch/BIN6012/scPrint-2/data/myocardial_infarction_subset.h5ad --ckpt_path /home/jordboul/scratch/BIN6012/scPrint-2/model_ckpt/small-v2.ckpt --output_filename grn.h5ad
#scprint2 gninfer --adata /home/jordboul/scratch/BIN6012/scPrint-2/data/myocardial_infarction_subset.h5ad --ckpt_path /home/jordboul/scratch/BIN6012/scPrint-2/model_ckpt/small-v2.ckpt --species 'NCBITaxon:9606' --cell_type 'cardiac muscle myoblast' --output_filename /home/jordboul/scratch/BIN6012/scPrint-2/grn.h5ad
#lamin init --storage ./testdb --name testdb --modules bionty
#scprint2 easy_setup
#python3 prepare_ontologies.py
