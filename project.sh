#!/bin/bash
#SBATCH --job-name=umap
#SBATCH --output=slurm_output/umap.%j.out
#SBATCH --error=slurm_output/umap.%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --account=def-gsarah

module load python/3.12

# Activer le virtual environnement qui contient scprint2
source ~/.venv/bin/activate

# ============= CUSTOM SCRIPTS =============

#python3 /home/jordboul/scratch/scPrint2/BIN6012/code/preprocess.py /home/jordboul/scratch/scPrint2/BIN6012/data/myocardial_infarction.h5ad --max_cells -1 #--cell_type "cardiac muscle myoblast"

python /home/jordboul/scratch/scPrint2/BIN6012/code/do_umap_pca.py /home/jordboul/scratch/scPrint2/BIN6012/data/myocardial_infarction_preprocessed_2.h5ad

#python3 /home/jordboul/scratch/scPrint2/BIN6012/code/grn.py /home/jordboul/scratch/scPrint2/BIN6012/model/small-v2.ckpt /home/jordboul/scratch/scPrint2/BIN6012/data/myocardial_infarction_preprocessed_2.h5ad
#python3 /home/jordboul/scratch/scPrint2/BIN6012/code/grn_plot.py --normal_h5ad "/home/jordboul/scratch/scPrint2/BIN6012/grn/myocardial_infarction_preprocessed_2_small-v2_cardiac_muscle_myoblast_normal_test.h5ad" --disease_h5ad "/home/jordboul/scratch/scPrint2/BIN6012/grn/myocardial_infarction_preprocessed_2_small-v2_cardiac_muscle_myoblast_myocardial_infarction_test.h5ad" --pathway all

#python3 /home/jordboul/scratch/scPrint2/BIN6012/code/grn_n_plot.py /home/jordboul/scratch/scPrint2/BIN6012/model/small-v2.ckpt /home/jordboul/scratch/scPrint2/BIN6012/data/myocardial_infarction_preprocessed_2.h5ad --cell_type "fibroblast of cardiac tissue"

# ============= scPRINT2 scripts =============
#scprint2 gninfer --adata /home/jordboul/scratch/BIN6012/scPrint-2/data/myocardial_infarction_subset.h5ad --ckpt_path /home/jordboul/scratch/BIN6012/scPrint-2/model_ckpt/small-v2.ckpt --output_filename grn.h5ad
#scprint2 gninfer --adata /home/jordboul/scratch/BIN6012/scPrint-2/data/myocardial_infarction_subset.h5ad --ckpt_path /home/jordboul/scratch/BIN6012/scPrint-2/model_ckpt/small-v2.ckpt --species 'NCBITaxon:9606' --cell_type 'cardiac muscle myoblast' --output_filename /home/jordboul/scratch/BIN6012/scPrint-2/grn.h5ad
#lamin init --storage ./testdb --name testdb --modules bionty
#scprint2 easy_setup
#python3 prepare_ontologies.py
