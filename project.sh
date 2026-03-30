#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=slurm_output/augment_data.%j.out
#SBATCH --error=slurm_output/augment_data.%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --account=def-gsarah

module load python/3.12
source ~/.venv/bin/activate

#python3 prepare_ontologies.py
#python3 augment_data.py
#scprint2 gninfer --adata /home/jordboul/scratch/BIN6012/scPrint-2/data/myocardial_infarction_subset.h5ad --ckpt_path /home/jordboul/scratch/BIN6012/scPrint-2/model_ckpt/small-v2.ckpt --output_filename grn.h5ad
scprint2 gninfer --adata /home/jordboul/scratch/BIN6012/scPrint-2/data/myocardial_infarction_subset.h5ad --ckpt_path /home/jordboul/scratch/BIN6012/scPrint-2/model_ckpt/small-v2.ckpt --species 'NCBITaxon:9606' --cell_type 'cardiac muscle myoblast' --output_filename grn.h5ad
# lamin init --storage ./testdb --name testdb --modules bionty
# scprint2 easy_setup