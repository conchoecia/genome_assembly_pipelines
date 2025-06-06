#!/bin/bash

#SBATCH --job-name=FCSgx   # This is the name of the parent job
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=512GB
#SBATCH --time=0-2:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# Thes are example parameters for running the GAP NCBI FCS adapter tool
#   with SLURM, including the expected RAM consumption and time.

# For this script we are running the NCBI FCS adapter tool
CORES=24

SNAKEFILE=/your/path/to/genome_assembly_pipelines/snakefiles/GAP_NCBI_fcs_adaptor_gx_trim_singularity

source ~/.bashrc

module load conda # something specific to your HPC environment
conda activate /lisc/user/schultz/miniconda3/envs/snakemake8 # Activate the conda environment with snakemake

snakemake --cores ${CORES} -p --rerun-incomplete --snakefile ${SNAKEFILE}
