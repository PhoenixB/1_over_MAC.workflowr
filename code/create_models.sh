#!/bin/zsh

#SBATCH --mail-user=p.jerney@students.unibe.ch
#SBATCH --mail-type=end,fail,array_tasks

#SBATCH --job-name="Run models with R"
#SBATCH --output=../output/create_models.%A.out

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --time=24:00:00

#SBATCH --partition=amd

#### Environment variables ####

R_MAX_VSIZE=256Gb
export R_MAX_VSIZE

#### Shell commands ####

Rscript --verbose --no-save --no-restore ./create_models.R
exit
