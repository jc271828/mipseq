#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=1 --mem=500000M --cpus-per-task=56
#SBATCH --mail-type=END --mail-user=jc750@duke.edu

eval "$(conda shell.bash hook)"
conda activate /data/baughlab/jingxian/conda/

Rscript /home/jc750/02_restrict_variants_to_list_strain_private_ones.R

