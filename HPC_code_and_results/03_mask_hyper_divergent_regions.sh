#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=1 –mem-per-cpu=1000M --cpus-per-task=4
#SBATCH --mail-type=END --mail-user=jc750@duke.edu
eval "$(conda shell.bash hook)"
conda activate /data/baughlab/jingxian/conda/
cd /data/baughlab/jingxian/mipseq_jc

Rscript ~/scripts/03_mask_hyper_divergent_regions.R
