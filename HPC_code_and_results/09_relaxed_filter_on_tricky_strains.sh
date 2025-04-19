#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=1 --mem-per-cpu=1000M --cpus-per-task=16
#SBATCH --mail-type=END --mail-user=jc750@duke.edu

eval "$(conda shell.bash hook)"
conda activate /data/baughlab/jingxian/conda/

cd /data/baughlab/jingxian/mipseq_jc/192IsotypeFiles/processed

Rscript ~/scripts/09_relaxed_filter_on_tricky_strains.R
