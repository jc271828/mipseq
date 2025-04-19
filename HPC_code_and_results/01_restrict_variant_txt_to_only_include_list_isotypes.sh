#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=1 –mem-per-cpu=9000M --cpus-per-task=56
#SBATCH --mail-type=END --mail-user=jc750@duke.edu
eval "$(conda shell.bash hook)"
conda activate /data/baughlab/jingxian/conda/
cd ~/scripts

Rscript ./01_restrict_variant_txt_to_only_include_list_isotypes.R
