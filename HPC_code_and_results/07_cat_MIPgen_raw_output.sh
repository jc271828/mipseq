#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=1 --mem-per-cpu=500M --cpus-per-task=2
#SBATCH --mail-type=END --mail-user=jc750@duke.edu

eval "$(conda shell.bash hook)"
conda activate /data/baughlab/jingxian/conda/

cd /data/baughlab/jingxian/mipseq_jc/192IsotypeFiles/processed

cat *picked* > MIPgen_output_raw_20210813.txt

Rscript ~/scripts/07_cat_and_tidy_MIPgen_raw_output.R
