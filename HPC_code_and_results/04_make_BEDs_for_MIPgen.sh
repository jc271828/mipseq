#!/bin/bash
#SBATCH --nodes=4 --ntasks-per-node=1 --mem=10000M --cpus-per-task=16 --mail-type=END --mail-user=jc750@duke.edu

eval "$(conda shell.bash hook)"
conda activate /data/baughlab/jingxian/conda/

cd /data/baughlab/jingxian/mipseq_jc
perl /home/jc750/scripts/04_make_BEDs_for_MIPgen.pl WI.20210121.soft-filter.isotype.jc.list-isotype-unique.StrainNames_hyp_div_masked.txt 5 250 WBcel235_genomic_chrom_length.txt '.2021.5bpflank.bed'
