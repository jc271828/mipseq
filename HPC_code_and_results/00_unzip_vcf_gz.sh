#!/bin/bash
#SBATCH --mail-type=END --mail-user=jc750@duke.edu
#SBATCH --nodes=1 --ntasks-per-node=1 --mem=500000M --cpus-per-task=56

cd /data/baughlab/jingxian/mipseq_jc

gunzip -c WI.20210121.soft-filter.isotype.vcf.gz > WI.20210121.soft-filter.isotype.vcf
