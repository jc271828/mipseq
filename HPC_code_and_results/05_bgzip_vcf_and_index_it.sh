#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=1 --mem=10000M --cpus-per-task=16 --mail-type=END --mail-user=jc750@duke.edu

eval "$(conda shell.bash hook)"
conda activate /data/baughlab/jingxian/conda/

cd /data/baughlab/jingxian/mipseq_jc

bgzip -c WI.20210121.soft-filter.isotype.jc.list-isotype-unique.StrainNames_hyp_div_masked.vcf > WI.20210121.soft-filter.isotype.jc.list-isotype-unique.StrainNames_hyp_div_masked.vcf.gz
tabix -p vcf WI.20210121.soft-filter.isotype.jc.list-isotype-unique.StrainNames_hyp_div_masked.vcf.gz
