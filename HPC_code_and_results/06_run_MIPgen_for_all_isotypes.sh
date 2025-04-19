#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=1 --mem-per-cpu=500M --cpus-per-task=48
#SBATCH --mail-type=END --mail-user=jc750@duke.edu

eval "$(conda shell.bash hook)"
conda activate /data/baughlab/jingxian/conda/

cd /data/baughlab/jingxian/mipseq_jc/192IsotypeFiles/to_process

for f in ./*.2021.5bpflank.bed; do
project_name=$(echo "$f" | awk -F'[./]' '{print $3}');
../../MIPgen/mipgen -regions_to_scan $f -project_name $project_name -min_capture_size 100 -max_capture_size 100 -tag_sizes 0,10 -bwa_genome_index ../../WS276_indices/c_elegans.PRJNA13758.WS276.genomic.fa -snp_file ../../WI.20210121.soft-filter.isotype.jc.list-isotype-unique.StrainNames_hyp_div_masked.vcf.gz &&
rm *all_*;
rm *collapsed_mips*;
rm *coverage*;
rm *_sequences*;
rm *local_*;
rm *oligo_*;
rm *snp_mips*;
rm *progress*
done
