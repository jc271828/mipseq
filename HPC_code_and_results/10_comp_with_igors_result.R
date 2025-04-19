library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)

setwd("/data/baughlab/jingxian/mipseq_jc/192IsotypeFiles/processed/results")

easy <- read_xlsx("mips_hard_filter_167EasyStrains_20210813.xlsx")
diff <- read_xlsx("mips_basic_filter_25DifficultStrains_20210813.xlsx")
tricky <- read_xlsx("mips_relaxed_filter_CX11314_N2_20210813.xlsx")

filtered_seqs <- union(easy$mip_sequence, diff$mip_sequence) %>% union(tricky$mip_sequence)
raw_seqs <- read_xlsx("MIPgen_output_raw_20210813.xlsx") %>% .$mip_sequence %>% unique()
igor_seqs <- read_xlsx("../../../igor_selected_v2_20180130.xlsx") %>% .$mip_sequence %>% unique()

intersect(igor_seqs, raw_seqs) %>% length() %>% paste("of Igor's unique MIPs are in my raw output", sep = " ") %>% print()
intersect(igor_seqs, filtered_seqs) %>% length() %>% paste("of Igor's unique MIPs are in my filtered outputput", sep = " ") %>% print()



