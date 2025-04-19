library(tidyverse)
library(magrittr)

setwd("/data/baughlab/jingxian/mipseq_jc")

jc_strain_names <- read_tsv("jc_strain_names.txt", col_names = TRUE)
lee2021 <- read_tsv("lee2021.txt", col_names = TRUE)

for (i in 1: nrow(jc_strain_names)) {
  pos <- jc_strain_names$POS[i]
  lee2021_filtered <- filter(lee2021, strain_chr == jc_strain_names$strain_chr[i])
  lee2021_filtered %<>% mutate(contain = (start_minus_20kb <= pos & end_plus_20kb >= pos))
  jc_strain_names[i, "to_remove"] <- any(lee2021_filtered$contain)
}

write_tsv(jc_strain_names, "jc_strain_names_with_to_remove_column.txt", col_names = TRUE)

jc_strain_names %>% filter(to_remove == FALSE) %>% .[, 1:5] %>% write_tsv("jc_strain_names_hyp_div_masked.txt", col_names = FALSE)

