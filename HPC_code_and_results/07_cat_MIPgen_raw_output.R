library(tidyverse)
library(magrittr)
library(writexl)

setwd("/data/baughlab/jingxian/mipseq_jc/192IsotypeFiles/processed")

raw <- read_tsv("MIPgen_output_raw_20210813.txt", col_names = TRUE)
colnames(raw)[1] <- "mip_key"

raw <- raw[str_detect(raw$logistic_score, "logistic_score", negate = TRUE), ]

raw$strain <- c(NA) %>% as.character()
for (i in 1: nrow(raw)) {
  raw[i, "strain"] <- raw$mip_name[i] %>%
    strsplit(split = "_", fixed = FALSE) %>%
    unlist() %>%
    .[1]
}

raw$logistic_score %<>% as.double()
write_xlsx(raw, "MIPgen_output_raw_20210813.xlsx", col_names = TRUE)
