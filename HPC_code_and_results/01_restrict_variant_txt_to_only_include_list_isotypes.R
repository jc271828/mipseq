library(tidyverse)
setwd("/data/baughlab/jingxian/mipseq_jc")

variants <- read_tsv("WI.20210121.soft-filter.isotype.txt")
list <- read_tsv("underground.gartersnake_list.txt")$strains
print("Finished loading files")

variants_final <- cbind(variants[, 1:9], variants[, which(colnames(variants) %in% list)])
print("Finished filtering to only include list isotypes")

rm(list = c("variants", "list")) # important for ensuring you have enough memory

write_tsv(variants_final, "WI.20210121.soft-filter.isotype.only.list.isotypes.txt", col_names = TRUE)
print("Finished writing the TXT file which only includes list isotypes")
