library(tidyverse)
library(matrixStats)

setwd("/data/baughlab/jingxian/mipseq_jc")

variants_final <- read_tsv("WI.20210121.soft-filter.isotype.only.list.isotypes.txt")
print("Finished loading files")

for (i in 10: ncol(variants_final)) variants_final[, i] <- lapply(variants_final[, i], str_detect, pattern = "1/1") %>% unlist()
print("Finished finding ALT homozygous isotypes for each variant")

temp <- as.matrix(variants_final)
print("Finished getting the matrix to be counted for the number of ALT HMZ isotypes per variant")

variants_final$num_of_ALT_HMZ_strains_per_SNV <- rowCounts(temp, value = TRUE, na.rm = TRUE)
print("Finished counting the number of ALT HMZ isotypes per variant")

rm(list = "temp") # important for ensuring you have enough memory

variants_final <- filter(variants_final, num_of_ALT_HMZ_strains_per_SNV == 1)
print("Finished filtering to only include variants present in one of the list isotypes")

variants_final$ALT_HMZ_strain <- NA
variants_final$ALT_HMZ_strain <- as.character(variants_final$ALT_HMZ_strain)

for (i in 1: nrow(variants_final)) variants_final[i, "ALT_HMZ_strain"] <- colnames(variants_final)[which(variants_final[i, 1:(ncol(variants_final)-2)] == TRUE)]
print("Finished matching isotypes to isotype-unique variants")

write_tsv(variants_final, "WI.20210121.soft-filter.isotype.jc.list-isotype-unique.txt", col_names = TRUE)
print("Finished writing the TXT file which only includes list isotypes and their unique variants")

