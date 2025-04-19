library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)
library(matrixStats)

setwd("/data/baughlab/jingxian/mipseq_jc/192IsotypeFiles/processed")

raw <- read_xlsx("MIPgen_output_raw_20210813.xlsx", col_names = TRUE)
paste(nrow(raw), "MIPs loaded", sep = " ") %>% print()

raw$var_pos <- c(NA) %>% as.double()
for (i in 1: nrow(raw)) {
  raw[i, "var_pos"] <- raw$mip_name[i] %>%
    strsplit(split = "[:]") %>%
    unlist() %>%
    .[2] %>%
    as.double()
}
raw %<>% mutate(var_in_arms = (lig_probe_start <= var_pos & lig_probe_stop >= var_pos) | (ext_probe_start <= var_pos & ext_probe_stop >= var_pos))
raw %<>% filter(var_in_arms == FALSE)
paste(nrow(raw), "MIPs after restricting arms to not span SNPs", sep = " ") %>% print()

raw %<>% select(-c("var_in_arms"))
raw %<>% filter(logistic_score > 0.97)
paste(nrow(raw), "MIPs after restricting score > 0.97", sep = " ") %>% print()

raw %<>% filter(failure_flags == "000")
paste(nrow(raw), "MIPs after restricting failure flag==000", sep = " ") %>% print()

raw$dis_from_UMI_end <- c(NA) %>% as.double()
raw$lig_probe_start %<>% as.double()
raw$lig_probe_stop %<>% as.double()
for (i in 1: nrow(raw)) {
  if (raw[i, "probe_strand"] == "-") {
    raw[i, "dis_from_UMI_end"] <- raw[i, "var_pos"] - raw[i, "lig_probe_start"] + 1
  }else {
    raw[i, "dis_from_UMI_end"] <- raw[i, "lig_probe_stop"] - raw[i, "var_pos"] + 1
  }
}
raw %<>% filter((dis_from_UMI_end > 0) & (dis_from_UMI_end <= 40))
paste(nrow(raw), "MIPs after restricting 0 < distance from end of UMI <= 40", sep = " ") %>% print()

raw$REF <- c(NA) %>% as.character()
raw$ALT <- c(NA) %>% as.character()
for (i in 1: nrow(raw)) {
  raw[i, "REF"] <- raw$mip_name[i] %>%
    strsplit(split = ":|->|_", fixed = FALSE) %>%
    unlist() %>%
    .[4]
  
  raw[i, "ALT"] <- raw$mip_name[i] %>%
    strsplit(split = ":|->|_", fixed = FALSE) %>%
    unlist() %>%
    .[5]
}

raw$ALT_num <- c(NA) %>% as.double()
for (i in 1: nrow(raw)) {
  raw[i, "ALT_num"] <- raw$ALT[i] %>%
    strsplit(split = ",") %>%
    unlist() %>%
    length()
}

raw$REF_and_ALT1_same_first_base <- c(NA) %>% as.logical()
for (i in 1: nrow(raw)) {
  if (substr(raw$REF[i], start = 1, stop = 1) == substr(raw$ALT[i], start = 1, stop = 1)) {
    raw[i, "REF_and_ALT1_same_first_base"] <- TRUE
  } else {
    raw[i, "REF_and_ALT1_same_first_base"] <- FALSE
  }
}

raw %<>% filter(REF_and_ALT1_same_first_base == FALSE) %>% select(-c("REF_and_ALT1_same_first_base"))
nrow(raw) %>% paste("MIPs after restricting REF and ALT1 to have different first bases", sep = " ") %>% print()

jc_list_isotypes <- read_tsv("../../WI.20210121.soft-filter.isotype.only.list.isotypes.txt", col_names = TRUE) %>% .[, c(1:2, 4:5, 10:ncol(.))]
colnames(jc_list_isotypes)[1] <- "chr"
colnames(jc_list_isotypes)[2] <- "var_pos"

merged <-
  raw[, c("chr", "var_pos", "strain")] %>%
  merge(jc_list_isotypes, by = c("chr", "var_pos"), all.x = TRUE, all.y = FALSE)

rm(list = "jc_list_isotypes")

temp1 <- merged
for (i in 6: ncol(temp1)) {
  temp1[, i] <- lapply(temp1[, i], str_detect, pattern = "0/1|1/0|\\./0|0/\\.|\\./1|1/\\.", negate = FALSE) %>% unlist()
}
temp1 %<>% add_column(how_many_0_or_1_hets = NA)
temp1$how_many_0_or_1_hets %<>% as.numeric()
temp1$how_many_0_or_1_hets <- rowCounts(as.matrix(temp1[, 6:(ncol(temp1) - 1)]), value = TRUE, na.rm = TRUE)

temp1 %<>% arrange(chr, var_pos)
raw %<>% arrange(chr, var_pos)
merged %<>% arrange(chr, var_pos)
raw$how_many_0_or_1_hets <- temp1$how_many_0_or_1_hets
merged$how_many_0_or_1_hets <- temp1$how_many_0_or_1_hets

rm("temp1")

raw %<>% filter(how_many_0_or_1_hets == 0)
nrow(raw) %>% paste("MIPs after removing 0/1, 1/0, ./1, 1/., ./0, and 0/.", sep = " ") %>% print()

merged %<>% filter(how_many_0_or_1_hets == 0)
raw %<>% select(-c("how_many_0_or_1_hets"))
merged %<>% select(-c("how_many_0_or_1_hets"))

raw %<>% arrange(chr, var_pos)
merged %<>% arrange(chr, var_pos)
merged$ALT_num <- raw$ALT_num
merged %<>% filter(ALT_num >= 2)
merged <- merged[, c(1:5, ncol(merged), 6: (ncol(merged) - 1))]


# detect "unordinary" genotypes (neither 1/1 nor 0/0 nor ./.)
temp2 <- merged
for (i in 7: ncol(temp2)) {
  temp2[, i] <- lapply(temp2[, i], str_detect, pattern = "1/1|0/0|\\./\\.", negate = TRUE) %>% unlist()
}
temp2$num_of_other_alleles <- rowCounts(as.matrix(temp2[, 7: ncol(temp2)]), value = TRUE, na.rm = TRUE)

temp2 %<>% arrange(chr, var_pos)
merged %<>% arrange(chr, var_pos)
merged$num_of_other_alleles <- temp2$num_of_other_alleles

rm(list = "temp2")

merged %<>% filter(num_of_other_alleles >= 2)
merged <- merged[, c(1:6, ncol(merged), 7: (ncol(merged) - 1))]

merged$which_strain <- c(NA) %>% as.character()
merged$allele <- c(NA) %>% as.character()

fun <- function(x) {
  str_split(x, pattern = ":") %>% unlist() %>% .[1]
}
for(i in 1: nrow(merged)) {
  col_num_minus7 <- lapply(merged[i, 8: ncol(merged)], str_detect, pattern = "1/1|0/0|\\./\\.", negate = TRUE) %>% unlist() %>% which()
  merged[i, "which_strain"] <-
    colnames(merged)[col_num_minus7 + 7] %>%
    paste(collapse = ",")
  merged[i, "allele"] <-
    merged[i, (col_num_minus7 + 7)] %>%
    lapply(fun) %>%
    unlist() %>%
    paste(collapse = ",")
}
merged <- merged[, c(1:7, (ncol(merged) - 1): ncol(merged), 8: (ncol(merged) - 2))]
nrow(merged) %>% paste("MIPs, out of the MIPs that meet previous criteria, have alleles that are not 1/1 or 0/0 or ./.", sep = " ") %>% print()

write_xlsx(raw, "MIPgen_output_preliminary_filter_20210813.xlsx", col_names = TRUE)
write_xlsx(merged, "has_non_0hmz_and_1hmz_and_NAhmz_alleles_20210813.xlsx", col_names = TRUE)
