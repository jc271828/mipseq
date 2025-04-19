library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)
library(matrixStats)

# keep filtering
# filter by:
# 1) no REF over any ALT present in our list
# 2) no ALT1 over any other ALT present in our list
# 3) if all other ALTs present in our list and REF differ in their first base => acceptable
# 4) if ALT X and REF have the same first base, but none of our list strains have homozygous ALT X in the genotype => acceptable
raw <- read_xlsx("MIPgen_output_preliminary_filter_20210813.xlsx", col_names = TRUE)
has_unord_alleles <- read_xlsx("has_non_0hmz_and_1hmz_and_NAhmz_alleles_20210813.xlsx", col_names = TRUE)

# to trash 1) REF heterozygosity and 2) ALT1 heterozygosity
has_unord_alleles$identifier <- paste(has_unord_alleles$strain, has_unord_alleles$chr, has_unord_alleles$var_pos, sep = "_")
raw$identifier <- paste(raw$strain, raw$chr, raw$var_pos, sep = "_")
ref_alt1_het_row_num <- has_unord_alleles$allele %>% lapply(str_detect, pattern = "0/|/0|1/|/1") %>% unlist() %>% which()
raw %<>% filter(!(identifier %in% has_unord_alleles$identifier[ref_alt1_het_row_num]))
nrow(raw) %>% paste("MIPs after filtering out REF over any ALTs and ALT1 over any other ALTs", sep = " ") %>% print()

# to filter based on 3) and 4)
# manually filter out rows where ALT X and REF have the same first base, and homozygous ALT X is present in our list strains
has_unord_alleles <- has_unord_alleles[-c(ref_alt1_het_row_num), ]
has_unord_alleles %<>% add_column(ref_other_alts_same_first_base = NA)
has_unord_alleles$ref_other_alts_same_first_base %<>% as.logical()
for (i in 1: nrow(has_unord_alleles)) {
  alt_first_bases <- has_unord_alleles$ALT[i] %>% str_split(pattern = ",") %>% unlist() %>% lapply(substr, start = 1, stop = 1) %>% unlist()
  has_unord_alleles[i, "ref_other_alts_same_first_base"] <- any(alt_first_bases == substr(has_unord_alleles$REF[i], start = 1, stop = 1))
}
has_unord_alleles %<>% filter(ref_other_alts_same_first_base == TRUE) %>% select(-c("ref_other_alts_same_first_base"))

has_unord_alleles %<>% add_column(alts_with_ref_first_base_present_in_list_isotypes = NA)
has_unord_alleles$alts_with_ref_first_base_present_in_list_isotypes %<>% as.logical()
for (i in 1: nrow(has_unord_alleles)) {
  alt_first_bases <- has_unord_alleles$ALT[i] %>% str_split(pattern = ",") %>% unlist() %>% lapply(substr, start = 1, stop = 1) %>% unlist()
  which_alt <- which(alt_first_bases == substr(has_unord_alleles$REF[i], start = 1, stop = 1))
  has_unord_alleles[i, "alts_with_ref_first_base_present_in_list_isotypes"] <- str_detect(has_unord_alleles$allele[i], pattern = paste(which_alt, collapse = "|"), negate = FALSE)
}
has_unord_alleles %<>% filter(alts_with_ref_first_base_present_in_list_isotypes == TRUE) %>% select(-c("alts_with_ref_first_base_present_in_list_isotypes"))
nrow(has_unord_alleles) %>% paste("MIPs to trash because they have ALTs with REF's first base present in list isotypes", sep = " ") %>% print()
raw %<>% filter(!(identifier %in% has_unord_alleles$identifier)) %>% select(-c("identifier"))
nrow(raw) %>% paste("MIPs after filtering out ALTs with REF's first base present in list isotypes", sep = " ") %>% print()

# basic filter:
# 533969 MIPs loaded from raw MIPgen output
# 519110 MIPs after restricting arms to not span SNPs
# 291601 MIPs after restricting score > 0.97
# 283478 MIPs after restricting failure flag==000
# 75556 MIPs after restricting 0 < distance from end of UMI <= 40
# 65428 MIPs after restricting REF and ALT1 to have different first bases
# 62998 MIPs after removing 0/1, 1/0, ./1, 1/., ./0, and 0/.
# filters up until this step are applied within the cluster
# 62566 MIPs after filtering out REF over any ALTs and ALT1 over any other ALTs
# 62518 MIPs after filtering out ALTs with REF's first base present in list isotypes

write_xlsx(raw, "mips_basic_filter_all_strains_20210813.xlsx", col_names = TRUE)


