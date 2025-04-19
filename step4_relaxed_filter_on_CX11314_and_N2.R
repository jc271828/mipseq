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
tricky <- read_xlsx("tricky_strains_preliminary_filter_20210813.xlsx", col_names = TRUE)
has_unord_alleles <- read_xlsx("tricky_strains_has_non_0hmz_and_1hmz_and_NAhmz_alleles_20210813.xlsx", col_names = TRUE)

# to trash 1) REF heterozygosity and 2) ALT1 heterozygosity
has_unord_alleles$identifier <- paste(has_unord_alleles$strain, has_unord_alleles$chr, has_unord_alleles$var_pos, sep = "_")
tricky$identifier <- paste(tricky$strain, tricky$chr, tricky$var_pos, sep = "_")
ref_alt1_het_row_num <- has_unord_alleles$allele %>% lapply(str_detect, pattern = "0/|/0|1/|/1") %>% unlist() %>% which()
tricky %<>% filter(!(identifier %in% has_unord_alleles$identifier[ref_alt1_het_row_num]))
nrow(tricky) %>% paste("MIPs of tricky strains after filtering out REF over any ALTs and ALT1 over any other ALTs", sep = " ") %>% print()

# Dataframe has_unord_alleles does not have remaining rows after removing REF and ALT1 heterozygous rows
# # to filter based on 3) and 4)
# # manually filter out rows where ALT X and REF have the same first base, and homozygous ALT X is present in our list strains
# has_unord_alleles <- has_unord_alleles[-c(ref_alt1_het_row_num), ]
# has_unord_alleles %<>% add_column(ref_other_alts_same_first_base = NA)
# has_unord_alleles$ref_other_alts_same_first_base %<>% as.logical()
# for (i in 1: nrow(has_unord_alleles)) {
#   alt_first_bases <- has_unord_alleles$ALT[i] %>% str_split(pattern = ",") %>% unlist() %>% lapply(substr, start = 1, stop = 1) %>% unlist()
#   has_unord_alleles[i, "ref_other_alts_same_first_base"] <- any(alt_first_bases == substr(has_unord_alleles$REF[i], start = 1, stop = 1))
# }
# has_unord_alleles %<>% filter(ref_other_alts_same_first_base == TRUE) %>% select(-c("ref_other_alts_same_first_base"))

# has_unord_alleles %<>% add_column(alts_with_ref_first_base_present_in_list_isotypes = NA)
# has_unord_alleles$alts_with_ref_first_base_present_in_list_isotypes %<>% as.logical()
# for (i in 1: nrow(has_unord_alleles)) {
#   alt_first_bases <- has_unord_alleles$ALT[i] %>% str_split(pattern = ",") %>% unlist() %>% lapply(substr, start = 1, stop = 1) %>% unlist()
#   which_alt <- which(alt_first_bases == substr(has_unord_alleles$REF[i], start = 1, stop = 1))
#   has_unord_alleles[i, "alts_with_ref_first_base_present_in_list_isotypes"] <- str_detect(has_unord_alleles$allele[i], pattern = paste(which_alt, collapse = "|"), negate = FALSE)
# }
# has_unord_alleles %<>% filter(alts_with_ref_first_base_present_in_list_isotypes == TRUE) %>% select(-c("alts_with_ref_first_base_present_in_list_isotypes"))
# nrow(has_unord_alleles) %>% paste("MIPs of tricky strains to trash because they have ALTs with REF's first base present in list isotypes", sep = " ") %>% print()
# tricky %<>% filter(!(identifier %in% has_unord_alleles$identifier)) %>% select(-c("identifier"))
# nrow(tricky) %>% paste("MIPs of tricky strains after filtering out ALTs with REF's first base present in list isotypes", sep = " ") %>% print()

# relaxed filter:
# 738 MIPs of tricky strains loaded
# 721 MIPs of tricky strains after restricting arms to not span SNPs
# 155 MIPs of tricky strains after restricting 0 < distance from end of UMI <= 40
# 122 MIPs of tricky strains after restricting REF and ALT1 to have different first bases
# 14 MIPs after removing 0/1, 1/0, ./1, 1/., ./0, and 0/.
# filters up until this step are applied within the cluster
# 12 MIPs of tricky strains after filtering out REF over any ALTs and ALT1 over any other ALTs
# 12 MIPs of tricky strains after filtering out ALTs with REF's first base present in list isotypes

mips_per_strain <- tricky$strain %>% table() %>% as.data.frame()
colnames(mips_per_strain) <- c("strain", "occ")
mips_per_strain %<>% arrange(desc(occ))

# save files
tricky %<>% mutate(mip_length = str_length(mip_sequence)) %>% select(-c("identifier"))
write_xlsx(tricky, "mips_relaxed_filter_CX11314_N2_20210813.xlsx", col_names = TRUE)

