library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)

processed <- read_xlsx("mips_basic_filter_all_strains_20210813.xlsx", col_names = TRUE)
nrow(processed) %>% paste("MIPs that have passed the basic filter are loaded", sep = " ") %>% print()

diff <- read_tsv("25difficult_strains_20210813.txt", col_names = FALSE)$X1
length(diff) %>% paste("strains have < 4 good MIPs", sep = " ") %>% print()

processed_diff <- filter(processed, strain %in% diff)
nrow(processed_diff) %>% paste("MIPs of difficult strains passed the basic filter", sep = " ") %>% print()

mips_per_strain <- processed_diff$strain %>% table() %>% as.data.frame()
colnames(mips_per_strain) <- c("strain", "occ")
mips_per_strain %<>% arrange(desc(occ))
filter(mips_per_strain, occ >= 4) %>% nrow() %>% paste("strains have >= 4 good MIPs") %>% print()
filter(mips_per_strain, occ == 3) %>% nrow() %>% paste("strains have 3 good MIPs") %>% print()
filter(mips_per_strain, occ == 2) %>% nrow() %>% paste("strains have 2 good MIPs") %>% print()
filter(mips_per_strain, occ == 1) %>% nrow() %>% paste("strains have 1 good MIP") %>% print()
no_good_mips <- setdiff(diff, mips_per_strain$strain)
length(no_good_mips) %>% paste("strains do not have good MIPs") %>% print()
paste("They are", paste(no_good_mips, collapse = " "), sep = " ") %>% print()

# save files
write_xlsx(processed_diff, "mips_basic_filter_25DifficultStrains_20210813.xlsx", col_names = TRUE)
write_xlsx(mips_per_strain, "mips_per_strain_basic_filter_25DifficultStrains_20210813.xlsx", col_names = TRUE)


# limit MIP lengths to 80 bases
processed_diff$mip_length <- str_length(processed_diff$mip_sequence)
processed_diff %<>% filter(mip_length == 80)
nrow(processed_diff) %>% paste("MIPs of difficult strains after restricting length == 80", sep = " ") %>% print()

mips_per_strain <- processed_diff$strain %>% table() %>% as.data.frame()
colnames(mips_per_strain) <- c("strain", "occ")
mips_per_strain %<>% arrange(desc(occ))
filter(mips_per_strain, occ >= 4) %>% nrow() %>% paste("strains have >= 4 good MIPs") %>% print()
filter(mips_per_strain, occ == 3) %>% nrow() %>% paste("strains have 3 good MIPs") %>% print()
filter(mips_per_strain, occ == 2) %>% nrow() %>% paste("strains have 2 good MIPs") %>% print()
filter(mips_per_strain, occ == 1) %>% nrow() %>% paste("strains have 1 good MIP") %>% print()
no_good_mips <- setdiff(diff, mips_per_strain$strain)
length(no_good_mips) %>% paste("strains do not have good MIPs") %>% print()
paste("They are", paste(no_good_mips, collapse = " "), sep = " ") %>% print()

# save files
write_xlsx(processed_diff, "mips_basic_filter_plus_80bp_limit_20DifficultStrains_20210813.xlsx", col_names = TRUE)
write_xlsx(mips_per_strain, "mips_per_strain_basic_filter_plus_80bp_limit_20DifficultStrains_20210813.xlsx", col_names = TRUE)

