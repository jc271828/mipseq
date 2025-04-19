library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)

processed <- read_xlsx("mips_basic_filter_all_strains_20210813.xlsx", col_names = TRUE)
nrow(processed) %>% paste("MIPs that have passed the basic filter are loaded", sep = " ") %>% print()

# limit logistic scores to >= 0.98
processed %<>% filter(logistic_score >= 0.98)
nrow(processed) %>% paste("MIPs after restricting score >= 0.98", sep = " ") %>% print()

# limit MIP lengths to 80 bases
processed$mip_length <- str_length(processed$mip_sequence)
processed %<>% filter(mip_length == 80)
nrow(processed) %>% paste("MIPs after restricting length == 80", sep = " ") %>% print()

# hard filter:
# # essential filter:
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
# # additional restrictions:
# 49277 MIPs after restricting score >= 0.98
# 7130 MIPs after restricting length == 80

mips_per_strain <- processed$strain %>% table() %>% as.data.frame()
colnames(mips_per_strain) <- c("strain", "occ")
mips_per_strain %<>% arrange(desc(occ))
filter(mips_per_strain, occ >= 4) %>% nrow() %>% paste("strains have >= 4 good MIPs") %>% print()
filter(mips_per_strain, occ == 3) %>% nrow() %>% paste("strains have 3 good MIPs") %>% print()
filter(mips_per_strain, occ == 2) %>% nrow() %>% paste("strains have 2 good MIPs") %>% print()
filter(mips_per_strain, occ == 1) %>% nrow() %>% paste("strains have 1 good MIP") %>% print()
no_good_mips <- read_tsv("underground.gartersnake_list.txt", col_names = TRUE)$strains %>% setdiff(mips_per_strain$strain)
length(no_good_mips) %>% paste("strains do not have good MIPs") %>% print()
paste("They are", paste(no_good_mips, collapse = " "), sep = " ") %>% print()

# save files
write_xlsx(mips_per_strain, "mips_per_strain_hard_filter_all_strains_20210813.xlsx", col_names = TRUE)

easy_strains <-
  processed %>%
  filter(strain %in% filter(mips_per_strain, occ >= 4)$strain) %>%
  arrange(strain, chr, var_pos)
easy_strains <- easy_strains[, c(26, 1:25, 27)]
write_xlsx(easy_strains, "mips_hard_filter_167EasyStrains_20210813.xlsx", col_names = TRUE)

c(as.character(filter(mips_per_strain, occ < 4)$strain), no_good_mips) %>%
  as.data.frame() %>%
  write_tsv("25difficult_strains_20210813.txt", col_names = FALSE)

