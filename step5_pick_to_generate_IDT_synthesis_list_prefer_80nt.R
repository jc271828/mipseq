library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)

setwd("filtered_MIPs/")

################################################################################################################
easy_strains <- read_xlsx("mips_hard_filter_167EasyStrains_20210813.xlsx")
easy_strains$logistic_score %<>% as.double()
easy_strains$mip_length %>% unique()
easy_strains %<>% mutate(ID = paste(strain, chr, sep = "_"))
easy_strains_IDs <- easy_strains$ID %>% unique()

picked_easy_strains <- data.frame()

for (i in 1:length(easy_strains_IDs)) {
  temp <- filter(easy_strains, ID == easy_strains_IDs[i])
  max_score_temp <- temp$logistic_score %>% max()
  picked_temp <- filter(temp, logistic_score == max_score_temp)
  picked_easy_strains %<>% rbind(picked_temp)
}

summary_picked_easy_strains <-
  picked_easy_strains %>%
  group_by(strain) %>%
  summarise(num_picked_mips = n()) %>%
  arrange(desc(num_picked_mips), strain)

picked_easy_strains_4MIPs_4CHRs <-
  picked_easy_strains %>%
  filter(strain %in% filter(summary_picked_easy_strains, num_picked_mips >= 4)$strain)

picked_easy_strains_4MIPs_4CHRs_new <- data.frame()

for(i in 1:length(unique(picked_easy_strains_4MIPs_4CHRs$strain))) {
  strain_temp = unique(picked_easy_strains_4MIPs_4CHRs$strain)[i]
  temp <- filter(picked_easy_strains_4MIPs_4CHRs, strain == strain_temp)
  max_score_temp = temp$logistic_score %>% sort(decreasing = TRUE) %>% .[1:4]
  picked_temp <- filter(temp, logistic_score %in% max_score_temp)
  picked_easy_strains_4MIPs_4CHRs_new %<>% rbind(picked_temp)
}

picked_easy_strains_4MIPs_4CHRs <- picked_easy_strains_4MIPs_4CHRs_new
rm(list = "picked_easy_strains_4MIPs_4CHRs_new")

write_xlsx(picked_easy_strains_4MIPs_4CHRs, "picked_easy_strains_4MIPs_4CHRs_20210907.xlsx", col_names = TRUE)
################################################################################################################
easy_strains %<>% filter(!(strain %in% unique(picked_easy_strains_4MIPs_4CHRs$strain)) & !(mip_key %in% picked_easy_strains$mip_key))
summary_picked_easy_strains %<>% filter(num_picked_mips < 4)
summary_picked_easy_strains %<>% mutate(still_need = 4 - num_picked_mips)
picked_easy_strains %<>% filter(!(strain %in% picked_easy_strains_4MIPs_4CHRs$strain))

for (i in 1:length(unique(easy_strains$strain))) {
  strain_temp <- unique(easy_strains$strain)[i]
  temp <- filter(easy_strains, strain == strain_temp)
  still_need_temp <- summary_picked_easy_strains %>% filter(strain == strain_temp) %>% .$still_need
  max_score_temp <- temp$logistic_score %>% sort(decreasing = TRUE) %>% .[1:still_need_temp]
  picked_temp <- filter(temp, logistic_score %in% max_score_temp)
  picked_easy_strains %<>% rbind(picked_temp)
}

picked_easy_strains %<>% arrange(strain, chr)
hh <- picked_easy_strains %>% group_by(strain) %>% summarise(count = n()) %>% arrange(count)
# QG4018 has two MIPs with the same score
hh <- picked_easy_strains %>% filter(strain == "QG4018")
# trash the QG4018 MIPs on Chr V with a score of 0.980352
picked_easy_strains %<>% filter(!(strain == "QG4018" & logistic_score == "0.980352" & chr == "V"))

write_xlsx(picked_easy_strains, "picked_easy_strains_4MIPs_some_on_same_CHRs_20210907.xlsx", col_names = TRUE)
################################################################################################################
rm(list = ls())
diff_strains <- read_xlsx("mips_basic_filter_25DifficultStrains_20210813.xlsx", col_names = TRUE)
diff_strains$logistic_score %<>% as.double()
diff_strains$mip_length <-str_length(diff_strains$mip_sequence)
diff_strains$mip_length %>% unique()
diff_strains %<>% mutate(ID = paste(strain, chr, sep = "_"))
diff_strains_IDs <- diff_strains$ID %>% unique()

diff_strains$logistic_score %>% range()
hh <- diff_strains %>% arrange(logistic_score)
# delete the very weird MIP with a logistic score of 6.16872e-05
diff_strains %<>% filter(logistic_score != 6.16872e-05)
diff_strains$logistic_score %>% range()

# pick as many 80nt as possible
picked_diff_strains <- filter(diff_strains, mip_length == 80)

for (i in 1:length(diff_strains_IDs)) {
  temp <- filter(diff_strains, ID == diff_strains_IDs[i])
  max_score_temp <- temp$logistic_score %>% max()
  picked_temp <- filter(temp, logistic_score == max_score_temp)
  picked_diff_strains %<>% rbind(picked_temp)
}

# remove duplicated rows (preserved by wanting more 80nt and also by picking max score per ID)
picked_diff_strains %<>% unique()

summary_picked_diff_strains <-
  picked_diff_strains %>%
  group_by(strain) %>%
  summarise(num_picked_mips = n()) %>%
  arrange(desc(num_picked_mips), strain)

# Except N2 and CX11314, all difficult strains have >= 4 MIPs
# Pick 4 for each while keeping all 80nt's
diff_strains %<>% filter(!(strain %in% c("CX11314", "N2")))
rm(list = c("picked_diff_strains", "summary_picked_diff_strains"))
diff_strains_IDs <- unique(diff_strains$ID)

at80nt_diff_strains <- filter(diff_strains, mip_length == 80)

summary_at80nt_diff_strains <-
  at80nt_diff_strains %>%
  group_by(strain) %>%
  summarise(num_80nt_mips = n()) %>%
  arrange(desc(num_80nt_mips), strain)

summary_at80nt_diff_strains[21:23, "strain"] <- setdiff(diff_strains$strain, summary_at80nt_diff_strains$strain)
summary_at80nt_diff_strains[21:23, "num_80nt_mips"] <- 0
summary_at80nt_diff_strains$still_need <- 4 - summary_at80nt_diff_strains$num_80nt_mips

diff_strains_IDs %<>% setdiff(at80nt_diff_strains$ID)

not80nt_diff_strains <- data.frame()

for (i in 1:length(diff_strains_IDs)) {
  temp <- filter(diff_strains, ID == diff_strains_IDs[i])
  max_score_temp <- temp$logistic_score %>% max()
  picked_temp <- filter(temp, logistic_score == max_score_temp)
  not80nt_diff_strains %<>% rbind(picked_temp)
}

summary_not80nt_diff_strains <-
  not80nt_diff_strains %>%
  group_by(strain) %>%
  summarise(num_not80nt_mips = n()) %>%
  arrange(desc(num_not80nt_mips), strain)

picked_diff_strains <- at80nt_diff_strains

for (i in 1:length(unique(not80nt_diff_strains$strain))) {
  strain_temp <- unique(not80nt_diff_strains$strain)[i]
  temp <- filter(not80nt_diff_strains, strain == strain_temp)
  still_need_temp <- summary_at80nt_diff_strains %>% filter(strain == strain_temp) %>% .$still_need
  max_score_temp <- temp$logistic_score %>% sort(decreasing = TRUE) %>% .[1:still_need_temp]
  picked_temp <- filter(temp, logistic_score %in% max_score_temp)
  picked_diff_strains %<>% rbind(picked_temp)
}

picked_diff_strains <-picked_diff_strains[, c(26, 1:25, 27:28)]
write_xlsx(picked_diff_strains, "picked_diff_strains_20210907.xlsx", col_names = TRUE)
################################################################################################################
rm(list = ls())
tricky_strains <- read_xlsx("mips_relaxed_filter_CX11314_N2_20210813.xlsx", col_names = TRUE)
tricky_strains$logistic_score %<>% as.double()

# no 80nt MIPs for tricky strains. Pick the highest scores.
tricky_strains$mip_length %>% unique()

tricky_strains %<>% mutate(ID = paste(strain, chr, sep = "_"))
# N2 has only three MIPs, so do the picking only for CX11314
CX11314_IDs <- filter(tricky_strains, !(strain == "N2"))$ID %>% unique()

picked_tricky_strains <- filter(tricky_strains, strain == "N2")

for (i in 1:length(CX11314_IDs)) {
  chr_temp <- CX11314_IDs[i] %>% str_split("_") %>% unlist() %>% .[2]
  temp <- filter(tricky_strains, strain == "CX11314" & chr == chr_temp)
  max_score_temp <- temp$logistic_score %>% max()
  picked_temp <- filter(temp, logistic_score == max_score_temp)
  picked_tricky_strains %<>% rbind(picked_temp)
}

min_score_CX11314 <- picked_tricky_strains %>% filter(strain == "CX11314") %>% .$logistic_score %>% sort() %>% .[1]

picked_tricky_strains %<>% filter(!(strain == "CX11314" & logistic_score == min_score_CX11314))

picked_tricky_strains <- picked_tricky_strains[, c(21, 1:20, 22:28)]
write_xlsx(picked_tricky_strains, "picked_tricky_strains_20210907.xlsx", col_names = TRUE)
################################################################################################################
# compile the results
rm(list = ls())

picked_easy_all_diff_CHRs <- read_xlsx("picked_easy_strains_4MIPs_4CHRs_20210907.xlsx", col_names = TRUE)
picked_easy_some_same_CHRs <- read_xlsx("picked_easy_strains_4MIPs_some_on_same_CHRs_20210907.xlsx", col_names = TRUE)

picked_diff_strains <- read_xlsx("picked_diff_strains_20210907.xlsx", col_names = TRUE)

picked_tricky_strains <- read_xlsx("picked_tricky_strains_20210907.xlsx", col_names = TRUE)

all_picked_MIPs <-
  rbind(picked_easy_all_diff_CHRs,
        picked_easy_some_same_CHRs,
        picked_diff_strains,
        picked_tricky_strains) %>%
  arrange(strain, chr)

# should have 191*4+3 = 767 MIPs. Numbers match.
nrow(all_picked_MIPs)

write_xlsx(list(all_picked_MIPs = all_picked_MIPs,
                picked_easy_all_diff_CHRs = picked_easy_all_diff_CHRs %>% arrange(strain, chr),
                picked_easy_some_same_CHRs = picked_easy_some_same_CHRs %>% arrange(strain, chr),
                picked_diff_strains = picked_diff_strains %>% arrange(strain, chr),
                picked_tricky_strains = picked_tricky_strains %>% arrange(strain, chr)),
           "picked_MIPs_for_IDT_synthesis_20210907.xlsx",
           col_names = TRUE)

#717 MIPs at 80 nt, to be synthesized in plates
nrow(filter(all_picked_MIPs, mip_length == 80))

# 50 MIPs over 80 nt, to be synthesized in tubes
nrow(filter(all_picked_MIPs, mip_length != 80))

write_xlsx(list(MIPs_at_80nt = filter(all_picked_MIPs, mip_length == 80) %>% arrange(strain, chr),
                MIPs_longer_than_80nt = filter(all_picked_MIPs, mip_length > 80) %>% arrange(strain, chr)),
           "picked_MIPs_for_IDT_synthesis_grouped_by_length_20210907_corrected.xlsx",
           col_names = TRUE)

