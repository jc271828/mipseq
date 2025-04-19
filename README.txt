# of variants in the original VCF: 6295334
# of list-unique variants (Igor filtered): 221099
# of list-unique variants (Jingxian filtered): 1043837
# of variants after masking strain-private hyper divergent region + 20kb on each end in Jingxian filtered VCF: 679401

106 of Igor's unique MIPs are in my raw output
99 of Igor's unique MIPs are in my filtered output

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

# hard filter:
# # basic filter:
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

# relaxed filter:
# 738 MIPs of tricky strains loaded
# 721 MIPs of tricky strains after restricting arms to not span SNPs
# 155 MIPs of tricky strains after restricting 0 < distance from end of UMI <= 40
# 122 MIPs of tricky strains after restricting REF and ALT1 to have different first bases
# 14 MIPs after removing 0/1, 1/0, ./1, 1/., ./0, and 0/.
# filters up until this step are applied within the cluster
# 12 MIPs of tricky strains after filtering out REF over any ALTs and ALT1 over any other ALTs
# 12 MIPs of tricky strains after filtering out ALTs with REF's first base present in list isotypes