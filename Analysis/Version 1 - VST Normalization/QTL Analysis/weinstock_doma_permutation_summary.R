### Options and libraries
options(stringsAsFactors = FALSE)
library(qtl2)







### Load viewer
load('~/Desktop/Weinstock_DOMA/Viewer/Version 1 - VST and RankZ/weinstock_doma_viewer_v1.Rdata')








### Load in additive permutations
otu_perm    <- readRDS('~/Desktop/Weinstock_DOMA/Viewer/Version 1 - VST and RankZ/Permutations 1000/doma_otu_additive_perm_1000_v1.rds')
genus_perm  <- readRDS('~/Desktop/Weinstock_DOMA/Viewer/Version 1 - VST and RankZ/Permutations 1000/doma_genus_additive_perm_1000_v1.rds')
family_perm <- readRDS('~/Desktop/Weinstock_DOMA/Viewer/Version 1 - VST and RankZ/Permutations 1000/doma_family_additive_perm_1000_v1.rds')
order_perm  <- readRDS('~/Desktop/Weinstock_DOMA/Viewer/Version 1 - VST and RankZ/Permutations 1000/doma_order_additive_perm_1000_v1.rds')
class_perm  <- readRDS('~/Desktop/Weinstock_DOMA/Viewer/Version 1 - VST and RankZ/Permutations 1000/doma_class_additive_perm_1000_v1.rds')
phylum_perm <- readRDS('~/Desktop/Weinstock_DOMA/Viewer/Version 1 - VST and RankZ/Permutations 1000/doma_phylum_additive_perm_1000_v1.rds')









### Permutation summaries
otu_summary    <- summary_scan1perm(object = otu_perm,    alpha = c(0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.67))
genus_summary  <- summary_scan1perm(object = genus_perm,  alpha = c(0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.67))
family_summary <- summary_scan1perm(object = family_perm, alpha = c(0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.67))
order_summary  <- summary_scan1perm(object = order_perm,  alpha = c(0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.67))
class_summary  <- summary_scan1perm(object = class_perm,  alpha = c(0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.67))
phylum_summary <- summary_scan1perm(object = phylum_perm, alpha = c(0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.67))









### Save additive permutations
dataset.doma.otu$perms    <- list(additive = otu_summary)
dataset.doma.genus$perms  <- list(additive = genus_summary)
dataset.doma.family$perms <- list(additive = family_summary)
dataset.doma.order$perms  <- list(additive = order_summary)
dataset.doma.class$perms  <- list(additive = class_summary)
dataset.doma.phylum$perms <- list(additive = phylum_summary)









### Save
rm(list = ls()[!grepl('dataset[.]|genoprobs|K|map|markers', ls())])
save.image(file = '~/Desktop/Weinstock_DOMA/Viewer/Version 1 - VST and RankZ/weinstock_doma_viewer_v1.Rdata')




