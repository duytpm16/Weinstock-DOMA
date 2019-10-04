### Options and libraries
options(stringsAsFactors = FALSE)
library(qtl2)





### Load viewer data and permutations
load('~/Desktop/weinstock_doma_viewer_v2.Rdata')
perms <- readRDS('~/Desktop/doma_otu_perm_1000_v2.rds')







### Permutation summaries
perms <- summary_scan1perm(object = perms, 
                           alpha = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.67))






### Store permutation summaries
write.csv(perms, file = '~/Desktop/doma_otu_additive_perms_1000_summary_v2.csv')
dataset.doma.otu$perms <- list(additive = perms)







### Save
rm(perms)
save.image(file = '~/Desktop/weinstock_doma_viewer_v2.Rdata')

