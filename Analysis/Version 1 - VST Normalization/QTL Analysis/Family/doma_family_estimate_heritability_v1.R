### Options and libraries
options(stringsAsFactors = FALSE)
library(qtl2)





### Load viewer data and permutations
load('~/Desktop/Weinstock_DOMA/Viewer/Version 1 - VST and RankZ/weinstock_doma_viewer_v1.Rdata')









### Extract data
family <- dataset.doma.family$data$rz
covar  <- dataset.doma.family$covar.matrix









### Estimate heritability per chromosome for each family
chr_herit <- as.data.frame(matrix(0, nrow = 20, ncol = ncol(family),
                                  dimnames = list(paste0('chr',c(1:19, 'X')), colnames(family))))

for(i in colnames(chr_herit)){
    for(j in rownames(chr_herit)){
        chr = gsub('chr', '', j)
        chr_herit[j, i] <- est_herit(pheno    = family[, i, drop = FALSE], 
                                     kinship  = K[[chr]],
                                     addcovar = covar,
                                     cores    = 0)[1]
    }
    print(i)
}









### Estimate overall heritability for each family
K_overall <- calc_kinship(probs = genoprobs, type = 'overall', cores = 0)

overall_herit <- as.data.frame(matrix(0, nrow = 1, ncol = ncol(family),
                                      dimnames = list('overall', colnames(family))))

for(i in colnames(overall_herit)){
    overall_herit[1, i] <- est_herit(pheno    = family[, i, drop = FALSE], 
                                     kinship  = K_overall,
                                     addcovar = covar,
                                     cores    = 0)[1]
    print(i)
}














### Save
dataset.doma.family$herit <- list(overall   = overall_herit, 
                                 chromosome = chr_herit)


rm(list = ls()[!grepl('dataset[.]|genoprobs|map|markers|K', ls())])
rm(K_overall)
save.image('~/Desktop/Weinstock_DOMA/Viewer/Version 1 - VST and RankZ/weinstock_doma_viewer_v1.Rdata')
