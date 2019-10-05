### Options and libraries
options(stringsAsFactors = FALSE)
library(qtl2)





### Load viewer data and permutations
load('weinstock_doma_viewer_v1.Rdata')





### Extract data
otu   <- dataset.doma.otu$data$rz
covar <- dataset.doma.otu$covar.matrix







### Estimate heritability per chromosome for each OTU
chr_herit <- as.data.frame(matrix(0, nrow = 20, ncol = ncol(otu),
                                  dimnames = list(paste0('chr',c(1:19, 'X')), colnames(otu))))

for(i in colnames(chr_herit)){
    for(j in rownames(chr_herit)){
        chr = gsub('chr', '', j)
        chr_herit[j, i] <- est_herit(pheno    = otu[, i, drop = FALSE], 
                                     kinship  = K[[chr]],
                                     addcovar = covar,
                                     cores    = 0)[1]
    }
    print(i)
}





### Estimate overall heritability for each OTU
K_overall <- calc_kinship(probs = genoprobs, type = 'overall', cores = 0)

overall_herit <- as.data.frame(matrix(0, nrow = 1, ncol = ncol(otu),
                                      dimnames = list('overall', colnames(otu))))

for(i in colnames(overall_herit)){
    overall_herit[1, i] <- est_herit(pheno = otu[, i, drop = FALSE], 
                                     kinship  = K_overall,
                                     addcovar = covar,
                                     cores    = 0)[1]
    print(i)
}









### Save
dataset.doma.otu$herit <- list(overall = overall_herit, chromosome = chr_herit)


rm(list = ls()[!grepl('dataset[.]|genoprobs|map|markers|K', ls())])
rm(K_overall)
save.image('weinstock_doma_viewer_v1.Rdata')
