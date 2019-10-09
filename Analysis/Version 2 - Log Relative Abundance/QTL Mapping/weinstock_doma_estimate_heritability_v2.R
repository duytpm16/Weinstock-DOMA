### Options and libraries
options(stringsAsFactors = FALSE)
library(qtl2)





### Load viewer data and permutations
load("~/Desktop/Weinstock_DOMA/Viewer/Version 2 - Centered Log Ratio/weinstock_doma_viewer_v2.Rdata")
K_overall <- calc_kinship(probs = genoprobs, type = 'overall', cores = 0)





### Estimate heritability begins
for(d in grep('dataset[.]', ls(), value = TRUE)){
  
  ### Extract data
  ds <- get(d)
  expr  <- ds$data$norm
  covar <- ds$covar.matrix
  
  
  ### Estimate heritability per chromosome
  chr_herit <- as.data.frame(matrix(0, nrow = 20, ncol = ncol(expr),
                                    dimnames = list(paste0('chr',c(1:19, 'X')), colnames(expr))))
  
  for(i in colnames(chr_herit)){
    for(j in rownames(chr_herit)){
      chr = gsub('chr', '', j)
      chr_herit[j, i] <- est_herit(pheno    = expr[, i, drop = FALSE], 
                                   kinship  = K[[chr]],
                                   addcovar = covar,
                                   cores    = 0)[1]
    }
    print(i)
  }
  
  
  
  
  
  ### Estimate overall heritability 
  overall_herit <- as.data.frame(matrix(0, nrow = 1, ncol = ncol(expr),
                                        dimnames = list('overall', colnames(expr))))
  
  for(i in colnames(overall_herit)){
    overall_herit[1, i] <- est_herit(pheno    = expr[, i, drop = FALSE], 
                                     kinship  = K_overall,
                                     addcovar = covar,
                                     cores    = 0)[1]
    print(i)
  }
  
  
  
  
  ### Store to dataset
  ds$herit <- list(overall = overall_herit, chromosome = chr_herit)
  assign(d, ds)
}







### Save
rm(list = ls()[!grepl('dataset[.]|genoprobs|map|markers|K', ls())])
rm(K_overall)
save.image('~/Desktop/Weinstock_DOMA/Viewer/Version 2 - Centered Log Ratio/weinstock_doma_viewer_v2.Rdata')
