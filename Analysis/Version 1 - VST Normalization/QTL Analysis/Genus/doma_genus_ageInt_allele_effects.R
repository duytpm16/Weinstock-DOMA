### Options and libraries
options(stringsAsFactors = FALSE)
library(qtl2)
library(dplyr)




### Load data
load('~/Desktop/Weinstock_DOMA/Viewer/Version 1 - VST and RankZ/weinstock_doma_viewer_v1.Rdata')






### Extract data
genus <- dataset.doma.genus$data$rz
covar <- dataset.doma.genus$covar.matrix
peaks <- dataset.doma.genus$lod.peaks$age_int
samps <- dataset.doma.genus$annot.samples





### Stratify by sex
stopifnot(samps$mouse.id == rownames(genus))
m6_genus  <- genus[samps$mouse.id[samps$cohort.age == 6], ]
m12_genus <- genus[samps$mouse.id[samps$cohort.age == 12], ]
m18_genus <- genus[samps$mouse.id[samps$cohort.age == 18], ]


covar <- covar %>% as.data.frame() %>% select(-cohort.age) %>% as.matrix()









### Get allele effects for each month seperately
peaks[,paste0(LETTERS[1:8], '_m6')]  <- 0
peaks[,paste0(LETTERS[1:8], '_m12')] <- 0
peaks[,paste0(LETTERS[1:8], '_m18')] <- 0

for(i in 1:nrow(peaks)){
  
    gp <- genoprobs[,peaks$qtl.chr[i]]
    gp[[1]] <- gp[[1]][,,peaks$marker.id[i], drop = FALSE]
  
  
    m6_blup  <- scan1blup(genoprobs = gp,
                          pheno     = m6_genus[,peaks$data.name[i], drop = FALSE],
                          kinship   = K[[peaks$qtl.chr[i]]],
                          addcovar  = covar)
  
    m12_blup <- scan1blup(genoprobs = gp,
                          pheno     = m12_genus[,peaks$data.name[i], drop = FALSE],
                          kinship   = K[[peaks$qtl.chr[i]]],
                          addcovar  = covar)
    
    m18_blup <- scan1blup(genoprobs = gp,
                          pheno     = m18_genus[,peaks$data.name[i], drop = FALSE],
                          kinship   = K[[peaks$qtl.chr[i]]],
                          addcovar  = covar)
  
  
    peaks[i,paste0(LETTERS[1:8], '_m6')]  <- m6_blup[1,LETTERS[1:8]]
    peaks[i,paste0(LETTERS[1:8], '_m12')] <- m12_blup[1,LETTERS[1:8]]
    peaks[i,paste0(LETTERS[1:8], '_m18')] <- m18_blup[1,LETTERS[1:8]]  
}




dataset.doma.genus$lod.peaks$age_int <- peaks







### Save
rm(list = ls()[!grepl('dataset[.]|genoprobs|map|markers|K', ls())])
save.image('~/Desktop/Weinstock_DOMA/Viewer/Version 1 - VST and RankZ/weinstock_doma_viewer_v1.Rdata')
