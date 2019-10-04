### Options and libraries
options(stringsAsFactors = FALSE)
library(qtl2)





### Load data
load('~/Desktop/weinstock_doma_viewer_v2.Rdata')
 





### Extract data
otu   <- dataset.doma.otu$data$norm
covar <- dataset.doma.otu$covar.matrix
peaks <- dataset.doma.otu$lod.peaks$sex_int
samps <- dataset.doma.otu$annot.samples





### Stratify by sex
stopifnot(samps$mouse.id == rownames(otu))
f_otu <- otu[samps$mouse.id[samps$sex == 'F'], ]
m_otu <- otu[samps$mouse.id[samps$sex == 'M'], ]

covar <- covar %>% as.data.frame() %>% select(-sex) %>% as.matrix()









### Get allele effects for males and females seperately
peaks[,paste0(LETTERS[1:8], '_female')] <- 0
peaks[,paste0(LETTERS[1:8], '_male')] <- 0

for(i in 1:nrow(peaks)){
  
    gp <- genoprobs[,peaks$qtl.chr[i]]
    gp[[1]] <- gp[[1]][,,peaks$marker.id[i], drop = FALSE]
    
    
    f_blup <- scan1blup(genoprobs = gp,
                        pheno     = f_otu[,peaks$data.name[i], drop = FALSE],
                        kinship   = K[[peaks$qtl.chr[i]]],
                        addcovar  = covar)
    
    m_blup <- scan1blup(genoprobs = gp,
                        pheno     = m_otu[,peaks$data.name[i], drop = FALSE],
                        kinship   = K[[peaks$qtl.chr[i]]],
                        addcovar  = covar)
    
    
    peaks[i,paste0(LETTERS[1:8], '_female')] <- f_blup[1,LETTERS[1:8]]
    peaks[i,paste0(LETTERS[1:8], '_male')]   <- m_blup[1,LETTERS[1:8]]
    
}




dataset.doma.otu$lod.peaks$sex_int <- peaks





### Save
rm(list = ls()[!grepl('dataset[.]|genoprobs|map|markers|K', ls())])
save.image('~/Desktop/weinstock_doma_viewer_v2.Rdata')
