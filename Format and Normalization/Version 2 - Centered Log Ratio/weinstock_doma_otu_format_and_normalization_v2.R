### Options and libraries
options(stringsAsFactors = FALSE)
library(qtl2convert)
library(tidyverse)
library(qtl2)








### Load data
load('~/Desktop/Weinstock_DOMA/Genotypes/DODB/qtl2/JAC_megaMUGA_genoprobs_qtl2.RData')
samples <- read.csv('~/Desktop/Weinstock_DOMA/Phenotypes/do_mice/DO_CrossSectional_Population.csv')
chrY_M  <- read.csv('~/Desktop/Weinstock_DOMA/Phenotypes/do_mice/JAC_crosssectional_sex_chrM_Y_20180618.csv')
otu  <- readRDS('~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/Modified/0.5/otu_raw_count_0.5_cleaned.rds')
taxa <- readRDS('~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/Modified/0.5/otu_taxa_table_0.5_cleaned.rds')









### Removing OTUs with <= 5% prevalence according to Hoan
#     otu: 403 x 384
otu  <- otu[, colSums(otu > 0) > (nrow(otu) * 0.05)]
taxa <- taxa %>% 
          filter(OTU %in% colnames(otu))







### Create sample annotations
samples <- samples %>%
              filter(Cohort == 'Cross-Sectional') %>%
              filter(!grepl('HDO-', Mouse.ID)) %>%
              select(Mouse.ID, Sex, DOB, Generation, Cohort.Age, Coat.Color, Wean.Date, DOD, Age.Death, From.Cage, JCMS.Cage, Cage, Ear.Notch, Cohort) %>%
              mutate(Mouse.ID = gsub('-', '.', Mouse.ID)) %>%
              left_join(y = chrY_M %>% select(X, chrM, chrY) %>% dplyr::rename(Mouse.ID = X), by = 'Mouse.ID') %>%
              filter(Mouse.ID %in% rownames(otu)) %>%
              mutate(Sex = factor(Sex), Generation = factor(Generation)) %>%
              `colnames<-`(tolower(colnames(.))) %>%
              distinct()








### Normalize OTU using Hoan / Dong-binh's code
codaSeq.clr <- function(x, samples.by.row=TRUE){
  
  if(min(x) == 0) stop("0 values must be replaced, estimated, or eliminated")
  if(samples.by.row == TRUE){margin=1}
  if(samples.by.row == FALSE){margin=2}

  return ( t(apply(x, margin, function(x){log(x) - mean(log(x))})) )
}
norm <- apply(otu, 2, function(x) x / sum(x))
norm <- norm + 1
norm <- codaSeq.clr(x = norm, samples.by.row = TRUE)











### Covariates
covar <- model.matrix(~ sex + generation + cohort.age, data = samples)[, -1, drop = FALSE]
covar <- covar[,c('sexM', 'cohort.age', grep('generation', colnames(covar), value = TRUE))]
colnames(covar)[1] <- 'sex'
rownames(covar) <- samples$mouse.id



covar.info <- data.frame(sample.column = c('sex', 'cohort.age', 'generation'),
                         covar.column  = c('sex', 'cohort.age', 'generation'),
                         display.name  = c('Sex', 'Cohort Age', 'Generation'),
                         interactive   = c(TRUE, TRUE, FALSE),
                         primary       = c(TRUE, TRUE, FALSE),
                         lod.peaks     = c('sex_int', 'age_int', NA))












### Annotate columns
annot.phenotype <- data.frame(data.name   = c(colnames(samples), taxa$OTU),
                              short.name  = c(colnames(samples), taxa$OTU),
                              R.name      = c(colnames(samples), taxa$OTU),
                              description = c('Mouse identifier', 'Gender of mouse: F (female) or M (male)', 'Mouse date-of-birth', 'Generation of mouse', 'Age of mouse in months',
                                              'Color coat of mouse', 'Date of wean', 'Date of death', 'Age at death in days', 'Cage of mouse', 'JCMS cage',
                                              'Cage of mouse', 'Ear notch', 'JAC Cohort', 'chrM', 'chrY', apply(taxa, 1, paste0, collapse = ' - ')),
                              units       = NA,
                              category    = c(rep('Demographic', ncol(samples)), rep('Microbiome', nrow(taxa))),
                              R.category  = c(rep('Demographic', ncol(samples)), rep('Microbiome', nrow(taxa))),
                              is.id       = c(TRUE, rep(FALSE, ncol(samples) + nrow(taxa) - 1)),
                              is.numeric  = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
                                              rep(TRUE, nrow(taxa))),
                              is.date     = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, rep(FALSE, ncol(samples) + nrow(taxa) - 8)),
                              is.factor   = c(FALSE, TRUE, FALSE, TRUE, rep(FALSE, ncol(samples) + nrow(taxa) - 4)),
                              factor.levels = c(NA, 'F:M', NA, '8:9:10:11:12', rep(FALSE, ncol(samples) + nrow(taxa) - 4)),
                              is.covar    = c(FALSE, TRUE, FALSE, TRUE, TRUE, rep(FALSE, ncol(samples) + nrow(taxa) - 5)),
                              is.pheno    = c(rep(FALSE, ncol(samples)), rep(TRUE, nrow(taxa))),
                              is.derived  = FALSE,
                              omit        = FALSE,
                              use.covar   = c(rep(NA, ncol(samples)), rep('sex:cohort.age:generation', nrow(taxa))),
                              transformation = c(rep(NA, ncol(samples)), rep('rankz', nrow(taxa))))













### QTL viewer format
dataset.doma.otu <- list(annot.phenotype = as_tibble(annot.phenotype),
                         annot.samples   = as_tibble(samples),
                         covar.matrix    = as.matrix(covar),
                         covar.info      = as_tibble(covar.info),
                         data            = list(norm = as.matrix(norm),
                                                raw  = as.matrix(otu)),
                         datatype        = 'phenotype',
                         display.name    = 'DOMA OTU Abundance',
                         lod.peaks       = list(),
                         taxa            = as_tibble(taxa))







### Markers as tibble
markers <- as_tibble(markers)








### Reducing genoprobs for lower memory. Then compute kinship
genoprobs <- probs_qtl2_to_doqtl(probs = genoprobs)
genoprobs <- genoprobs[dimnames(genoprobs)[[1]] %in% samples$mouse.id,,]
genoprobs <- probs_doqtl_to_qtl2(probs = genoprobs, map = as.data.frame(markers), marker_column = 'marker.id', pos_column = 'pos')

K <- calc_kinship(probs = genoprobs, type = 'loco', cores = 0)









### Save
rm(list = ls()[!grepl('dataset[.]|K|map|markers|genoprobs', ls())])
save.image(file = '~/Desktop/Weinstock_DOMA/Viewer/Version 2 - Centered Log Ratio/weinstock_doma_viewer_v2.Rdata')
