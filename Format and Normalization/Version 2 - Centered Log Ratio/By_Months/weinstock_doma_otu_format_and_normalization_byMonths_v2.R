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







### Divide samples into motnhs
weeks <- split(samples, samples$Cohort.Age)









### Get samples by month
otu_m6 <- otu[rownames(otu) %in% weeks[['6']]$JAC.Mouse.ID, ]
otu_m6 <- otu_m6[, colSums(otu_m6 > 0) > (nrow(otu_m6) * 0.05)]


otu_m12 <- otu[rownames(otu) %in% weeks[['12']]$JAC.Mouse.ID, ]
otu_m12 <- otu_m12[, colSums(otu_m12 > 0) > (nrow(otu_m12) * 0.05)]


otu_m18 <- otu[rownames(otu) %in% weeks[['18']]$JAC.Mouse.ID, ]
otu_m18 <- otu_m18[, colSums(otu_m18 > 0) > (nrow(otu_m18) * 0.05)]









### Normalize OTU using Hoans'/Dong-bing's normalization. 
codaSeq.clr <- function(x, samples.by.row=TRUE){
  
  if(min(x) == 0) stop("0 values must be replaced, estimated, or eliminated")
  if(samples.by.row == TRUE){margin=1}
  if(samples.by.row == FALSE){margin=2}
  
  return ( t(apply(x, margin, function(x){log(x) - mean(log(x))})) )
}
norm_m6  <- apply(otu_m6, 2, function(x) x / sum(x))
norm_m6  <- norm_m6 + 1
norm_m6  <- codaSeq.clr(x = norm_m6, samples.by.row = TRUE)

norm_m12 <- apply(otu_m12, 2, function(x) x / sum(x))
norm_m12 <- norm_m12 + 1
norm_m12 <- codaSeq.clr(x = norm_m12, samples.by.row = TRUE)

norm_m18 <- apply(otu_m18, 2, function(x) x / sum(x))
norm_m18 <- norm_m18 + 1
norm_m18 <- codaSeq.clr(x = norm_m18, samples.by.row = TRUE)









### Get taxa for each month
taxa_m6  <- taxa %>% filter(OTU %in% colnames(otu_m6))
taxa_m12 <- taxa %>% filter(OTU %in% colnames(otu_m12))
taxa_m18 <- taxa %>% filter(OTU %in% colnames(otu_m18))







### Create sample annotations
samples <- samples %>%
              filter(Cohort == 'Cross-Sectional') %>%
              filter(!grepl('HDO-', Mouse.ID)) %>%
              select(Mouse.ID, Sex, DOB, Generation, Cohort.Age, Coat.Color, Wean.Date, DOD, Age.Death, From.Cage, JCMS.Cage, Cage, Ear.Notch, Cohort) %>%
              mutate(Mouse.ID = gsub('-', '.', Mouse.ID)) %>%
              left_join(y = chrY_M %>% select(X, chrM, chrY) %>% dplyr::rename(Mouse.ID = X), by = 'Mouse.ID') %>%
              mutate(Sex = factor(Sex), Generation = factor(Generation)) %>%
              `colnames<-`(tolower(colnames(.))) %>%
              filter(mouse.id %in% c(rownames(otu_m6), rownames(otu_m12), rownames(otu_m18))) %>%
              distinct()

samples_m6  <- samples %>% filter(mouse.id %in% rownames(otu_m6))
samples_m12 <- samples %>% filter(mouse.id %in% rownames(otu_m12))
samples_m18 <- samples %>% filter(mouse.id %in% rownames(otu_m18))









### Covariates
covar <- model.matrix(~ sex + generation, data = samples)[, -1, drop = FALSE]
colnames(covar)[1] <- 'sex'
rownames(covar) <- samples$mouse.id

covar_m6  <- covar[rownames(covar) %in% rownames(class_m6),]
covar_m12 <- covar[rownames(covar) %in% rownames(class_m12),]
covar_m18 <- covar[rownames(covar) %in% rownames(class_m18),]



covar.info <- data.frame(sample.column = c('sex', 'generation'),
                         covar.column  = c('sex', 'generation'),
                         display.name  = c('Sex', 'Generation'),
                         interactive   = c(TRUE, FALSE),
                         primary       = c(TRUE, FALSE),
                         lod.peaks     = c('sex_int', NA))











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

annot.pheno.m6  <- annot.phenotype %>% filter(data.name %in% c(colnames(samples), colnames(otu_m6)))
annot.pheno.m12 <- annot.phenotype %>% filter(data.name %in% c(colnames(samples), colnames(otu_m12)))
annot.pheno.m18 <- annot.phenotype %>% filter(data.name %in% c(colnames(samples), colnames(otu_m18)))











### QTL viewer format
dataset.otu.m6 <- list(annot.phenotype = as_tibble(annot.pheno.m6),
                       annot.samples   = as_tibble(samples_m6),
                       covar.matrix    = as.matrix(covar_m6),
                       covar.info      = as_tibble(covar.info),
                       data            = list(norm = as.matrix(norm_m6),
                                              raw  = as.matrix(otu_m6)),
                       datatype        = 'phenotype',
                       display.name    = 'DOMA OTU Abundance Month 6',
                       lod.peaks       = list(),
                       taxa            = as_tibble(taxa_m6))


dataset.otu.m12 <- list(annot.phenotype = as_tibble(annot.pheno.m12),
                        annot.samples   = as_tibble(samples_m12),
                        covar.matrix    = as.matrix(covar_m12),
                        covar.info      = as_tibble(covar.info),
                        data            = list(norm = as.matrix(norm_m12),
                                               raw  = as.matrix(otu_m12)),
                        datatype        = 'phenotype',
                        display.name    = 'DOMA OTU Abundance Month 12',
                        lod.peaks       = list(),
                        taxa            = as_tibble(taxa_m12))


dataset.otu.m18 <- list(annot.phenotype = as_tibble(annot.pheno.m18),
                        annot.samples   = as_tibble(samples_m18),
                        covar.matrix    = as.matrix(covar_m18),
                        covar.info      = as_tibble(covar.info),
                        data            = list(norm = as.matrix(norm_m18),
                                               raw  = as.matrix(otu_m18)),
                        datatype        = 'phenotype',
                        display.name    = 'DOMA OTU Abundance Month 18',
                        lod.peaks       = list(),
                        taxa            = as_tibble(taxa_m18))




### Markers as tibble
markers <- as_tibble(markers)








### Reducing genoprobs for lower memory. Then compute kinship
genoprobs <- probs_qtl2_to_doqtl(probs = genoprobs)
genoprobs <- genoprobs[dimnames(genoprobs)[[1]] %in% samples$mouse.id,,]
genoprobs <- probs_doqtl_to_qtl2(probs = genoprobs, map = as.data.frame(markers), marker_column = 'marker.id', pos_column = 'pos')

K <- calc_kinship(probs = genoprobs, type = 'loco', cores = 0)









### Save
rm(list = ls()[!grepl('dataset[.]|K|map|markers|genoprobs', ls())])
save.image(file = '~/Desktop/Weinstock_DOMA/Viewer/Version 2 - Centered Log Ratio/weinstock_doma_by_months_viewer_v2.Rdata')




