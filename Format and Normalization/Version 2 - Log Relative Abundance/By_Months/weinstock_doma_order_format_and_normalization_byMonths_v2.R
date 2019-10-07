### Options and libraries
options(stringsAsFactors = FALSE)
library(qtl2convert)
library(microbiome)
library(tidyverse)
library(qtl2)








### Load data
samples <- read.csv('~/Desktop/Weinstock_DOMA/Phenotypes/do_mice/DO_CrossSectional_Population.csv')
chrY_M  <- read.csv('~/Desktop/Weinstock_DOMA/Phenotypes/do_mice/JAC_crosssectional_sex_chrM_Y_20180618.csv')
load('~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/Modified/0.5/order_raw_count_and_taxa_0.5.Rdata')



order <- order_counts
taxa   <- order_taxa
taxa   <- taxa[,ncol(taxa):1]



### Divide samples into motnhs
weeks <- split(samples, samples$Cohort.Age)









### Get samples by month
order_m6 <- order[rownames(order) %in% weeks[['6']]$JAC.Mouse.ID, ]
order_m6 <- order_m6[, colSums(order_m6 > 0) > (nrow(order_m6) * 0.05)]


order_m12 <- order[rownames(order) %in% weeks[['12']]$JAC.Mouse.ID, ]
order_m12 <- order_m12[, colSums(order_m12 > 0) > (nrow(order_m12) * 0.05)]


order_m18 <- order[rownames(order) %in% weeks[['18']]$JAC.Mouse.ID, ]
order_m18 <- order_m18[, colSums(order_m18 > 0) > (nrow(order_m18) * 0.05)]









### Normalize order. 
norm_m6  <- t(transform(t(order_m6),  transform = 'clr', target = 'sample'))
norm_m12 <- t(transform(t(order_m12), transform = 'clr', target = 'sample'))
norm_m18 <- t(transform(t(order_m18), transform = 'clr', target = 'sample'))








### Get taxa for each month
taxa_m6  <- taxa %>% filter(Order %in% colnames(order_m6))
taxa_m12 <- taxa %>% filter(Order %in% colnames(order_m12))
taxa_m18 <- taxa %>% filter(Order %in% colnames(order_m18))







### Create sample annotations
samples <- samples %>%
              filter(Cohort == 'Cross-Sectional') %>%
              filter(!grepl('HDO-', Mouse.ID)) %>%
              select(Mouse.ID, Sex, DOB, Generation, Cohort.Age, Coat.Color, Wean.Date, DOD, Age.Death, From.Cage, JCMS.Cage, Cage, Ear.Notch, Cohort) %>%
              mutate(Mouse.ID = gsub('-', '.', Mouse.ID)) %>%
              left_join(y = chrY_M %>% select(X, chrM, chrY) %>% dplyr::rename(Mouse.ID = X), by = 'Mouse.ID') %>%
              mutate(Sex = factor(Sex), Generation = factor(Generation)) %>%
              `colnames<-`(tolower(colnames(.))) %>%
              filter(mouse.id %in% c(rownames(order_m6), rownames(order_m12), rownames(order_m18))) %>%
              distinct()

samples_m6  <- samples %>% filter(mouse.id %in% rownames(order_m6))
samples_m12 <- samples %>% filter(mouse.id %in% rownames(order_m12))
samples_m18 <- samples %>% filter(mouse.id %in% rownames(order_m18))









### Covariates
covar <- model.matrix(~ sex + generation + cohort.age, data = samples)[, -1, drop = FALSE]
covar <- covar[,c('sexM', 'cohort.age', grep('generation', colnames(covar), value = TRUE))]
colnames(covar)[1] <- 'sex'
rownames(covar) <- samples$mouse.id

covar_m6  <- covar[rownames(covar) %in% rownames(order_m6),]
covar_m12 <- covar[rownames(covar) %in% rownames(order_m12),]
covar_m18 <- covar[rownames(covar) %in% rownames(order_m18),]



covar.info <- data.frame(sample.column = c('sex', 'cohort.age', 'generation'),
                         covar.column  = c('sex', 'cohort.age', 'generation'),
                         display.name  = c('Sex', 'Cohort Age', 'Generation'),
                         interactive   = c(TRUE, TRUE, FALSE),
                         primary       = c(TRUE, TRUE, FALSE),
                         lod.peaks     = c('sex_int', 'age_int', NA))












### Annotate columns
annot.phenotype <- data.frame(data.name   = c(colnames(samples), taxa[,'Order']),
                              short.name  = c(colnames(samples), taxa[,'Order']),
                              R.name      = c(colnames(samples), taxa[,'Order']),
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

annot.pheno.m6  <- annot.phenotype %>% filter(data.name %in% c(colnames(samples), colnames(order_m6)))
annot.pheno.m12 <- annot.phenotype %>% filter(data.name %in% c(colnames(samples), colnames(order_m12)))
annot.pheno.m18 <- annot.phenotype %>% filter(data.name %in% c(colnames(samples), colnames(order_m18)))











### QTL viewer format
dataset.order.m6 <- list(annot.phenotype = as_tibble(annot.pheno.m6),
                         annot.samples   = as_tibble(samples_m6),
                         covar.matrix    = as.matrix(covar_m6),
                         covar.info      = as_tibble(covar.info),
                         data            = list(norm = as.matrix(norm_m6),
                                                raw  = as.matrix(order_m6)),
                         datatype        = 'phenotype',
                         display.name    = 'DOMA Order Abundance Month 6',
                         lod.peaks       = list(),
                         taxa            = as_tibble(taxa_m6))


dataset.order.m12 <- list(annot.phenotype = as_tibble(annot.pheno.m12),
                          annot.samples   = as_tibble(samples_m12),
                          covar.matrix    = as.matrix(covar_m12),
                          covar.info      = as_tibble(covar.info),
                          data            = list(norm = as.matrix(norm_m12),
                                                 raw  = as.matrix(order_m12)),
                          datatype        = 'phenotype',
                          display.name    = 'DOMA Order Abundance Month 12',
                          lod.peaks       = list(),
                          taxa            = as_tibble(taxa_m12))


dataset.order.m18 <- list(annot.phenotype = as_tibble(annot.pheno.m18),
                          annot.samples   = as_tibble(samples_m18),
                          covar.matrix    = as.matrix(covar_m18),
                          covar.info      = as_tibble(covar.info),
                          data            = list(norm = as.matrix(norm_m18),
                                                 raw  = as.matrix(order_m18)),
                          datatype        = 'phenotype',
                          display.name    = 'DOMA Order Abundance Month 18',
                          lod.peaks       = list(),
                          taxa            = as_tibble(taxa_m18))









### Save
rm(list = ls()[!grepl('dataset[.]|K|map|markers|genoprobs', ls())])
load('~/Desktop/Weinstock_DOMA/Viewer/Version 2 - Centered Log Ratio/weinstock_doma_by_months_viewer_v2.Rdata')
save.image(file = '~/Desktop/Weinstock_DOMA/Viewer/Version 2 - Centered Log Ratio/weinstock_doma_by_months_viewer_v2.Rdata')
