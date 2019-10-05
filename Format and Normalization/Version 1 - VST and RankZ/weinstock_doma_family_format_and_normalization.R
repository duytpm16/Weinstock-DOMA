### Options and libraries
options(stringsAsFactors = FALSE)
library(qtl2convert)
library(tidyverse)
library(DESeq2)
library(qtl2)






### Load data
samples <- read.csv('~/Desktop/Weinstock_DOMA/Phenotypes/do_mice/DO_CrossSectional_Population.csv')
chrY_M  <- read.csv('~/Desktop/Weinstock_DOMA/Phenotypes/do_mice/JAC_crosssectional_sex_chrM_Y_20180618.csv')
load('family_raw_count_and_taxa.Rdata')





### Removing family with <= 5% prevalence according to Hoan
#     family: 403 x 380
family_counts <- family_counts[, colSums(family_counts > 0) > (nrow(family_counts) * 0.05)]
family_taxa   <- family_taxa[,ncol(family_taxa):1]
family_taxa   <- family_taxa %>% filter(Family %in% colnames(family_counts))



### Create sample annotations
samples <- samples %>%
              filter(Cohort == 'Cross-Sectional') %>%
              filter(!grepl('HDO-', Mouse.ID)) %>%
              select(Mouse.ID, Sex, DOB, Generation, Cohort.Age, Coat.Color, Wean.Date, DOD, Age.Death, From.Cage, JCMS.Cage, Cage, Ear.Notch, Cohort) %>%
              mutate(Mouse.ID = gsub('-', '.', Mouse.ID)) %>%
              left_join(y = chrY_M %>% select(X, chrM, chrY) %>% dplyr::rename(Mouse.ID = X), by = 'Mouse.ID') %>%
              filter(Mouse.ID %in% rownames(family_counts)) %>%
              mutate(Sex = factor(Sex), Generation = factor(Generation)) %>%
              `colnames<-`(tolower(colnames(.))) %>%
              distinct()








### Normalize family. 
form <- formula(~ 1)
dds  <- DESeqDataSetFromMatrix(countData = t(family_counts), colData  = samples, design = form) 
vst  <- t(as.matrix(assay(varianceStabilizingTransformation(dds))))







### Rank Z transformation function
rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
} # rankZ()


rz <- apply(vst,  2, rankZ)








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
annot.phenotype <- data.frame(data.name   = c(colnames(samples), family_taxa[,'Family']),
                              short.name  = c(colnames(samples), family_taxa[,'Family']),
                              R.name      = c(colnames(samples), family_taxa[,'Family']),
                              description = c('Mouse identifier', 'Gender of mouse: F (female) or M (male)', 'Mouse date-of-birth', 'Generation of mouse', 'Age of mouse in months',
                                              'Color coat of mouse', 'Date of wean', 'Date of death', 'Age at death in days', 'Cage of mouse', 'JCMS cage',
                                              'Cage of mouse', 'Ear notch', 'JAC Cohort', 'chrM', 'chrY', apply(family_taxa, 1, paste0, collapse = ' - ')),
                              units       = NA,
                              category    = c(rep('Demographic', ncol(samples)), rep('Microbiome', nrow(family_taxa))),
                              R.category  = c(rep('Demographic', ncol(samples)), rep('Microbiome', nrow(family_taxa))),
                              is.id       = c(TRUE, rep(FALSE, ncol(samples) + nrow(family_taxa) - 1)),
                              is.numeric  = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
                                              rep(TRUE, nrow(family_taxa))),
                              is.date     = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, rep(FALSE, ncol(samples) + nrow(family_taxa) - 8)),
                              is.factor   = c(FALSE, TRUE, FALSE, TRUE, rep(FALSE, ncol(samples) + nrow(family_taxa) - 4)),
                              factor.levels = c(NA, 'F:M', NA, '8:9:10:11:12', rep(FALSE, ncol(samples) + nrow(family_taxa) - 4)),
                              is.covar    = c(FALSE, TRUE, FALSE, TRUE, TRUE, rep(FALSE, ncol(samples) + nrow(family_taxa) - 5)),
                              is.pheno    = c(rep(FALSE, ncol(samples)), rep(TRUE, nrow(family_taxa))),
                              is.derived  = FALSE,
                              omit        = FALSE,
                              use.covar   = c(rep(NA, ncol(samples)), rep('sex:cohort.age:generation', nrow(family_taxa))),
                              transformation = c(rep(NA, ncol(samples)), rep('rankz', nrow(family_taxa))))









### QTL viewer format
dataset.doma.family <- list(annot.phenotype = as_tibble(annot.phenotype),
                            annot.samples   = as_tibble(samples),
                            covar.matrix    = as.matrix(covar),
                            covar.info      = as_tibble(covar.info),
                            data            = list(norm = as.matrix(vst),
                                                   raw  = as.matrix(family_counts),
                                                   rz   = as.matrix(rz)),
                            datatype        = 'phenotype',
                            display.name    = 'DOMA Family Abundance',
                            lod.peaks       = list(),
                            taxa            = as_tibble(family_taxa))









### Save
rm(list = ls()[!grepl('dataset[.]', ls())])
load('~/Desktop/weinstock_doma_viewer_v1.Rdata')
save.image(file = '~/Desktop/weinstock_doma_viewer_v1.Rdata')
