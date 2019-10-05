### Options and libraries
options(stringsAsFactors = FALSE)
library(tidyverse)
library(qtl2)





### Read in OTU data
taxa_0.8 <- read.delim(file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/otu_taxa_rdpc0.8.tsv')
taxa_0.5 <- read.delim(file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/otu_taxa_rdpc0.5.tsv')
otu_0.8  <- read.delim(file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/otu_table_by_otutab_0.8.tsv')
otu_0.5  <- read.delim(file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/DOMA_otu_403_0.5.txt')







### Editing taxa table
taxa_0.8 <- taxa_0.8 %>% 
              dplyr::rename(OTU = X) %>%
              mutate(OTU = gsub('_', '-', OTU)) %>%
              select(OTU, genus, family, order, class, phylum, domain)
taxa_0.8 <- apply(taxa_0.8, 2, function(x) gsub('_', ' ', x))
taxa_0.8 <- as.data.frame(taxa_0.8)



taxa_0.5 <- taxa_0.5 %>% 
              dplyr::rename(OTU = X) %>%
              mutate(OTU = gsub('_', '-', OTU)) %>%
              select(OTU, genus, family, order, class, phylum, domain)
taxa_0.5 <- apply(taxa_0.5, 2, function(x) gsub('_', ' ', x))
taxa_0.5 <- as.data.frame(taxa_0.5)










### Editing OTU 0.8 table
otu_0.8 <- otu_0.8 %>%
              mutate(X.OTU.ID = gsub('_', '-', X.OTU.ID)) %>%
              column_to_rownames('X.OTU.ID')


orig.id <- colnames(otu_0.8)
colnames(otu_0.8) <- gsub('DOMA_J00[0-9][A-Z][A-Z]_1_ST_T0_B0_0000_|DOMA_J00[0-9][A-Z][0-9]_1_ST_T0_B0_0000_|DOMA_J00[0-9][0-9][A-Z]_1_ST_T0_B0_0000_','',colnames(otu_0.8))
colnames(otu_0.8) <- gsub('(DO_[0-9][0-9][0-9][0-9]).*', '\\1', colnames(otu_0.8))
colnames(otu_0.8) <- gsub('(D0_[0-9][0-9][0-9][0-9]).*', '\\1', colnames(otu_0.8))
colnames(otu_0.8) <- gsub('D0', 'DO', colnames(otu_0.8))



### Removing duplicates. Keep the one with with highest abundance according to Dong-Binh
dups <- unique(colnames(otu_0.8)[duplicated(colnames(otu_0.8))])
keep <- rep(TRUE, ncol(otu_0.8))
for(i in dups){
  
    wh <- which(colnames(otu_0.8) == i)
    lessZeros <- which.min(colSums(otu_0.8[,wh] == 0))
    totalSum  <- which.max(colSums(otu_0.8[,wh]))
  
    keep[wh[-totalSum]] <- FALSE
}
otu_0.8 <- otu_0.8[,keep]
orig.id <- orig.id[keep]











### Editing OTU 0.5 table
otu_0.5 <- otu_0.5 %>% 
             dplyr::mutate(SampleID = gsub('_', '-', SampleID)) %>%  
             arrange(SampleID) %>%
             column_to_rownames('SampleID')
colnames(otu_0.5) <- gsub('(DO_[0-9][0-9][0-9][0-9]).*', "\\1", colnames(otu_0.5))
colnames(otu_0.5) <- gsub('_', ".", colnames(otu_0.5))








### Reorder
otu_0.8 <- otu_0.8[taxa_0.8$OTU, order(colnames(otu_0.8))]
colnames(otu_0.8) <- gsub('_', '.', colnames(otu_0.8))
otu_0.8 <- t(otu_0.8)


otu_0.5 <- otu_0.5[taxa_0.5$OTU, order(colnames(otu_0.5))]
otu_0.5 <- t(otu_0.5)








### Save
saveRDS(taxa_0.8, file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/Modified/otu_taxa_table_0.8_cleaned.rds')
saveRDS(taxa_0.5, file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/Modified/otu_taxa_table_0.5_cleaned.rds')
saveRDS(otu_0.5, file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/Modified/otu_raw_count_0.5_cleaned.rds')
saveRDS(otu_0.8, file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/Modified/otu_raw_count_0.8_cleaned.rds')
