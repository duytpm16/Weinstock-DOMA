### Options and libraries
options(stringsAsFactors = FALSE)
library(tidyverse)
library(qtl2)





### Read in OTU data
taxa_0.8 <- read.delim(file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/otu_taxa_rdpc0.8.tsv')
taxa_0.5 <- read.delim(file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/otu_taxa_rdpc0.5.tsv')
otu  <- read.delim(file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/otu_table_by_otutab.tsv')







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










### Editing OTU table
otu <- otu %>%
          mutate(X.OTU.ID = gsub('_', '-', X.OTU.ID)) %>%
          column_to_rownames('X.OTU.ID')


orig.id <- colnames(otu)
colnames(otu) <- gsub('DOMA_J00[0-9][A-Z][A-Z]_1_ST_T0_B0_0000_|DOMA_J00[0-9][A-Z][0-9]_1_ST_T0_B0_0000_|DOMA_J00[0-9][0-9][A-Z]_1_ST_T0_B0_0000_','',colnames(otu))
colnames(otu) <- gsub('(DO_[0-9][0-9][0-9][0-9]).*', '\\1', colnames(otu))
colnames(otu) <- gsub('(D0_[0-9][0-9][0-9][0-9]).*', '\\1', colnames(otu))
colnames(otu) <- gsub('D0', 'DO', colnames(otu))








### Removing duplicates. Keep the one with with highest abundance according to Dong-Binh
dups <- unique(colnames(otu)[duplicated(colnames(otu))])
keep <- rep(TRUE, ncol(otu))
for(i in dups){
  
    wh <- which(colnames(otu) == i)
    lessZeros <- which.min(colSums(otu[,wh] == 0))
    totalSum  <- which.max(colSums(otu[,wh]))
  
    keep[wh[-totalSum]] <- FALSE
}
otu <- otu[,keep]
orig.id <- orig.id[keep]










### Reorder
otu <- otu[taxa_0.8$OTU, order(colnames(otu))]
colnames(otu) <- gsub('_', '.', colnames(otu))
otu <- t(otu)











### Save
saveRDS(taxa_0.8, file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/Modified/otu_taxa_table_0.8_cleaned.rds')
saveRDS(taxa_0.5, file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/Modified/otu_taxa_table_0.5_cleaned.rds')
saveRDS(otu, file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/Modified/otu_raw_count_cleaned')
