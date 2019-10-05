### Options and libraries
options(stringsAsFactors = FALSE)
library(tidyverse)
library(phyloseq)





### Read in OTU and taxa data
otu  <- readRDS('~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/Modified/0.5/otu_raw_count_0.5_cleaned.rds')
taxa <- readRDS('~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/Modified/0.5/otu_taxa_table_0.5_cleaned.rds')








### Rename taxa columns and set rownames as OTU
taxa <- taxa %>% 
          column_to_rownames('OTU') %>%
          dplyr::rename(Genus = genus, 
                        Family = family, 
                        Order = order, 
                        Class = class, 
                        Phylum = phylum, 
                        Domain = domain) %>%
          select(Domain, Phylum, Class, Order, Family, Genus) %>%
          as.matrix()











### Convert to phyloseq table and then merge to phyloseq object
phylo_otu  <- otu_table(otu, taxa_are_rows = FALSE)
phylo_taxa <- tax_table(taxa)

phyloseq <- merge_phyloseq(phylo_otu, phylo_taxa)










### Get genus level counts
genus <- subset_taxa(physeq = phyloseq, Genus != 'NA')
genus <- tax_glom(physeq = genus, taxrank = "Genus")

genus_counts <- as(otu_table(genus), "matrix")
genus_taxa   <- as.data.frame(as(tax_table(genus), "matrix")) %>% select(Genus)

colnames(genus_counts) <- genus_taxa$Genus[match(colnames(genus_counts), rownames(genus_taxa))]
genus_taxa <- as.data.frame(taxa) %>% 
                  filter(Genus %in% colnames(genus_counts)) %>% 
                  distinct() %>% 
                  arrange(Genus)

save(genus_counts, genus_taxa, file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/Modified/0.5/genus_raw_count_and_taxa_0.5.Rdata')









### Get family level counts
family <- subset_taxa(physeq = phyloseq, Family != 'NA')
family <- tax_glom(physeq = family, taxrank = "Family")

family_counts <- as(otu_table(family), "matrix")
family_taxa   <- as.data.frame(as(tax_table(family), "matrix")) %>% select(Family)

colnames(family_counts) <- family_taxa$Family[match(colnames(family_counts), rownames(family_taxa))]
family_taxa <- as.data.frame(taxa) %>% 
                    filter(Family %in% colnames(family_counts)) %>% 
                    select(-Genus) %>%
                    distinct() %>% 
                    arrange(Family)

save(family_counts, family_taxa, file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/Modified/0.5/family_raw_count_and_taxa_0.5.Rdata')











### Get order level counts
order <- subset_taxa(physeq = phyloseq, Order != 'NA')
order <- tax_glom(physeq = order, taxrank = "Order")

order_counts <- as(otu_table(order), "matrix")
order_taxa   <- as.data.frame(as(tax_table(order), "matrix")) %>% select(Order)

colnames(order_counts) <- order_taxa$Order[match(colnames(order_counts), rownames(order_taxa))]
order_taxa <- as.data.frame(taxa) %>% 
                  filter(Order %in% colnames(order_counts)) %>% 
                  select(-Genus, -Family) %>%
                  distinct() %>% 
                  arrange(Order)

save(order_counts, order_taxa, file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/Modified/0.5/order_raw_count_and_taxa_0.5.Rdata')










### Get class level counts
class <- tax_glom(physeq = phyloseq, taxrank = "Class")

class_counts <- as(otu_table(class), "matrix")
class_taxa   <- as.data.frame(as(tax_table(class), "matrix")) %>% select(Class)

colnames(class_counts) <- class_taxa$Class[match(colnames(class_counts), rownames(class_taxa))]
class_taxa <- as.data.frame(taxa) %>% 
                  filter(Class %in% colnames(class_counts)) %>% 
                  distinct(Class) %>% 
                  arrange(Class)
class_taxa <- merge(class_taxa, taxa[,c('Class', 'Phylum', 'Domain')], by = 'Class')
class_taxa <- class_taxa %>% distinct()

save(class_counts, class_taxa, file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/Modified/0.5/class_raw_count_and_taxa_0.5.Rdata')











### Get phylum level counts
phylum <- tax_glom(physeq = phyloseq, taxrank = "Phylum")

phylum_counts <- as(otu_table(phylum), "matrix")
phylum_taxa   <- as.data.frame(as(tax_table(phylum), "matrix")) %>% select(Phylum)

colnames(phylum_counts) <- phylum_taxa$Phylum[match(colnames(phylum_counts), rownames(phylum_taxa))]
phylum_taxa <- as.data.frame(taxa) %>% 
                  filter(Phylum %in% colnames(phylum_counts)) %>% 
                  select(-Genus, -Family, -Order, -Class) %>%
                  distinct(Phylum) %>% 
                  arrange(Phylum)
phylum_taxa$Domain <- 'Bacteria'

save(phylum_counts, phylum_taxa, file = '~/Desktop/Weinstock_DOMA/Phenotypes/doma_otu_16s_data/Modified/0.5/phylum_raw_count_and_taxa_0.5.Rdata')
