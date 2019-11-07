###############################################################################
#
#  This script is used to format the genoprobs and markers downloaded from dodb
#     for qtl2. I downloaded the genoprobs on 10/03/2019
#
#  Input:
#    1.) genoprobs downloaded from dodb
#    2.) physical markers from dodb
#
#  Output:
#    1.) .Rdata containing genoprobs, map list, markers dataframe, and kinship list
#
#
#
#  Author: Duy Pham
#  E-mail: duy.pham@jax.org
#
################################################################################


### Options and libraries
options(stringsAsFactors = FALSE)
library(qtl2convert)
library(dplyr)
library(readr)       
library(qtl2)





### Load data
#    1.)Physical markers
#    2.) Read in DODB genoprobs data and convert to doqtl format for easy editing
markers   <- read.csv('~/Desktop/Weinstock_DOMA/Genotypes/DODB/Original/Churchill-046_JAC_DO_Aging-MegaMUGA_pmap.csv')
genoprobs <- read_rds('~/Desktop/Weinstock_DOMA/Genotypes/DODB/Original/Churchill_046_JAC_DO_Aging__genoprobs_8state_MegaMUGA.rds')
genoprobs <- probs_qtl2_to_doqtl(probs = genoprobs)








### Edit markers and convert to map list for qtl2
markers <- markers %>%
             filter(chr %in% c(1:19,'X')) %>%
             filter(pos > 0) %>%
             dplyr::rename(marker.id = marker) %>%
             mutate(chr = factor(chr, levels = c(1:19,'X')))

map <- map_df_to_list(map = markers, chr_column = 'chr', pos_column = 'pos', marker_column = 'marker.id')









### Remove JDO and HDO mice. Remove extra characters in sample IDs
genoprobs <- genoprobs[grep('_DO-', dimnames(genoprobs)[[1]]),,]

orig.id <- dimnames(genoprobs)[[1]]
dimnames(genoprobs)[[1]] <- gsub('Mouse_03oct2013_|Jackson_Laboratory_MEGMUGV01_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_', '', dimnames(genoprobs)[[1]])
dimnames(genoprobs)[[1]] <- gsub('Jackson_Laboratory_Churchill_MEGMUGV01_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_', '', dimnames(genoprobs)[[1]])
dimnames(genoprobs)[[1]] <- gsub("(DO-[0-9][0-9][0-9][0-9]).*", "\\1", dimnames(genoprobs)[[1]])
dimnames(genoprobs)[[1]] <- gsub("-", ".", dimnames(genoprobs)[[1]])






### Removing the duplicated samples that were genotyped at an earlier time.
# Will keep the most recent of the duplicated samples:
# [1] "Jackson_Laboratory_Churchill_MEGMUGV01_20150116_DO-0420_F5" 
# [2] "Jackson_Laboratory_Churchill_MEGMUGV01_20150116_DO-0392_B8" 
# [3] "Jackson_Laboratory_Churchill_MEGMUGV01_20150116_DO-0393_E7" 
# [4] "Jackson_Laboratory_Churchill_MEGMUGV01_20150116_DO-0397_H6" 
# [5] "Jackson_Laboratory_Churchill_MEGMUGV01_20150116_DO-0404_G10"
# [6] "Jackson_Laboratory_Churchill_MEGMUGV01_20150116_DO-0176_G5" 
# [7] "Jackson_Laboratory_Churchill_MEGMUGV01_20150116_DO-0177_B7" 
# [8] "Jackson_Laboratory_Churchill_MEGMUGV01_20150116_DO-0208_D7" 
# [9] "Jackson_Laboratory_Churchill_MEGMUGV01_20150116_DO-0408_H1" 
# [10] "Jackson_Laboratory_Churchill_MEGMUGV01_20150116_DO-0409_H2" 
# [11] "Jackson_Laboratory_MEGMUGV01_20140627_DO-0914_C10"          
# [12] "Jackson_Laboratory_MEGMUGV01_20140627_DO-0940_B11"
# This can be done as on line 92 due to the ordering of the samples.
genoprobs <- genoprobs[!duplicated(dimnames(genoprobs)[[1]]),,]








### Genoprobs back to qtl2 format
genoprobs <- probs_doqtl_to_qtl2(probs = genoprobs, map = markers, chr_column = 'chr', pos_column = 'pos', marker_column = 'marker.id')
K <- calc_kinship(probs = genoprobs, type = 'loco', cores = 0)








### Save
save.image(file = 'JAC_megaMUGA_genoprobs_qtl2.RData')



