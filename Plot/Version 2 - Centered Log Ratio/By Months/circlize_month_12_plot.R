### Options and libraries
library(tidyverse)
library(circlize)




### Load data
load('~/Desktop/Weinstock_DOMA/Viewer/Version 2 - Centered Log Ratio/weinstock_doma_by_months_viewer_v2.Rdata')






### Helper edit peaks function
edit_peaksTable <- function(peaks, lod_thres){
  
  peaks %>% 
    filter(lod > lod_thres) %>%
    select(qtl.chr, qtl.pos, lod) %>%
    mutate(qtl.chr = paste0('chr', qtl.chr),
           qtl.pos = round(qtl.pos * 1e6),
           end     = qtl.pos) %>%
    dplyr::rename(chr    = qtl.chr, 
                  start  = qtl.pos,
                  value1 = lod) %>%
    select(chr, start, end, value1)
  
}




### Filter peaks to LOD > 7.5 and edit for circlize plotting
otu_m12    <- edit_peaksTable(dataset.otu.m12$lod.peaks$additive,    lod_thres = 7.5)
genus_m12  <- edit_peaksTable(dataset.genus.m12$lod.peaks$additive,  lod_thres = 7.5)
family_m12 <- edit_peaksTable(dataset.family.m12$lod.peaks$additive, lod_thres = 7.5)
order_m12  <- edit_peaksTable(dataset.order.m12$lod.peaks$additive,  lod_thres = 7.5)
class_m12  <- edit_peaksTable(dataset.class.m12$lod.peaks$additive,  lod_thres = 7.5)
phylum_m12 <- edit_peaksTable(dataset.phylum.m12$lod.peaks$additive, lod_thres = 7.5)








### Initialize circos plot
circos.par("start.degree" = 70, "gap.degree" = c(rep(0, 19), 40))
circos.initializeWithIdeogram(species = 'mm10', chromosome.index = paste0('chr',c(1:19,'X')))





### OTU track
circos.genomicTrackPlotRegion(otu_m12, track.height = .1, bg.border = NA, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, cex = 1.2, pch = 21, col = 'grey', bg = 'black', ...)
})

circos.yaxis(side = "left", at = seq(0,20,4), sector.index = 'chr1', labels.cex = .6)
for(i in paste0('chr', c(1:19,'X'))){
    circos.xaxis(h = 'bottom', sector.index = i, labels = FALSE)
}





### Genus track
circos.genomicTrackPlotRegion(genus_m12, track.height = .1, bg.border = NA, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, cex = 1.2, pch = 21, col = 'grey', bg = 'darkseagreen3', ...)
})

circos.yaxis(side = "left", at = seq(0,24,4), sector.index = 'chr1', labels.cex = .6)
for(i in paste0('chr', c(1:19,'X'))){
    circos.xaxis(h = 'bottom', sector.index = i, labels = FALSE)
}





### Family track
circos.genomicTrackPlotRegion(family_m12, track.height = .1, bg.border = NA,, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, cex = 1.2, pch = 21, col = 'grey', bg = 'firebrick3', ...)
})

circos.yaxis(side = "left", at = seq(0,12,2), sector.index = 'chr1', labels.cex = .6)
for(i in paste0('chr', c(1:19,'X'))){
    circos.xaxis(h = 'bottom', sector.index = i, labels = FALSE)
}





### Order track
circos.genomicTrackPlotRegion(order_m12, track.height = .1, bg.border = NA, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, cex = 1.2, pch = 21, col = 'grey', bg = 'dodgerblue3', ...)
})

circos.yaxis(side = "left", at = seq(0,12,1), sector.index = 'chr1', labels.cex = .6)
for(i in paste0('chr', c(1:19,'X'))){
    circos.xaxis(h = 'bottom', sector.index = i, labels = FALSE)
}






### Class track
circos.genomicTrackPlotRegion(class_m12, track.height = .1, bg.border = NA,  panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, cex = 1.2, pch = 21, col = 'grey', bg = 'darkorange', ...)
})

circos.yaxis(side = "left", at = seq(0,12,2), sector.index = 'chr1', labels.cex = .6)
for(i in paste0('chr', c(1:19,'X'))){
    circos.xaxis(h = 'bottom', sector.index = i, labels = FALSE)
}






### Phylum track
circos.genomicTrackPlotRegion(phylum_m12, track.height = .1, bg.border = NA, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, cex = 1.2, pch = 21, col = 'grey', bg = 'darkorchid3', ...)
})

circos.yaxis(side = "left", at = seq(0,10,2), sector.index = 'chr1', labels.cex = .6)
for(i in paste0('chr', c(1:19,'X'))){
    circos.xaxis(h = 'bottom', sector.index = i, labels = FALSE)
}






### Add text
text(0, 0.25, 'Phylum', col = 'darkorchid3',   font = 2)
text(0, 0.36, 'Class',  col = 'darkorange',    font = 2)
text(0, 0.47, 'Order',  col = 'dodgerblue3',   font = 2)
text(0, 0.70, 'Genus',  col = 'darkseagreen3', font = 2)
text(0, 0.58, 'Family', col = 'firebrick3',    font = 2)
text(0, 0.80, 'OTU',    col = 'black',         font = 2)
text(0, 1,   'Month 12', cex = 3)
