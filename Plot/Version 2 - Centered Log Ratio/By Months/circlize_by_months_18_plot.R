### Library
library(tidyverse)
library(circlize)





### Helper edit function
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






### Filter peaks to QTLs with LOD > 7.5 and edit for circlize
otu_m18    <- edit_peaksTable(dataset.otu.m18$lod.peaks$additive,    lod_thres = 7.5)
genus_m18  <- edit_peaksTable(dataset.genus.m18$lod.peaks$additive,  lod_thres = 7.5)
family_m18 <- edit_peaksTable(dataset.family.m18$lod.peaks$additive, lod_thres = 7.5)
order_m18  <- data.frame(chr = paste0('chr', 1), start = 0, end = 0, value1 = c(6.5,6.7)) # Dummy for plotting. Wont show up in plot
class_m18  <- data.frame(chr = paste0('chr', 1), start = 0, end = 0, value1 = c(6.5,6.7)) # Dummy for plotting. Wont show up in plot
phylum_m18 <- edit_peaksTable(dataset.phylum.m18$lod.peaks$additive, lod_thres = 7.2)






### Initialize circos
circos.par("start.degree" = 70, "gap.degree" = c(rep(0, 19), 40))
circos.initializeWithIdeogram(species = 'mm10', chromosome.index = paste0('chr',c(1:19,'X')))





### OTU track
circos.genomicTrackPlotRegion(otu_m18, track.height = .1, bg.border = NA, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, cex = 1.2, pch = 21, col = 'grey', bg = 'black', ...)
})

circos.yaxis(side = "left", at = seq(0,20,4), sector.index = 'chr1', labels.cex = .6)
for(i in paste0('chr', c(1:19,'X'))){
    circos.xaxis(h = 'bottom', sector.index = i, labels = FALSE)
}





### Genus track
circos.genomicTrackPlotRegion(genus_m18, track.height = .1, bg.border = NA, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, cex = 1.2, pch = 21, col = 'grey', bg = 'darkseagreen3', ...)
})

circos.yaxis(side = "left", at = seq(0,24,4), sector.index = 'chr1', labels.cex = .6)
for(i in paste0('chr', c(1:19,'X'))){
    circos.xaxis(h = 'bottom', sector.index = i, labels = FALSE)
}




### Family track
circos.genomicTrackPlotRegion(family_m18, track.height = .1, bg.border = NA,, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, cex = 1.2, pch = 21, col = 'grey', bg = 'firebrick3', ...)
})

circos.yaxis(side = "left", at = seq(0,12,2), sector.index = 'chr1', labels.cex = .6)
for(i in paste0('chr', c(1:19,'X'))){
    circos.xaxis(h = 'bottom', sector.index = i, labels = FALSE)
}






### Order track (colors are white because there were no peaks > 7.5)
circos.genomicTrackPlotRegion(order_m18, track.height = .1, bg.border = NA, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, cex = 1.2, pch = 21, col = 'white', bg = 'white', ...)
})

circos.yaxis(side = "left", at = seq(0,12,1), sector.index = 'chr1', labels.cex = .6)
for(i in paste0('chr', c(1:19,'X'))))
    circos.xaxis(h = 'bottom', sector.index = i, labels = FALSE)
}





### Class track (colors are white because there were no peaks > 7.5)
circos.genomicTrackPlotRegion(class_m18, track.height = .1, bg.border = NA,  panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, cex = 1.2, pch = 21, col = 'white', bg = 'white', ...)
})

circos.yaxis(side = "left", at = seq(0,12,2), sector.index = 'chr1', labels.cex = .6)
for(i in paste0('chr', c(1:19,'X'))){
    circos.xaxis(h = 'bottom', sector.index = i, labels = FALSE)
}






### Phylum track
circos.genomicTrackPlotRegion(phylum_m18, track.height = .1, bg.border = NA, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, cex = 1.2, pch = 21, col = 'grey', bg = 'darkorchid3', ...)
})

circos.yaxis(side = "left", at = seq(0,10,2), sector.index = 'chr1', labels.cex = .6)
for(i in paste0('chr', c(1:19,'X'))){
    circos.xaxis(h = 'bottom', sector.index = i, labels = FALSE)
}





### Add text
text(0,.25, 'Phylum', col = 'darkorchid3', font = 2)
text(0,.36, 'Class',  col = 'darkorange', font = 2)
text(0,.47, 'Order',  col = 'dodgerblue3', font = 2)
text(0,.70, 'Genus',  col = 'darkseagreen3', font = 2)
text(0,.58, 'Family', col = 'firebrick3', font = 2)
text(0,.80, 'OTU', font = 2)
text(0,1, 'Month 18', cex = 3)
