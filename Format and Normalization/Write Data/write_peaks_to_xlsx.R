library(openxlsx)


load('~/Desktop/Weinstock_DOMA/Viewer/Version 1 - VST and RankZ/weinstock_doma_viewer_v1.Rdata')







### Write peaks to .xlsx begin
for(i in grep('dataset[.]', ls(), value = TRUE)){
  
    # Get peaks
    age_peaks <- get(i)$lod.peaks$age_int
    sex_peaks <- get(i)$lod.peaks$sex_int
    add_peaks <- get(i)$lod.peaks$additive
    
    
    
    
    # Create workbook
    wb <- createWorkbook()
    
    addWorksheet(wb, 'additive')
    writeData(wb, 'additive', add_peaks)
    addWorksheet(wb, 'age_int')
    writeData(wb, 'age_int', age_peaks)
    addWorksheet(wb, 'sex_int')
    writeData(wb, 'sex_int', sex_peaks) 
    
    
    
    
    # Save
    taxa <- strsplit(i, split = '.', fixed = TRUE)[[1]][3]
    saveWorkbook(wb, file = paste0('doma_', taxa, '_qtl_peaks_summary.xlsx'))
    rm(wb)
    

}




