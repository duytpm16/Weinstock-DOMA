library(openxlsx)


#load('~/Desktop/Weinstock_DOMA/Viewer/Version 1 - VST and RankZ/weinstock_doma_viewer_v1.Rdata')







### Write heritability results to .xlsx begin
for(i in grep('dataset[.]', ls(), value = TRUE)){
  
    # Get peaks
    overall <- get(i)$herit$overall
    chromosome <- get(i)$herit$chromosome
    
    
    
    
    # Create workbook
    wb <- createWorkbook()
    
    addWorksheet(wb, 'Overall')
    writeData(wb, 'Overall', overall, rowNames = TRUE)
    addWorksheet(wb, 'By_Chromosome')
    writeData(wb, 'By_Chromosome', chromosome, rowNames = TRUE)
    
    
    
    
    # Save
    taxa <- strsplit(i, split = '.', fixed = TRUE)[[1]][3]
    saveWorkbook(wb, file = paste0('doma_', taxa, '_heritability_summary.xlsx'))
    rm(wb)
    

}




