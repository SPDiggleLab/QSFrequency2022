rm(list=ls())
## 02.22.2020 - CYZ
## This script goes through the raw inputs from KAO's BLAST search and standardizes the format.

require(readxl)

if(!grepl('Raw', getwd())){setwd("Raw Inputs/")}

files = list.files(pattern = '.xlsx')
ipcd_files = files[grepl('IPCD', files)]
lab_files = files[!grepl('IPCD', files)]

#add colnames to ipcds
ipcd_colnames = c('INFO', 'SEQUENCE', 'MISMATCH', 'GAPS', 'PCT_ID', 'QUERY_COV', 'BIT_SCORE', 'E-VALUE')
for(file in ipcd_files){
  temp = read_xlsx(file, col_names = FALSE)
  temp = data.frame(temp, stringsAsFactors = F)
  colnames(temp) = ipcd_colnames
  
  temp_outName = paste('genes/', strsplit(file, '\\.')[[1]][1], '.csv', sep='', collapse = '')
  write.table(temp, temp_outName, sep='\t', row.names = F)
}

#add colnames to labs
lab_colnames = c('STRAINID', 'COUNTRY', 'ANIMAL', 'SOURCE', 'ENV', 'PI', 'DETAILS')
for(file in lab_files){
  temp = read_xlsx(file, col_names = FALSE)
  temp = data.frame(temp, stringsAsFactors = F)
  if(ncol(temp)==10){
    temp = temp[,-c(1,3,9)]
  }else{
    temp = temp[,-c(1,3)]
  }
  colnames(temp) = lab_colnames
  
  temp_outName = paste('labs/', strsplit(file, '\\.')[[1]][1], '.csv', sep='', collapse = '')
  write.table(temp, temp_outName, sep='\t', row.names = F)
}
