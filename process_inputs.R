rm(list=ls())
## 04.11.2020 - CYZ
## This script goes through the Raw Inputs from KAO's BLAST search and standardizes the format.

require(readxl)



# (0) Find Inputs ---------------------------------------------------------

# Obtain list of available files. 
files = list.files('INPUT/Raw Inputs/', pattern = '.xlsx')

# Files denoted '_from_IPCD.xlsx' are sequences
ipcd_files = files[grepl('IPCD', files)]

# All others are metadata files denoted '[lab]-[environment].xlsx'
lab_files = files[!grepl('IPCD', files)]


# (1) Handle IPCD Files ---------------------------------------------------

# Load IPCD files with column names
warning("IPCD COLUMN NAMES ARE BACK INFERRED AS OF 04.11.2020")
ipcd_colnames = c('INFO', 'SEQUENCE', 'MISMATCH', 'GAPS', 'PCT_ID', 'QUERY_COV', 'BIT_SCORE', 'E-VALUE')
for(file in ipcd_files){
  temp = read_xlsx(paste('INPUT/Raw Inputs/', file, sep='', collapse=''), col_names = FALSE)
  temp = data.frame(temp, stringsAsFactors = F)
  colnames(temp) = ipcd_colnames
  
  # Write to a new file
  temp_outName = paste('INPUT/genes/', strsplit(file, '\\.')[[1]][1], '.tsv', sep='', collapse = '')
  write.table(temp, temp_outName, sep='\t', row.names = F)
}


# (2) Handle Lab Files ----------------------------------------------------

# Load Lab files with column names
warning("LAB COLUMN NAMES ARE BACK INFERRED AS OF 04.11.2020")
lab_colnames = c('STRAINID', 'COUNTRY', 'ANIMAL', 'SOURCE', 'ENV', 'PI', 'DETAILS')
for(file in lab_files){
  temp = read_xlsx(paste('INPUT/Raw Inputs/', file, sep='', collapse=''), col_names = FALSE)
  temp = data.frame(temp, stringsAsFactors = F)
  if(ncol(temp)==10){
    temp = temp[,-c(1,3,9)]
  }else{
    temp = temp[,-c(1,3)]
  }
  colnames(temp) = lab_colnames
  
  temp_outName = paste('INPUT/labs/', strsplit(file, '\\.')[[1]][1], '.tsv', sep='', collapse = '')
  write.table(temp, temp_outName, sep='\t', row.names = F)
}
