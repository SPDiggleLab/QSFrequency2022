rm(list=ls())
## 04.11.2020 - CYZ
## This script goes through the RAW_INPUTS from KAO's BLAST search and standardizes the format.
## Outputs are saved in '/INPUT/genes' and '/INPUT/labs'

require(readxl)


# (0) Find Inputs ---------------------------------------------------------

# Obtain list of available excel files. 
files = list.files('INPUT/RAW_INPUTS/', pattern = '.xlsx')

# Files denoted '_from_IPCD.xlsx' are sequences
ipcd_files = files[grepl('IPCD', files)]

# All others are metadata files denoted '[lab]-[environment].xlsx'
lab_files = files[!grepl('IPCD', files) & !grepl('METADATA', files)] # one RAW_INPUT file is 'STRAIN_METADATA'

if(!dir.exists('INPUT/genes')){
  dir.create('INPUT/genes')
}
if(!dir.exists('INPUT/labs')){
  dir.create('INPUT/labs')
}

# (1) Handle IPCD Files ---------------------------------------------------

# Load IPCD files with column names
warning("IPCD COLUMN NAMES ARE BACK INFERRED AS OF 04.11.2020")
ipcd_colnames = c('INFO', 'SEQUENCE', 'MISMATCH', 'GAPS', 'PCT_ID', 'QUERY_COV', 'BIT_SCORE', 'E-VALUE')
for(file in ipcd_files){
  temp = read_xlsx(paste('INPUT/RAW_INPUTS/', file, sep='', collapse=''), col_names = FALSE)
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
  temp = read_xlsx(paste('INPUT/RAW_INPUTS/', file, sep='', collapse=''), col_names = FALSE)
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


# (3) Cleanup -------------------------------------------------------------

rm(temp, file, files, ipcd_colnames, ipcd_files, lab_colnames, lab_files, temp_outName)
