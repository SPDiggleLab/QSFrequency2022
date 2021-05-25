## 05.25.2020 - MM
## This script loads all nucleotide sequence inputs, extracts the unique identifiers of mutated strains (ACCESSION IDs), and sorts to remove duplicate entries.
## A table is generated to output the numberand percentage of strains with nucleotide mutations for each gene.
## The sorted lists can then be uploaded into InteractiVenn to generate comparative diagrams.

library(readr)

# Load all nucleotide input files-----------------------------------------
dir.create('INPUT/nucleotide_inputs/accessions_only')
nuc_dir = "INPUT/nucleotide_inputs"
nuc_files = list.files(path=nuc_dir, pattern="*.csv", full.names=TRUE)


# Edit strain lists for each gene to work with----------------------------
for(nuc_file in nuc_files){
  
  # Load input file for each gene
  nuc_df = data.frame(read_csv(nuc_file, col_names = F), stringsAsFactors = FALSE)

  # Remove all but first 4 characters from each line in each file---------
  nuc_df$X1 <- substr(nuc_df$X1, 0, 4)

  # Remove duplicate IDs from each file-----------------------------------
  nuc_df<-unique(nuc_df)
  
  # Append gene name to beginning of file---------------------------------
  gene_name <- basename(nuc_file)
  gene <- strsplit(gene_name, '\\.')[[1]][1]

  # Save edited list as temp file-----------------------------------------
  temp_outName = paste('INPUT/nucleotide_inputs/accessions_only/', gene,'.csv', sep='', collapse = '')
  write.table(nuc_df, temp_outName, col.names = F, row.names = F, sep = ",")
}



# Load edited nucleotide input files--------------------------------------
dir.create('OUTPUT/nucleotide_output')
nuc_tempdir = "INPUT/nucleotide_inputs/accessions_only"
nuc_tempfiles = list.files(path=nuc_tempdir, pattern="*.csv", full.names=TRUE)


# Create a data frame to input gene data----------------------------------
strains_df <- data.frame(x = character(0), y = numeric(0), y =numeric(0), stringsAsFactors = FALSE)
colnames(strains_df) <- c("GENE", "NUMBER_OF_STRAINS", "PERCENT_OF_DB")


for(nuc_tempfile in nuc_tempfiles){
  
  # read in data
  nuc_tempdf <- data.frame(read_csv(nuc_tempfile, col_names = F), stringsAsFactors = FALSE)
  
  # Get the gene name from file name--------------------------------------
  gene_name <- basename(nuc_tempfile)
  gene <- strsplit(gene_name, '\\.')[[1]][1]
  
  
  # Count no. of lines in temp file = no. of strains----------------------
  strain <- nrow(nuc_tempdf)
  
  # calculate % of DB for each gene---------------------------------------
  percent <- ((strain)/852)*100
  
  # Combine gene name, no. of strains, and percent into one table---------
  temp_df <- data.frame("GENE" = gene, "NUMBER_OF_STRAINS" = strain, "PERCENT_OF_DB" = percent, stringsAsFactors = TRUE)
  temp_df$NUMBER_OF_STRAINS <- strain
  temp_df$PERCENT_OF_DB <- percent
  temp_df$GENE <- gene
  
  # Iterate through each file and combine data into one frame-------------
  strains_df <- rbind(temp_df, strains_df)
  
}

out_Table = paste('OUTPUT/nucleotide_output/', "nucleotide_table",'.csv', sep='', collapse = '')
write.table(strains_df, out_Table, sep=',', row.names = F)


  
# Cleanup-----------------------------------------------------------------
rm(nuc_df, nuc_tempdf, strains_df, temp_df, gene, gene_name, nuc_dir, nuc_file, nuc_files, nuc_tempfile, nuc_tempfiles, nuc_tempdir,out_Table, percent, strain, temp_outName)


