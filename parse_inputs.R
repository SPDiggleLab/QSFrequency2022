rm(list=ls())
## 04.11.2020 - CYZ
## This script loads all lab inputs and creates a master input data file.


#  (0) Helpers ------------------------------------------------------------

Fix_Known_Database_Discrepancies = function(strains){
  ret = strains
  ret = as.vector(gsub('MCF', 'AL', ret))
  return(ret)
}


# (1) Second Method -------------------------------------------------------

# Idea: Should be gene-focused. For each gene, add on the metadata information for the gene.
temp.ls_out = list()

# alternate labs
df.lab_strain_meta = read.table('INPUT/strain_meta.tsv', sep='\t', header = T, stringsAsFactors = F)

# if there are query duplicates, remove both.
v.duplicates = df.lab_strain_meta$STRAIN_ID[duplicated(df.lab_strain_meta$STRAIN_ID)]
df.lab_strain_meta = df.lab_strain_meta[!(df.lab_strain_meta$STRAIN_ID %in% v.duplicates),]

# Extract gene information
gene_files = list.files('INPUT/genes/', pattern='.tsv')
for(genefile in gene_files){
  gene = strsplit(genefile, '_')[[1]][1]
  temp.in = read.table(paste('INPUT/genes/', genefile, sep='', collapse=''), header=T, stringsAsFactors=F)
  
  # there are duplicate entries for fragments
  temp.in$INFO[duplicated(temp.in$INFO)] = paste(temp.in$INFO[duplicated(temp.in$INFO)], '_frag', sep='')
  rownames(temp.in) = temp.in$INFO 
  
  # raw_strain corrections: word 'strain' missing from info
  b.no_strain = !grepl('Pseudomonas aeruginosa strain', temp.in$INFO)
  temp.in$INFO[b.no_strain] = gsub('Pseudomonas aeruginosa', 'Pseudomonas aeruginosa strain', temp.in$INFO[b.no_strain])
  
  # Parse INFO for strain name
  raw_strains = temp.in$INFO
  
  # Obtain lab-reported strain ID
  parsed_strains = sapply(raw_strains, function(x){
    if(grepl('isolate', x)){
      split1 = strsplit(x, 'isolate ')[[1]][2]
      split2 = strsplit(split1, ' IPC')[[1]][1]
    }else{
      split1 = strsplit(x, 'strain ')[[1]][2]
      split2 = strsplit(split1, ' IPC')[[1]][1]
    }
    return(split2)
  })
  
  # Obtain NCBI ID
  alt_strain_id = sapply(raw_strains, function(x){
    split1 = strsplit(x, ' ')[[1]][1]
    return(split1)
  })
  
  # Obtain IPC ID
  IPCD_strain_id = sapply(raw_strains, function(x){
    split1 = strsplit(x, ',')[[1]][1]
    split2 = strsplit(split1, ' IPC')[[1]][2]
    return(paste('IPC', split2, sep='', collapse=''))
  })
  
  # Some parsed_strains have issues as well; this is a manual rule application
  fixed_strains = Fix_Known_Database_Discrepancies(parsed_strains)
  
  # Create a summary df for debugging
  df.strains_ref_by_gene = data.frame('STRAIN' = fixed_strains
                                      , 'ALT_ID' = alt_strain_id
                                      , 'IPCD_ID' = IPCD_strain_id
                                      , 'GENE' = gene
                                      , 'raw_ID' = raw_strains
                                      , 'parsed_ID' = parsed_strains
                                      , stringsAsFactors = F)
  rownames(df.strains_ref_by_gene) = df.strains_ref_by_gene$raw_ID
  
  # Combine 
  temp.in = cbind(temp.in, df.strains_ref_by_gene[temp.in$INFO, c('STRAIN', 'ALT_ID', 'IPCD_ID', 'GENE')])
  
  # Setup new columns for combination
  temp.in[,colnames(df.lab_strain_meta)] = NA
  
  errlist = cbind('gene' = rep(NA, nrow(temp.in)), 'errtype' = rep(NA, nrow(temp.in)))
  
  # query lab_strain_meta the careful way
  for(i in 1:nrow(temp.in)){
    # Expect only one query per temp.in entry
    i_strain = temp.in[i,'STRAIN']
    b.strain_meta = (df.lab_strain_meta$STRAIN_ID == i_strain)
    
    # If no exact match, use grep to identify single entry with closest match
    if(sum(b.strain_meta)==0){
      b.strain_meta = grepl(i_strain, df.lab_strain_meta$STRAIN_ID)
      errlist[i,] = c(i_strain, 'grepl used')
    }
    if(sum(b.strain_meta)== 1){
      temp.in[i,colnames(df.lab_strain_meta)] = df.lab_strain_meta[b.strain_meta,]
    }else{
      errlist[i,] = c(i_strain, 'multi-hit')
    }
  }
  
  # Save info
  temp.ls_out[[gene]] = temp.in
}

df.masterList = do.call('rbind', temp.ls_out)

# output
write.csv(df.masterList,  paste('OUTPUT/',Sys.Date(), '_master-table.csv', sep='', collapse=''), row.names = F, na = '')
