rm(list=ls())

Parse_Strain = function(strains){
  ret = strains
  ret = as.vector(gsub('MCF', 'AL', ret))
  b.iso = ret[grepl('isolate', ret)]
  
  print(ret[b.iso])
  
  temp.ret = sapply(ret[b.iso], function(x){
    ret2 = strsplit(x, ' isolate ')[[1]][2]
    print(ret2)
    return(ret2)
  })
  
  print(temp.ret)
  
  ret[b.iso] = temp.ret
  
  return(ret)
}

# load genes and labs
lab_files = list.files('INPUT/labs/', pattern='.tsv')

# Create df.labs
ls.labs = list()
for(labfile in lab_files){
  temp.in = read.table(sep='\t', (paste('INPUT/labs/', labfile, sep='', collapse='')), header=T, stringsAsFactors=F)
  temp.in$DETAILS = temp.in$SOURCE
  temp.in$SOURCE = strsplit(substr(labfile, start=1, stop=nchar(labfile)-4), '-')[[1]][2]
  
  ls.labs[[labfile]] = temp.in
}
df.labs = do.call('rbind', ls.labs)
df.labs = df.labs[!duplicated(df.labs$STRAINID),]
rownames(df.labs) = df.labs$STRAINID


# alternate labs
df.labs2 = read.table('INPUT/strain_meta.tsv', sep='\t', header = T, stringsAsFactors = F)
df.labs2[grepl('MCF', df.labs2$STRAIN_ID),]

# Extract gene information
gene_files = list.files('INPUT/genes/', pattern='.tsv')
for(genefile in gene_files){
  gene = strsplit(genefile, '_')[[1]][1]
  temp.in = read.table(paste('INPUT/genes/', genefile, sep='', collapse=''), header=T, stringsAsFactors=F)

  # Subset only for strains in gene files
  strains = sapply(temp.in$INFO, function(x){
    split1 = strsplit(x, 'strain ')[[1]][2]
    split2 = strsplit(split1, ' IPC')[[1]][1]
  })
  
  strains = Parse_Strain(strains)
  
  names(strains) = NULL
  
  
  
  # Remove NA strains
  temp.in = temp.in[!is.na(strains),]
  strains = strains[!is.na(strains)]
  temp.in = temp.in[!duplicated(strains),]
  strains = strains[!duplicated(strains)]
  rownames(temp.in) = strains
  
  # Save info
  df.labs = cbind(df.labs, NA)
  colnames(df.labs)[ncol(df.labs)] = gene
  df.labs[strains,gene] = temp.in[strains, 'SEQUENCE']
}
df.labs = df.labs[!is.na(df.labs$lasI),]

# output
write.csv(df.labs, 'OUTPUT/master_table.csv', row.names = F, na = '')