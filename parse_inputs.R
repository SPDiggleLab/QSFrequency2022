rm(list=ls())

# load genes and labs
lab_files = list.files('INPUT/labs/', pattern='.tsv')

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

gene_files = list.files('INPUT/genes/', pattern='.tsv')
for(genefile in gene_files){
  gene = strsplit(genefile, '_')[[1]][1]
  temp.in = read.table(paste('INPUT/genes/', genefile, sep='', collapse=''), header=T, stringsAsFactors=F)

  strains = sapply(temp.in$INFO, function(x){
    split1 = strsplit(x, 'strain ')[[1]][2]
    split2 = strsplit(split1, ' IPC')[[1]][1]
  })
  names(strains) = NULL
  
  temp.in = temp.in[!is.na(strains),]
  strains = strains[!is.na(strains)]
  temp.in = temp.in[!duplicated(strains),]
  strains = strains[!duplicated(strains)]
  
  rownames(temp.in) = strains
  
  df.labs = cbind(df.labs, NA)
  colnames(df.labs)[ncol(df.labs)] = gene
  
  df.labs[strains,gene] = temp.in[strains, 'SEQUENCE']
}
df.labs = df.labs[!is.na(df.labs$lasI),]

write.csv(df.labs, 'master_table.csv', row.names = F, na = '')
