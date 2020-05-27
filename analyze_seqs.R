rm(list=ls())
## 04.11.2020 - CYZ
## This script analyzes the master table

require(Biostrings)
require(factoextra)
require(ape)
require(vegan)

TRUNCATION_CUTOFF = 0

# (0) Helpers -------------------------------------------------------------

ALIGN_VS_REF_BIOSTRING = function(REF=lasI_ref, ALL=lasI_all, gapOpen=11, gapExt=1){
  ret = lapply(ALL, function(x){
    return (pairwiseAlignment(REF
                              , x
                              , substitutionMatrix = "BLOSUM80"
                              , gapOpening=gapOpen
                              , gapExtension=gapExt)@score) 
    })
  
  ret = do.call('rbind', ret)
  return(ret)
}#test_ALIGN_VS_REF_BIOSTRING = ALIGN_VS_REF_BIOSTRING(ls.aastring_map[[1]], ls.aastring_map)

ALIGN_VS_SELF_BIOSTRING = function(ALL=lasI_all){
  
  # to strings
  v.seqs = unlist(lapply(ALL, as.character)) # char
  v.seq_nos = names(ALL) # char

  names(v.seqs) = v.seq_nos
  
  nSeqs = length(v.seqs)

  # find uniques
  v.unique_seqs = v.seqs[!duplicated(v.seqs)]
  
  # Reduce to unique AA sequences
  ls_unique_seqs = ALL[names(v.unique_seqs)] # WARNING: the 'unique' function doesn't work with lists

  # Calculate distances on AA sequences BIOSTRINGS
  ls.uniqueAADist = lapply(ls_unique_seqs, function(x){
    ret = ALIGN_VS_REF_BIOSTRING(x, ls_unique_seqs)
  })
  
  # Format distance matrix of unique AA sequences
  df.uniqueAADist = do.call('cbind', ls.uniqueAADist)
  colnames(df.uniqueAADist) = v.unique_seqs
  rownames(df.uniqueAADist) = v.unique_seqs
  
  # Format distance matrix of unique NT sequences
  df.fullAADist = data.frame(matrix(rep(0, nSeqs^2), ncol = nSeqs))
  colnames(df.fullAADist) = v.seq_nos
  rownames(df.fullAADist) = v.seq_nos
  
  # Fill distance matrix of unique NT sequences
  for(i in 1:nSeqs){
    for(j in 1:nSeqs){
      i.seq = v.seqs[v.seq_nos[i]]
      j.seq = v.seqs[v.seq_nos[j]]
      
      if(!j.seq %in% v.unique_seqs){
        print(j.seq)
      }
      
      df.fullAADist[i,j] = df.uniqueAADist[i.seq, j.seq]
    }
  }
  
  return(df.fullAADist)
} #test_ALIGN_VS_SELF_BIOSTRING = ALIGN_VS_SELF_BIOSTRING(ls.aastring_map)

# (1) Load ----------------------------------------------------------------

# Get truncation statistics
if(!file.exists('OUTPUT/gene_trunc_stats.R')){
  
  # Read in sequence master table
  df.in = read.table('OUTPUT/2020-04-11_master-table.csv', header = T, sep=',', stringsAsFactors = F)
  df.in = df.in[,c('STRAIN', 'GENE', 'HOST', 'SOURCE', 'ENV', 'PI'
                   , 'SEQUENCE' , 'MISMATCH', 'GAPS', 'PCT_ID', 'QUERY_COV', 'BIT_SCORE', 'E.VALUE'
                   , 'INCLUDE')]
  
  # Fill in Blanks
  df.in[df.in==''] = NA
  
  # Must have PI
  b.noENV = is.na(df.in$ENV)
  b.noPI = is.na(df.in$PI)
  
  df.error = df.in[b.noENV | b.noPI,]
  df.in = df.in[!(b.noENV | b.noPI),]
  
  # Split by gene
  ls.in = split(df.in, df.in$GENE)
  ls.truncations = list()
  ls.complete_gene_MetaSeq = list()
  
  for(gene in names(ls.in)){
    # load Ref
    ref_gene = read.table(paste('INPUT/ref_genes/', gene, '_PAO1_protein.txt', sep='', collapse='')
                          , stringsAsFactors = F)$V1
    
    # query gene data
    df.gene = ls.in[[gene]]
    
    # Calculate reference gene length
    gene_length = nchar(ref_gene)*3 # protein sequence * 3
    
    # Get the number of genes used in the analysis
    n_genes = nrow(df.gene)
    
    # Get the number of unique genes used in the analysis
    n_unique_genes = n_genes - sum(duplicated(df.gene$STRAIN))
    
    # Get number of entries with '-' gaps
    b.gaps = grepl('-', df.gene$SEQUENCE)
    n_gaps = sum(b.gaps)
    
    # Remove gap-aligned sequences
    df.gene = df.gene[!b.gaps,]
    
    # Get number of truncations
    v.gene_lens = sapply(df.gene$SEQUENCE, function(x){nchar(x)}) 
    b.ref_length = (v.gene_lens < TRUNCATION_CUTOFF*gene_length)
    n_trunc = sum(b.ref_length)
    mean_length = mean(v.gene_lens)
    
    # Remove 
    df.gene = df.gene[!b.ref_length,]
    
    # Get number of 
    ls.truncations[[gene]] = c('gene' = gene
                               , 'n_base' = gene_length
                               , 'n_genes' = n_genes
                               , 'n_unique' = n_unique_genes
                               , 'p_incomplete' = (n_trunc+n_gaps)/n_genes
                               , 'mean_len' = mean_length/gene_length) # save stats
    ls.complete_gene_MetaSeq[[gene]] = df.gene # save t
  }
  
  df.truncations = data.frame(do.call('rbind', ls.truncations), stringsAsFactors = F)
  df.truncations[,-1] = data.matrix(df.truncations[,-1]) # fix datatype
  
  write.csv(df.truncations, file = 'OUTPUT/gene_truncated_stats.csv')
  
  save(ls.complete_gene_MetaSeq
       , df.truncations
       , TRUNCATION_CUTOFF
       , ALIGN_VS_REF_BIOSTRING
       , ALIGN_VS_SELF_BIOSTRING
       , file='OUTPUT/gene_trunc_stats.R')
  
}else{
  load('OUTPUT/gene_trunc_stats.R')
}



# (3) Calculate Distance Matrices -----------------------------------------

if(!file.exists('OBJECTS/list_of_strain_distMats_by_gene.R')){
  ls.distMats_by_gene = list()
  
  for(gene in names(ls.complete_gene_MetaSeq)){
    print(gene)
    
    temp.in = ls.complete_gene_MetaSeq[[gene]]
    temp.in$STRAIN[duplicated(temp.in$STRAIN)] = paste(temp.in$STRAIN[duplicated(temp.in$STRAIN)], '_2', sep='', collapse='')
    
    temp.in = temp.in[(sapply(temp.in$SEQUENCE, function(x){substr(x, 1,3)=='ATG'})),] # only keep entries that start with 'atg'    
    temp.in = temp.in[!duplicated(temp.in$STRAIN),] # remove duplicates
    
    # Sequence-to-Strain map
    v.STRAIN_TO_SEQ = temp.in$SEQUENCE
    names(v.STRAIN_TO_SEQ) = temp.in$STRAIN
    
    nStrains = length(v.STRAIN_TO_SEQ)
    
    # Remove redundant calculations & create Sequence-to-SeqNO map
    v.seq_map = unique(temp.in$SEQUENCE)
    v.seq_names = paste('SEQ_', 1:length(v.seq_map), sep='')
    names(v.seq_map) = v.seq_names
    
    ls.aastring_map = list()
    for(i in 1:length(v.seq_map)){
      temp.seq_name = v.seq_names[i]
      temp.str_in = v.seq_map[i]
      temp.dnastring = DNAString(temp.str_in)
      temp.aastring = translate(temp.dnastring, if.fuzzy.codon = c('solve', 'X'), no.init.codon = T)
      
      ls.aastring_map[[temp.seq_name]] = temp.aastring
    }
    
    df.distMat_by_seq = ALIGN_VS_SELF_BIOSTRING(ls.aastring_map)
    rownames(df.distMat_by_seq) = v.seq_map[rownames(df.distMat_by_seq)]
    colnames(df.distMat_by_seq) = v.seq_map[colnames(df.distMat_by_seq)]
    
    df.distMat_by_strain = data.frame(matrix(rep(0, nStrains^2), ncol = nStrains))
    rownames(df.distMat_by_strain) = names(v.STRAIN_TO_SEQ)
    colnames(df.distMat_by_strain) = names(v.STRAIN_TO_SEQ)
    
    # go by row to make it easier
    for(i_distMat_by_strain in 1:nrow(df.distMat_by_strain)){
      df.distMat_by_strain[i_distMat_by_strain,names(v.STRAIN_TO_SEQ)] = df.distMat_by_seq[v.STRAIN_TO_SEQ[i_distMat_by_strain], v.STRAIN_TO_SEQ]
    }
    
    ls.distMats_by_gene[[gene]] = df.distMat_by_strain
  }
  
  save(ls.distMats_by_gene, file='OBJECTS/list_of_strain_distMats_by_gene.R')
}else{
  load('OBJECTS/list_of_strain_distMats_by_gene.R')
}
