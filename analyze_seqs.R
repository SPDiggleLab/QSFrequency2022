rm(list=ls())
## 04.11.2020 - CYZ
## This script analyzes the master table

require(Biostrings)
require(factoextra)
require(ape)
require(vegan)

TRUNCATION_CUTOFF = 0 #0 for all data. 1 for no truncations allowed.

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
if(!file.exists('OBJECTS/gene_trunc_stats.R')){
  
  # Find the most recent master table file
  temp.fileList = list.files('OUTPUT/')
  temp.fileDirs = temp.fileList[grep('_master-table.csv', temp.fileList)]
  temp.fileDir = sort(temp.fileDirs, decreasing = T)[1]
  
  # Read in sequence master table
  df.in = read.table(paste('OUTPUT/', temp.fileDir, sep=''), header = T, sep=',', stringsAsFactors = F)
  df.in = df.in[,c('STRAIN', 'GENE', 'HOST', 'SOURCE', 'ENV', 'PI'
                   , 'SEQUENCE' , 'MISMATCH', 'GAPS', 'PCT_ID', 'QUERY_COV', 'BIT_SCORE', 'E.VALUE'
                   , 'INCLUDE')]
  rm(temp.fileList, temp.fileDirs, temp.fileDir)
  
  # Fill in Blanks
  df.in[df.in==''] = NA
  
  # Split by gene
  ls.in = split(df.in, df.in$GENE)
  ls.truncations = list()
  ls.complete_gene_MetaSeq = list()
  for(gene in names(ls.in)){
    # load Ref
    ref_aminos = read.table(paste('INPUT/ref_genes/', gene, '_PAO1_protein.txt', sep='', collapse='')
                          , stringsAsFactors = F)$V1
    ref_aminos_length = nchar(ref_aminos)
    
    # query gene data
    df.gene = ls.in[[gene]]
    
    # How many genes pass quality filter?
    n.allGenes = nrow(df.gene)
    
    # Must have: PI, (ENV|(SOURCE|HOST))
    b.noENV = is.na(df.gene$ENV) # envcan't be NA
    b.noPI = is.na(df.gene$PI)
    df.error = df.gene[b.noENV | b.noPI,]
    df.gene = df.gene[!(b.noENV | b.noPI),]
    
    # If Not ENV, must have SOURCE|HOST
    b.ENV = df.gene$ENV
    b.noSOURCE = is.na(df.gene$SOURCE)
    b.noHOST = is.na(df.gene$HOST)
    df.gene = df.gene[(b.ENV|!(b.noSOURCE&b.noHOST)),]
    
    # Remove entries that do not start with 'ATG'
    b.atg = sapply(df.gene$SEQUENCE, function(x){substr(x, 1, 3)!='ATG'})
    df.error = rbind(df.error, df.gene[b.atg,])
    df.gene = df.gene[!b.atg,]
    
    # Remove entries with '-' gaps
    b.gaps = grepl('-', df.gene$SEQUENCE)
    n.gaps = sum(b.gaps)
    df.error = rbind(df.error, df.gene[b.gaps,])
    df.gene = df.gene[!b.gaps,]
    
    # Calculate reference gene length
    gene_length = nchar(ref_aminos)*3 # protein sequence * 3
    
    # Get the number of genes used in the analysis
    n.genes = nrow(df.gene)
    
    # Get the number of unique strains/seqs used in the analysis
    n.unique_strains = n.genes - sum(duplicated(df.gene$STRAIN))
    n.unique_sequences = n.genes - sum(duplicated(df.gene$SEQUENCE))

    v.gene_lens = sapply(df.gene$SEQUENCE, function(x){nchar(x)}) 
    b.ref_length_passed_cutoff = (v.gene_lens > TRUNCATION_CUTOFF*gene_length)
    n.trunc = sum(!b.ref_length_passed_cutoff)
    mean_length = mean(v.gene_lens)
    
    # Remove 
    df.gene = df.gene[b.ref_length_passed_cutoff,]
    
    # Get number of 
    ls.truncations[[gene]] = c('gene' = gene
                               , 'n_base' = gene_length
                               , 'n_allGenes' = n.allGenes
                               , 'n_genes' = n.genes
                               , 'n_unique_strains' = n.unique_strains
                               , 'n_unique_seqs' = n.unique_sequences
                               , 'p_incomplete' = (n.trunc+n.gaps)/n.genes
                               , 'mean_len' = mean_length/gene_length) # save stats
    ls.complete_gene_MetaSeq[[gene]] = df.gene # save t
  }
  
  df.truncations = data.frame(do.call('rbind', ls.truncations), stringsAsFactors = F)
  df.truncations[,-1] = apply(df.truncations[,-1], 2, as.numeric) # fix datatype
  
  write.csv(df.truncations, file = 'OUTPUT/gene_truncated_stats.csv')
  
  save(ls.complete_gene_MetaSeq
       , df.truncations
       , TRUNCATION_CUTOFF
       , ALIGN_VS_REF_BIOSTRING
       , ALIGN_VS_SELF_BIOSTRING
       , file='OBJECTS/gene_trunc_stats.R')
  
}else{
  load('OBJECTS/gene_trunc_stats.R')
}


# (3) Calculate Distance Matrices -----------------------------------------

if(!file.exists('OBJECTS/list_of_strain_distMats_by_gene.R')){
  ls.distMats_by_gene = list()
  ls.refDist_by_gene = list() 
  
  for(gene in names(ls.complete_gene_MetaSeq)){
    print(gene)
    
    # load Ref
    ref_aminos = read.table(paste('INPUT/ref_genes/', gene, '_PAO1_protein.txt', sep='', collapse='')
                          , stringsAsFactors = F)$V1
    ref_aminos = AAString(ref_aminos) # read in as biocstring object
    
    # Query for the relevant gene
    temp.in = ls.complete_gene_MetaSeq[[gene]]
    
    # For all strain duplicates, append "_2" to name.
    max_duplicated = max(table(temp.in$STRAIN))
    for(i.dupl in 1:max_duplicated){
      b.duplicated = duplicated(temp.in$STRAIN)
      temp.in$STRAIN[b.duplicated] = paste(temp.in$STRAIN[b.duplicated], '_', i.dupl+1, sep='', collapse='')
    }

    # Only keep entries that start with start codon.
    temp.in = temp.in[(sapply(temp.in$SEQUENCE, function(x){substr(x, 1,3)=='ATG'})),] # should be redundant.  
    temp.in = temp.in[!duplicated(temp.in$STRAIN),] # should be redundant
    
    # Sequence-to-Strain map
    v.STRAIN_TO_SEQ = temp.in$SEQUENCE
    names(v.STRAIN_TO_SEQ) = temp.in$STRAIN

    nStrains = length(v.STRAIN_TO_SEQ)
    
    # Remove redundant calculations & create Sequence-to-SeqNO map
    v.seq_map = unique(temp.in$SEQUENCE)
    v.seq_names = paste('SEQ_', 1:length(v.seq_map), sep='')
    names(v.seq_map) = v.seq_names
    
    # Reduce to unique calculations to avoid O(N^2) time
    ls.aastring_map = list()
    for(i in 1:length(v.seq_map)){
      temp.seq_name = v.seq_names[i]
      temp.str_in = v.seq_map[i]
      temp.dnastring = DNAString(temp.str_in)
      temp.aastring = translate(temp.dnastring, if.fuzzy.codon = c('solve', 'X'), no.init.codon = T)
      
      ls.aastring_map[[temp.seq_name]] = temp.aastring
    }
    
    # Calculate distances to ref
    df.refDist = ALIGN_VS_REF_BIOSTRING(REF = ref_aminos, ALL = ls.aastring_map)
    rownames(df.refDist) = v.seq_map[rownames(df.refDist)]
    
    df.refDist_by_strain = df.refDist[v.STRAIN_TO_SEQ,1]
    names(df.refDist_by_strain) = names(v.STRAIN_TO_SEQ)
    
    ls.refDist_by_gene[[gene]] = df.refDist_by_strain
    
    # Calculate pairwise distances
    df.distMat_by_seq = ALIGN_VS_SELF_BIOSTRING(ls.aastring_map)
    rownames(df.distMat_by_seq) = v.seq_map[rownames(df.distMat_by_seq)]
    colnames(df.distMat_by_seq) = v.seq_map[colnames(df.distMat_by_seq)]
    
    df.distMat_by_strain = data.frame(matrix(rep(0, nStrains^2), ncol = nStrains))
    rownames(df.distMat_by_strain) = names(v.STRAIN_TO_SEQ)
    colnames(df.distMat_by_strain) = names(v.STRAIN_TO_SEQ)
    
    for(i_distMat_by_strain in 1:nrow(df.distMat_by_strain)){
      # go by row to make it easier
      df.distMat_by_strain[i_distMat_by_strain,names(v.STRAIN_TO_SEQ)] = df.distMat_by_seq[v.STRAIN_TO_SEQ[i_distMat_by_strain], v.STRAIN_TO_SEQ]
    }
    
    ls.distMats_by_gene[[gene]] = df.distMat_by_strain
  }
  
  save(ls.distMats_by_gene
       , ls.refDist_by_gene
       , file='OBJECTS/list_of_strain_distMats_by_gene.R')
}else{
  load('OBJECTS/list_of_strain_distMats_by_gene.R')
}
