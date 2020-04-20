rm(list=ls())

# CYZ
#   This script creates PCAs to show LasR variance between lab strains of 
#   PA from the IPCD database. 

# (0) Libraries -----------------------------------------------------------
require(Biostrings)
require(factoextra)
require(ape)
require(vegan)

# (1) Useful Functions ----------------------------------------------------
ALIGN_VS_REF = function(REF=lasI_ref, ALL=lasI_all, gapOpen=11, gapExt=1){
  sapply(ALL, function(x){
    return (pairwiseAlignment(REF
                              , x
                              , substitutionMatrix = "BLOSUM80"
                              , gapOpening=gapOpen
                              , gapExtension=gapExt)@score)
  })
}

CLEAN_AA_SEQS = function(ALL){
  ret = sapply(ALL, function(x){gsub('-', '', x)})
  ret = sapply(ret, function(x){gsub('\\?', 'X', x)})
  return(ret)
}

ALIGN_VS_SELF = function(ALL=lasI_all){
  nSeq = length(ALL)
  
  # Calculate distances
  df.uniqueAADist = sapply(unique(ALL), function(x){
    ALIGN_VS_REF(x, unique(ALL))
  })
  
  # Format and Fill full matrix
  df.fullAADist = data.frame(matrix(rep(0, nSeq^2), ncol = nSeq))
  colnames(df.fullAADist) = paste('s', 1:nSeq, sep='')
  rownames(df.fullAADist) = paste('s', 1:nSeq, sep='')
  
  for(i in 1:length(ALL)){
    for(j in 1:length(ALL)){
      df.fullAADist[i,j] = df.uniqueAADist[ALL[i], ALL[j]]
    }
  }
  
  return(df.fullAADist)
}

# (2) Analyses ------------------------------------------------------------

# Load Data
v.inputFileNames = list.files('Lab_compare/')

ls.in = lapply(v.inputFileNames, function(x){
  df.in = read.table(paste('Lab_compare/', x, sep='', collapse=''), header = T, sep = '\t', stringsAsFactors = F)
  df.in = data.frame(df.in, stringsAsFactors = F)
  df.in$lab = strsplit(x, '_')[[1]][1]
  
  return(df.in)
})


df.in = do.call('rbind', ls.in)

# Align vs self
v.protSeqs = CLEAN_AA_SEQS(df.in$Protein)
temp.dist = ALIGN_VS_SELF(v.protSeqs)

fviz_pca_ind(prcomp(temp.dist), geom.ind = 'point', habillage = df.in$lab) + ggtitle('Variance Among Labs')
