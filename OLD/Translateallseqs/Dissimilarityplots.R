rm(list=ls())

# Conan Y Zhao - 2019.10.16
#   This code does protein alignment and gives protein alignment scores for LasI
#   and LasR genes obtained from the IPCD. Input data is from Kathleen O'Connor,
#   received 02.04.2020

require(Biostrings)
require(factoextra)
require(ape)
require(vegan)

# Notes


# (0) Helpful Functions ---------------------------------------------------
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
# (1) Inputs --------------------------------------------------------------

# Load Inputs
lasI_ref=read.table('inputs/lasI_PAO1_protein.txt', header = F, stringsAsFactors = F)[,1]
lasI_all=read.table('inputs/lasI.txt', header = F, stringsAsFactors = F)[,1]

lasR_ref=read.table('inputs/lasR_PAO1_protein.txt', header = F, stringsAsFactors = F)[,1]
lasR_all=read.table('inputs/lasR.txt', header = F, stringsAsFactors=F)[,1]

rhlR_ref=read.table('inputs/rhlR_PAO1_protein.txt', header = F, stringsAsFactors = F)[,1]
rhlR_all=read.table('inputs/rhlR.txt', header = F, stringsAsFactors=F)[,1]

rhlI_ref=read.table('inputs/rhlI_PAK_protein.txt', header = F, stringsAsFactors = F)[,1]
rhlI_all=read.table('inputs/rhlI.txt', header = F, stringsAsFactors=F)[,1]

rhlI_ref=read.table('inputs/rhlI_PAK_protein.txt', header = F, stringsAsFactors = F)[,1]
rhlI_all=read.table('inputs/rhlI.txt', header = F, stringsAsFactors=F)[,1]

pqsA_ref=read.table('inputs/pqsA_PAO1_protein.txt', header = F, stringsAsFactors = F)[,1]
pqsA_all=read.table('inputs/pqsA.txt', header = F, stringsAsFactors=F)[,1]

pqsB_ref=read.table('inputs/pqsB_PAO1_protein.txt', header = F, stringsAsFactors = F)[,1]
pqsB_all=read.table('inputs/pqsB.txt', header = F, stringsAsFactors=F)[,1]

pqsC_ref=read.table('inputs/pqsC_PAO1_protein.txt', header = F, stringsAsFactors = F)[,1]
pqsC_all=read.table('inputs/pqsC.txt', header = F, stringsAsFactors=F)[,1]

pqsD_ref=read.table('inputs/pqsD_PAK_protein.txt', header = F, stringsAsFactors = F)[,1]
pqsD_all=read.table('inputs/pqsD.txt', header = F, stringsAsFactors=F)[,1]

pqsE_ref=read.table('inputs/pqsE_PAO1_protein.txt', header = F, stringsAsFactors = F)[,1]
pqsE_all=read.table('inputs/pqsE.txt', header = F, stringsAsFactors=F)[,1]

pqsR_ref=read.table('inputs/pqsR_PAO1_protein.txt', header = F, stringsAsFactors = F)[,1]
pqsR_all=read.table('inputs/pqsR.txt', header = F, stringsAsFactors=F)[,1]

pqsH_ref=read.table('inputs/pqsH_PAO1_protein.txt', header = F, stringsAsFactors = F)[,1]
pqsH_all=read.table('inputs/pqsH.txt', header = F, stringsAsFactors=F)[,1]

pqsL_ref=read.table('inputs/pqsL_PAO1_protein.txt', header = F, stringsAsFactors = F)[,1]
pqsL_all=read.table('inputs/pqsL.txt', header = F, stringsAsFactors=F)[,1]

# Cleanup
lasI_all = CLEAN_AA_SEQS(lasI_all)
lasR_all = CLEAN_AA_SEQS(lasR_all)
rhlR_all = CLEAN_AA_SEQS(rhlR_all)
rhlI_all = CLEAN_AA_SEQS(rhlI_all)
pqsA_all = CLEAN_AA_SEQS(pqsA_all)
pqsB_all = CLEAN_AA_SEQS(pqsB_all)
pqsC_all = CLEAN_AA_SEQS(pqsC_all)
pqsD_all = CLEAN_AA_SEQS(pqsD_all)
pqsE_all = CLEAN_AA_SEQS(pqsE_all)
pqsR_all = CLEAN_AA_SEQS(pqsR_all)
pqsH_all = CLEAN_AA_SEQS(pqsH_all)
pqsL_all = CLEAN_AA_SEQS(pqsL_all)

# Ensure first entry is the reference 
lasI_all = c(lasI_ref, lasI_all)
lasR_all = c(lasR_ref, lasR_all)
rhlR_all = c(rhlR_ref, rhlR_all)
rhlI_all = c(rhlI_ref, rhlI_all)
pqsA_all = c(pqsA_ref, pqsA_all)
pqsB_all = c(pqsB_ref, pqsB_all)
pqsC_all = c(pqsC_ref, pqsC_all)
pqsD_all = c(pqsD_ref, pqsD_all)
pqsR_all = c(pqsR_ref, pqsR_all)
pqsH_all = c(pqsH_ref, pqsH_all)
pqsL_all = c(pqsL_ref, pqsL_all)

# Uniques
lasI_unique = unique(lasI_all)
lasR_unique = unique(lasR_all)
rhlR_unique = unique(rhlR_all)
rhlI_unique = unique(rhlI_all)
pqsA_unique = unique(pqsA_all)
pqsB_unique = unique(pqsB_all)
pqsC_unique = unique(pqsC_all)
pqsD_unique = unique(pqsD_all)
pqsR_unique = unique(pqsR_all)
pqsH_unique = unique(pqsH_all)
pqsL_unique = unique(pqsL_all)

# (2) Analyses ------------------------------------------------------------

if(!file.exists('Full_dfRefPlot_Scores.R')){
  # Step 1: Align against a reference genome
  v.lasI_ref_scores=ALIGN_VS_REF(lasI_ref, lasI_all)
  v.lasR_ref_scores=ALIGN_VS_REF(lasR_ref, lasR_all)
  v.rhlR_ref_scores=ALIGN_VS_REF(rhlR_ref, rhlR_all)
  v.rhlI_ref_scores=ALIGN_VS_REF(rhlI_ref, rhlI_all)
  v.pqsA_ref_scores=ALIGN_VS_REF(pqsA_ref, pqsA_all)
  v.pqsB_ref_scores=ALIGN_VS_REF(pqsB_ref, pqsB_all)
  v.pqsC_ref_scores=ALIGN_VS_REF(pqsC_ref, pqsC_all)
  v.pqsD_ref_scores=ALIGN_VS_REF(pqsD_ref, pqsD_all)
  v.pqsE_ref_scores=ALIGN_VS_REF(pqsE_ref, pqsE_all)
  v.pqsR_ref_scores=ALIGN_VS_REF(pqsR_ref, pqsR_all)
  v.pqsH_ref_scores=ALIGN_VS_REF(pqsH_ref, pqsH_all)
  v.pqsL_ref_scores=ALIGN_VS_REF(pqsL_ref, pqsL_all)
  
  # 1b: normalize
  v.lasI_ref_scores = v.lasI_ref_scores/v.lasI_ref_scores[1]
  v.lasR_ref_scores = v.lasR_ref_scores/v.lasR_ref_scores[1]
  v.rhlR_ref_scores = v.rhlR_ref_scores/v.rhlR_ref_scores[1]
  v.rhlI_ref_scores = v.rhlI_ref_scores/v.rhlI_ref_scores[1]
  v.pqsA_ref_scores = v.pqsA_ref_scores/v.pqsA_ref_scores[1]
  v.pqsB_ref_scores = v.pqsB_ref_scores/v.pqsB_ref_scores[1]
  v.pqsC_ref_scores = v.pqsC_ref_scores/v.pqsC_ref_scores[1]
  v.pqsD_ref_scores = v.pqsD_ref_scores/v.pqsD_ref_scores[1]
  v.pqsE_ref_scores = v.pqsE_ref_scores/v.pqsE_ref_scores[1]
  v.pqsR_ref_scores = v.pqsR_ref_scores/v.pqsR_ref_scores[1]
  v.pqsH_ref_scores = v.pqsH_ref_scores/v.pqsH_ref_scores[1]
  v.pqsL_ref_scores = v.pqsL_ref_scores/v.pqsL_ref_scores[1]
  
  # 1c: plot
  df.refplots = rbind(data.frame(gene = 'lasR', score = v.lasR_ref_scores, seq = names(v.lasR_ref_scores))
                      , data.frame(gene = 'lasI', score = v.lasI_ref_scores, seq = names(v.lasI_ref_scores))
                      , data.frame(gene = 'rhlR', score = v.rhlR_ref_scores, seq = names(v.rhlR_ref_scores))
                      , data.frame(gene = 'rhlI', score = v.rhlI_ref_scores, seq = names(v.rhlI_ref_scores))
                      , data.frame(gene = 'pqsA', score = v.pqsA_ref_scores, seq = names(v.pqsA_ref_scores))
                      , data.frame(gene = 'pqsB', score = v.pqsB_ref_scores, seq = names(v.pqsB_ref_scores))
                      , data.frame(gene = 'pqsC', score = v.pqsC_ref_scores, seq = names(v.pqsC_ref_scores))
                      , data.frame(gene = 'pqsD', score = v.pqsD_ref_scores, seq = names(v.pqsD_ref_scores))
                      , data.frame(gene = 'pqsE', score = v.pqsE_ref_scores, seq = names(v.pqsE_ref_scores))
                      , data.frame(gene = 'pqsH', score = v.pqsH_ref_scores, seq = names(v.pqsH_ref_scores))
                      , data.frame(gene = 'pqsL', score = v.pqsL_ref_scores, seq = names(v.pqsL_ref_scores))
                      , data.frame(gene = 'pqsR', score = v.pqsR_ref_scores, seq = names(v.pqsR_ref_scores)))
  save(df.refplots, file = 'Full_dfRefPlot_Scores.R')
}else{
  print('WARNING: Loading pre-calculated scores')
  load('Full_dfRefPlot_Scores.R')
}

# Calculate means
temp.genes = as.character(unique(df.refplots$gene))
temp.out = data.frame('GENE' = NA, 'mean'=NA, 'N' = NA, 'sd' = NA)
for(gene_in in temp.genes){
  temp.data = df.refplots$score[df.refplots$gene == gene_in]
  temp.mean = mean(temp.data)
  temp.N = length(temp.data)
  temp.sd = sd(temp.data)
  temp.df = data.frame('GENE' = gene_in, 'mean'=temp.mean, 'N' = temp.N, 'sd' = temp.sd)
  temp.out = rbind(temp.out, temp.df)
}
temp.out = temp.out[-1,]
write.csv(temp.out, 'SUMMARY_STATS_PER_GENE.csv')

ggplot(df.refplots, aes(x = gene, y= 1-score, color = gene)) + 
  geom_boxplot() +
  ylab("Dissimilarity") #+
  #scale_y_log10()

# Step 2: Align against each other, and create a PCoA
if(!file.exists('geneFullDists.R')){
  df.lasI_fullDist = ALIGN_VS_SELF(lasI_all)
  df.rhlR_fullDist = ALIGN_VS_SELF(rhlR_all)
  df.lasR_fullDist = ALIGN_VS_SELF(lasR_all)
  df.rhlI_fullDist = ALIGN_VS_SELF(rhlI_all)
  save(df.lasI_fullDist, df.rhlR_fullDist, df.rhlI_fullDist, df.lasR_fullDist, file = 'geneFullDists.R')
}else{
  print('WARNING: Loading precalculated Distances')
  load('geneFullDists.R')
}
# 
# # Step 3: Ref distances by seq position
# data("BLOSUM80")
# 
# df.refPos = 3


# PCOA of lasI
fviz_pca_ind(prcomp(df.lasI_fullDist), geom.ind = 'point') + ggtitle('lasI Proteins PCA')
# PCoA of LasR
fviz_pca_ind(prcomp(df.lasR_fullDist), geom.ind = 'point') + ggtitle('lasR Proteins PCA')

fviz_pca_ind(prcomp(df.rhlR_fullDist), geom.ind = 'point') + ggtitle('rhlR Proteins PCA')

fviz_pca_ind(prcomp(df.rhlI_fullDist), geom.ind = 'point') + ggtitle('rhlI Proteins PCA')
