require(factoextra)
require(ape)

# PCA by environment
for(gene in names(ls.distMats_by_gene)){
  # [1] "lasI" "lasR" "pqsA" "pqsB" "pqsC" "pqsD" "pqsE" "pqsH" "pqsL" "pqsR" "rhlI" "rhlR"

  # Load Data
  df.dist = ls.distMats_by_gene[[gene]]
  meta = ls.complete_gene_MetaSeq[[gene]][,c('STRAIN', 'GENE', 'HOST', 'SOURCE', 'ENV', 'PI')]
  
  # Clean Data
  meta = meta[!duplicated(meta$STRAIN),] # remove duplicates
  rownames(meta) = meta$STRAIN
  meta = meta[rownames(df.dist),] # order w.r.t. df.dist
  
  b.has_SOURCE = !is.na(meta$SOURCE)
  
  df.dist = data.matrix(df.dist)
  df.dist = df.dist[b.has_SOURCE, b.has_SOURCE]

  # PCA on a reduced set.
  tryCatch(expr = {
    p_pca = fviz_pca_ind(princomp(df.dist), label = '', habillage = meta[rownames(df.dist), 'SOURCE'], invisible="quali", pointsize = 2)
    ggsave(filename = paste('./OUTPUT/PCAs/pca_', gene, '.png', sep='', collapse = '')
           , plot = p_pca
           , device = 'png')
  }
  , warning = function(cond) {
    p_pca = fviz_pca_ind(princomp(df.dist), label = '', habillage = meta[rownames(df.dist), 'SOURCE'], invisible="quali", pointsize = 2)
    ggsave(filename = paste('./OUTPUT/PCAs/pca_', gene, '.png', sep='', collapse = '')
           , plot = p_pca
           , device = 'png')
  }
  , error = function(cond) {
    df.dist_in = unique(t(unique(df.dist)))
    meta_in = meta[rownames(meta),]
    p_pca = fviz_pca_ind(princomp(df.dist_in), label = '', invisible="quali", pointsize = 2)
    ggsave(filename = paste('./OUTPUT/PCAs/pca_', gene, '.png', sep='', collapse = '')
           , plot = p_pca
           , device = 'png')
  })
}




# Assume distances are euclidean and create a pcoa