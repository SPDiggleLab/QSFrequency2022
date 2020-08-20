require(factoextra)
require(ape)

# PCA by environment
ls.pcas = list()
for(gene in names(ls.distMats_by_gene)){
  # [1] "lasI" "lasR" "pqsA" "pqsB" "pqsC" "pqsD" "pqsE" "pqsH" "pqsL" "pqsR" "rhlI" "rhlR"

  # Load Data
  df.dist = ls.distMats_by_gene[[gene]]
  meta = ls.complete_gene_MetaSeq[[gene]][,c('STRAIN', 'GENE', 'HOST', 'SOURCE', 'ENV', 'PI')]
  
  # Clean Data
  meta = meta[!duplicated(meta$STRAIN),] # remove duplicates
  rownames(meta) = meta$STRAIN
  meta = meta[rownames(df.dist),] # order w.r.t. df.dist
  
  meta$SOURCE[meta$ENV] = 'ENV'
  meta$SOURCE[meta$SOURCE == 'non-CF'] = NA # Remove non_CF (08.20.2020)
  b.has_SOURCE = !is.na(meta$SOURCE)
    
  df.dist = data.matrix(df.dist)
  df.dist = df.dist[b.has_SOURCE, b.has_SOURCE]

  # PCA on a reduced set.
  tryCatch(expr = {
    p_pca = fviz_pca_ind(princomp(df.dist), label = '', habillage = meta[rownames(df.dist), 'SOURCE']
                         , invisible="quali", pointsize = 2)
    ls.pcas[[gene]] = p_pca
    ggsave(filename = paste('./OUTPUT/PCAs/pca_', gene, '.png', sep='', collapse = '')
           , plot = p_pca
           , height = 5
           , width = 5
           , device = 'png')
  }
  , warning = function(cond) {
    p_pca = fviz_pca_ind(princomp(df.dist), label = '', habillage = meta[rownames(df.dist), 'SOURCE']
                         , invisible="quali", pointsize = 2)
    ls.pcas[[gene]] = p_pca
    ggsave(filename = paste('./OUTPUT/PCAs/pca_', gene, '.png', sep='', collapse = '')
           , plot = p_pca
           , height = 5
           , width = 5
           , device = 'png')
  }
  , error = function(cond) {
    df.dist_in = unique(t(unique(df.dist)))
    meta_in = meta[rownames(df.dist_in),]
    p_pca = fviz_pca_ind(princomp(df.dist_in), label = '', invisible="quali", pointsize = 2, habillage = meta_in$SOURCE)
    temp.pcas = list()
    temp.pcas[[gene]] = p_pca
    save(temp.pcas, file = paste('./temp/PCA_', gene, '.R', sep=''))
    ggsave(filename = paste('./OUTPUT/PCAs/pca_', gene, '.png', sep='', collapse = '')
           , plot = p_pca
           , height = 5
           , width = 5
           , device = 'png')
  })
}

for(gene_pca in names(ls.distMats_by_gene)[!names(ls.distMats_by_gene) %in% names(ls.pcas)]){
  load(paste('./temp/PCA_', gene_pca, '.R', sep=''))
  ls.pcas[[gene_pca]] = temp.pcas[[gene_pca]]
}


save(ls.pcas, file = 'OBJECTS/PCA_VIZ.R')

# Assume distances are euclidean and create a pcoa