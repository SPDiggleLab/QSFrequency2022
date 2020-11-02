rm(list=ls())
require(ggplot2)
require(reshape2)
require(Biostrings)
require(cowplot)

# Load Data
load('OBJECTS/gene_trunc_stats.R')
load('OBJECTS/list_of_strain_distMats_by_gene.R')
load('OBJECTS/PCA_VIZ.R')

df.truncations$similarity_mean = NA
df.truncations$similarity_sd = NA

temp.truncDist = list()

for(gene in names(ls.refDist_by_gene)){
  # Get Reference
  ref_aminos = read.table(paste('INPUT/ref_genes/', gene, '_PAO1_protein.txt', sep='', collapse='')
                          , stringsAsFactors = F)$V1
  ref_aminos = AAString(ref_aminos) # read in as biocstring object
  
  # Calculate self-distance score
  self_align = pairwiseAlignment(ref_aminos
                                 , ref_aminos
                                 , substitutionMatrix = "BLOSUM80"
                                 , gapOpening=11
                                 , gapExtension=1)@score
  
  # normalize to self
  v.dist = ls.refDist_by_gene[[gene]] / self_align
  
  # get simple stats
  mean_dist = mean(v.dist)
  sd_dist = sd(v.dist)
  
  df.truncations[gene, 'similarity_mean'] = mean_dist
  df.truncations[gene, 'similarity_sd'] = sd_dist
  
  # for plotting histogram
  temp.truncDist[[gene]] = data.frame('gene' = gene, 'dist' = v.dist, 'same' = 0+(v.dist!=1))
}

temp.truncDist_df = do.call('rbind', temp.truncDist)
temp.truncDist_df$sys = sapply(temp.truncDist_df$gene, function(x){substr(x,1,3)})

# For Fig 2a
df.lasR_meta = ls.complete_gene_MetaSeq[['lasR']][,c('GENE', 'HOST', 'SOURCE', 'ENV', 'SEQUENCE')]
lasR_ref = read.table(paste('INPUT/ref_genes/lasR_PAO1_protein.txt', sep='', collapse='')
                      , stringsAsFactors = F)$V1

# Get AAs
b.truncs = sapply(df.lasR_meta$SEQUENCE, function(x){
  seqDNA = toString(translate(DNAString(x), if.fuzzy.codon = c('solve', 'X'), no.init.codon = T))
  return(gregexpr('\\*', seqDNA)[[1]][1] == nchar(lasR_ref))
})

# Add column: TRUNCATED
nBases_lasR = nchar(lasR_ref)*3
df.lasR_meta$TRUNCATED = ifelse(as.vector(b.truncs), 'FULL', 'TRUNC')#sapply(df.lasR_meta$SEQUENCE, function(x){nchar(x) < nBases_lasR})
df.lasR_meta$SEQUENCE = NULL # remove SEQ for easier reading


# If SOURCE is NA and HOST is HUMAN, assign source to non-CF
df.lasR_meta$SOURCE[(is.na(df.lasR_meta$SOURCE) & df.lasR_meta$HOST == 'Human')] = 'non-CF'

# Animals are part of the environment
df.lasR_meta$ENV[df.lasR_meta$HOST == 'Animal'] = TRUE

# Add a GROUP column for plotting
df.lasR_meta$GROUP = ifelse(df.lasR_meta$ENV, 'ENV', df.lasR_meta$SOURCE)

# Remove Non-CF
df.lasR_meta = df.lasR_meta[!df.lasR_meta$GROUP == 'non-CF',]


# PLOT --------------------------------------------------------------------
fig1a = ggplot(df.truncations, aes(x = gene, y = similarity_mean)) + 
  theme_gray(base_size = 14) +
  geom_bar(stat = 'identity', color = 'black') +
  geom_errorbar(size =2, aes(ymin = similarity_mean-similarity_sd
                    , ymax = similarity_mean+similarity_sd)
                , width = .5) +
  ylim(c(0, 1.2)) +
  xlab('') + 
  ylab('Normalized Similarity Score') +
  scale_fill_brewer(palette = 'Set1') +
  theme(legend.position = 'none'
        , axis.text.x = element_text(angle = 45, hjust=1, vjust=1))
  
fig1a2 = ggplot(df.truncations, aes(x = n_base, y = similarity_mean)) +
  theme_gray(base_size=14) + 
  geom_point(size=2) + 
  geom_text(aes(label = ifelse(similarity_mean<0.993, gene, ''))
            , hjust=-0.3, vjust=0) + 
  ylab('Normalized Similarity Score') +
  xlim(600, 1600)

fig1a3 = ggplot(df.truncations, aes(x = gene)) + 
  theme_gray(base_size = 14) +
  geom_bar(color = 'grey30', stat='identity', alpha = 0.4, fill = 'grey10', aes(y = n_unique_strains)) + 
  geom_bar(color = 'black', stat='identity', alpha = 0.9, aes(y = n_unique_seqs, fill = sapply(df.truncations$gene, function(x){substr(x,1,3)}))) + 
  scale_fill_brewer(palette = 'Set1') +
  theme(legend.title=element_blank()) + 
  ylab('n sequences') + 
  xlab('') +
  theme(legend.position = 'none')
fig1a3

fig1a4 = ggplot(temp.truncDist_df, aes(x = gene, y = 1-dist, color = sys)) + 
  theme_gray(base_size = 14) + 
  #geom_point(alpha=0.7) +
  geom_jitter(alpha = 0.7) + 
  scale_color_brewer(palette = 'Set1') +
  theme(legend.title=element_blank()) + 
  xlab('') + 
  ylab('dissimilarity') +
  facet_grid(.~gene, scale = 'free_x') +
  theme(legend.position = 'none')
fig1a4


fig1b = ls.pcas[['lasR']] +
  theme_gray(base_size = 14) + 
  geom_point(color = 'black', size = 3) + 
  theme(legend.position = 'none') + 
  ggtitle('')

fig1c = ls.pcas[['lasI']] +
  theme_gray(base_size = 14)  + 
  geom_point(color = 'black', size = 3) + 
  theme(legend.position = 'none') + 
  ggtitle('')

fig2a = ggplot(df.lasR_meta, aes(x = GROUP, fill = TRUNCATED)) + 
  geom_bar(stat = 'count', color = 'black') + 
  xlab('')
  

fig2b = ls.pcas[['lasR']]

fig1_full = plot_grid(fig1a3, fig1b, fig1a4, fig1c 
                      , labels = c('A', 'C', 'B', 'D')
                      , nrow = 2
                      , rel_widths = c(1.5, 1, 1.5, 1))
fig1_full
fig2_full = plot_grid(fig2a, fig2b, labels = c('A', 'B')
                      , rel_widths = c(1, 2))
fig2_full

ggsave('2020-11-02_fig1-full.png', fig1_full, device = 'png', width = 11, height = 8)
ggsave('2020-11-02_fig2-full.png', fig2_full, device = 'png', width = 9, height = 5)


