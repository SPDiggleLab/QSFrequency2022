rm(list=ls())
# 04.19.2020
# CYZ
# This is a master file for our pipeline.




# (1) Process Raw Inputs --------------------------------------------------

source('process_inputs.R')
source('parse_metadata.R')
source('parse_inputs.R')


# (2) Analyze sequences ---------------------------------------------------

source('analyze_seqs.R')

# (3) Visualizations

source('visualize_distances.R')