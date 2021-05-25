rm(list=ls())
# 05.25.2020
# CYZ & MM
# This is a master file for our pipeline.

# NUCLEOTIDE ANALYSIS
source('nucleotide.R')

# PROTEIN ANALYSIS
# (1) Process Raw Inputs --------------------------------------------------

source('process_inputs.R')
source('parse_metadata.R')
source('parse_inputs.R')


# (2) Analyze sequences ---------------------------------------------------

source('analyze_seqs.R')

# (3) Visualizations ------------------------------------------------------

source('visualize_distances.R')
source('visualize_data_overview.R')
