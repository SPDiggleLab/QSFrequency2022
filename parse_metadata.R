rm(list=ls())
## 04.11.2020 - CYZ
# This is a custom script to generate a standardized metadata output by lab from STRAIN_METADATA.xlsx
# STRAIN_METADATA.xlsx is a handmade data table copied and pasted from the IPCD web browser

require(readxl)


# (0) Load STRAIN_METADATA.xlsx -------------------------------------------

temp.in = read_xlsx('INPUT/RAW_INPUTS/STRAIN_METADATA.xlsx')
temp.in = data.frame(temp.in, stringsAsFactors = F)

# Standardize Human Pathogen Codes
v.humanPathoCodes = c('Spine pressure sore' = 'non-CF'
                    , 'Wound' = 'WND'
                    , 'Ulcer' = 'WND'
                    , 'Burn' = 'WND'
                    , 'Cystic fibrosis' = 'CF'
                    , 'Bronchiectasis' = 'non-CF'
                    , 'Pneumonia' = 'non-CF'
                    , 'Keratitis' = 'non-CF'
                    , 'Bacteraemia' = 'non-CF'
                    , 'Ear infection' = 'non-CF'
                    , 'UTI' = 'non-CF'
                    , 'Non-CF' = 'non-CF'
                    , 'Leg ulcer' = 'non-CF'
                    , 'Clinical non-CF' = 'non-CF'
                    , 'Screen' = 'non-CF'
                    , 'Hyponeutremia' = 'non-CF'
                    , 'Primary Ciliary Dyskinesia' = 'non-CF'
                    , 'Intestinal Cancer' = 'non-CF'
                    , 'Acute infection' = 'non-CF'
                    , 'COPD' = 'non-CF'
                    , 'Cystis fibrosis' = 'CF'
                    , 'Cystic Fibrosis' = 'CF')

# Standardize country Codes
v.countryCodes = unique(temp.in$Isolate.Country)
names(v.countryCodes) = v.countryCodes
temp.countryCodes = c('United Kingdom' = 'UK'
                   , 'United States' = 'USA'
                   , 'Canada' = 'CAN')
v.countryCodes[names(temp.countryCodes)] = temp.countryCodes


# (1) Create Standardized df.meta -----------------------------------------

# Create New Metadata
df.meta = data.frame('STRAIN_ID' = temp.in$Original.ID
                     , 'ALIAS_1' = do.call('c', lapply(temp.in$Aliases, function(x){strsplit(x, ', ')[[1]][1]}))
                     , 'ALIAS_2' = do.call('c', lapply(temp.in$Aliases, function(x){strsplit(x, ', ')[[1]][2]}))
                     , 'COUNTRY' = v.countryCodes[temp.in$Isolate.Country]
                     , 'HOST' = temp.in$Host
                     , 'SOURCE' = v.humanPathoCodes[temp.in$Human.Pathology]
                     , 'ENV' = (temp.in$Environmental == 'Yes')
                     , 'PI' = temp.in$Researcher.Name
                     , 'INCLUDE' = T)


# (2) Save ----------------------------------------------------------------

write.table(df.meta, 'INPUT/strain_meta.tsv', sep='\t', row.names=F)


# (3) Cleanup -------------------------------------------------------------

rm(df.meta, temp.in, temp.countryCodes, v.countryCodes, v.humanPathoCodes)
