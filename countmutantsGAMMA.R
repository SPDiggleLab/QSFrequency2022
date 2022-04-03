rm(list=ls())
library(tidyverse)

df <- read.csv2("gammaalldata.csv", sep = ",", header = TRUE)

#make it easier to group strains by renaming them to remove contig info
df$Name <- str_extract_all(df$Contig, "[A-Z]+")

#remove all sequences that are at the contig edge
df <- df[!grepl("Contig Edge", df$Match_Type),]
df <- df[!grepl("contig edge", df$Match_Type),]

#(i)create dataframes for each gene (ii) count how many sequences are available for each gene -
#expect somewhere in the 800s (iii) remove all Native proteins to create mutant list
#(iv) count the number of mutants there are
lasR <- df[grep("PAO1_lasR", df$Gene), ]
nrow(lasR)
lasRmutant <- lasR[!grepl("Native", lasR$Match_Type),]
nrow(lasRmutant)

lasI <- df[grep("PAO1_lasI", df$Gene), ]
nrow(lasI)
lasImutant <- lasI[!grepl("Native", lasI$Match_Type),]
nrow(lasImutant)

rhlR <- df[grep("PAO1_rhlR", df$Gene), ]
nrow(rhlR)
rhlRmutant <- rhlR[!grepl("Native", rhlR$Match_Type),]
nrow(rhlRmutant)

rhlI <- df[grep("PAO1_rhlI", df$Gene), ]
nrow(rhlI)
rhlImutant <- rhlI[!grepl("Native", rhlI$Match_Type),]
nrow(rhlImutant)

rhlI_PAK <- df[grep("PAK_rhlI", df$Gene), ]
nrow(rhlI_PAK)
rhlImutant_PAK <- rhlI_PAK[!grepl("Native", rhlI_PAK$Match_Type),]
nrow(rhlImutant_PAK)

pqsA <- df[grep("PAO1_pqsA", df$Gene), ]
nrow(pqsA)
pqsAmutant <- pqsA[!grepl("Native", pqsA$Match_Type),]
nrow(pqsAmutant)

pqsB <- df[grep("PAO1_pqsB", df$Gene), ]
nrow(pqsB)
pqsBmutant <- pqsB[!grepl("Native", pqsB$Match_Type),]
nrow(pqsBmutant)

pqsC <- df[grep("PAO1_pqsC", df$Gene), ]
nrow(pqsC)
pqsCmutant <- pqsC[!grepl("Native", pqsC$Match_Type),]
nrow(pqsCmutant)

pqsD_PAK <- df[grep("PAK_pqsD", df$Gene), ]
nrow(pqsD_PAK)
pqsDmutant_PAK <- pqsD_PAK[!grepl("Native", pqsD_PAK$Match_Type),]
nrow(pqsDmutant_PAK)

pqsD <- df[grep("PAO1_pqsD", df$Gene), ]
nrow(pqsD)
pqsDmutant <- pqsD[!grepl("Native", pqsD$Match_Type),]
nrow(pqsDmutant)

pqsE <- df[grep("PAO1_pqsE", df$Gene), ]
nrow(pqsE)
pqsEmutant <- pqsE[!grepl("Native", pqsE$Match_Type),]
nrow(pqsEmutant)

pqsH <- df[grep("PAO1_pqsH", df$Gene), ]
nrow(pqsH)
pqsHmutant <- pqsH[!grepl("Native", pqsH$Match_Type),]
nrow(pqsHmutant)

pqsL <- df[grep("PAO1_pqsL", df$Gene), ]
nrow(pqsL)
pqsLmutant <- pqsL[!grepl("Native", pqsL$Match_Type),]
nrow(pqsLmutant)

pqsL_PAK <- df[grep("PAK_pqsL", df$Gene), ]
nrow(pqsL_PAK)
pqsLPAKmutant <- pqsL_PAK[!grepl("Native", pqsL_PAK$Match_Type),]
nrow(pqsLPAKmutant)

pqsR <- df[grep("PAO1_pqsR", df$Gene), ]
nrow(pqsR)
pqsRmutant <- pqsR[!grepl("Native", pqsR$Match_Type),]
nrow(pqsRmutant)

qscR <- df[grep("PAO1_qscR", df$Gene), ]
nrow(qscR)
qscRmutant <- qscR[!grepl("Native", qscR$Match_Type),]
nrow(qscRmutant)

rsaL <- df[grep("PAO1_rsaL", df$Gene), ]
nrow(rsaL)
rsaLmutant <- rsaL[!grepl("Native", rsaL$Match_Type),]
nrow(rsaLmutant)

vqsR <- df[grep("PAO1_vqsR", df$Gene), ]
nrow(vqsR)
vqsRmutant <- vqsR[!grepl("Native", vqsR$Match_Type),]
nrow(vqsRmutant)

mexT_MPAO1 <- df[grep("MPAO1_mexT", df$Gene), ]
nrow(mexT_MPAO1)
mexTMPAO1mutant <- mexT_MPAO1[!grepl("Native", mexT_MPAO1$Match_Type),]
nrow(mexTMPAO1mutant)

mexT <- df[grep("U0330A_mexT", df$Gene), ]
nrow(mexT)
mexTmutant <- mexT[!grepl("Native", mexT$Match_Type),]
nrow(mexTmutant)

psdR <- df[grep("PAO1_psdR", df$Gene), ]
nrow(psdR)
psdRmutant <- psdR[!grepl("Native", psdR$Match_Type),]
nrow(psdRmutant)

rpsL <- df[grep("PAO1_rpsL", df$Gene), ]
nrow(rpsL)
rpsLmutant <- rpsL[!grepl("Native", rpsL$Match_Type),]
nrow(rpsLmutant)

mucA <- df[grep("PAO1_mucA", df$Gene), ]
nrow(mucA)
mucAmutant <- mucA[!grepl("Native", mucA$Match_Type),]
nrow(mucAmutant)

#now we have the % mutation for each gene great! Now want some way to compare pqsL and
#lasR because we know they're obviously different
#do this by looking at the average Codon_Percent column  - want average of this
#the higher this number, the more similar

#lasR
listlasR <- lasR$Codon_Percent
listlasR <- as.numeric(as.character(listlasR))
mean(listlasR)

#lasI
listlasI <- lasI$Codon_Percent
listlasI <- as.numeric(as.character(listlasI))
mean(listlasI)

#rhlR
listrhlR <- rhlR$Codon_Percent
listrhlR <- as.numeric(as.character(listrhlR))
mean(listrhlR)

#rhlI
listrhlI <- rhlI$Codon_Percent
listrhlI <- as.numeric(as.character(listrhlI))
mean(listrhlI)

#rhlI_PAK
listrhlI_PAK <- rhlI_PAK$Codon_Percent
listrhlI_PAK <- as.numeric(as.character(listrhlI_PAK))
mean(listrhlI_PAK)

#mexT
listmexT <- mexT$Codon_Percent
listmexT <- as.numeric(as.character(listmexT))
mean(listmexT)

#mexT_MPAO1
listmexTMPAO1 <- mexT_MPAO1$Codon_Percent
listmexTMPAO1 <- as.numeric(as.character(listmexTMPAO1))
mean(listmexTMPAO1)

#pqsA
listpqsA <- pqsA$Codon_Percent
listpqsA <- as.numeric(as.character(listpqsA))
mean(listpqsA)

#pqsB
listpqsB <- pqsB$Codon_Percent
listpqsB <- as.numeric(as.character(listpqsB))
mean(listpqsB)

#pqsC
listpqsC <- pqsC$Codon_Percent
listpqsC <- as.numeric(as.character(listpqsC))
mean(listpqsC)

#pqsD
listpqsD_PAK <- pqsD_PAK$Codon_Percent
listpqsD_PAK <- as.numeric(as.character(listpqsD_PAK))
mean(listpqsD_PAK)

listpqsD <- pqsD$Codon_Percent
listpqsD <- as.numeric(as.character(listpqsD))
mean(listpqsD)

#pqsE
listpqsE <- pqsE$Codon_Percent
listpqsE <- as.numeric(as.character(listpqsE))
mean(listpqsE)

#pqsL
listpqsL <- pqsL$Codon_Percent
listpqsL <- as.numeric(as.character(listpqsL))
mean(listpqsL)

#pqsH
listpqsH <- pqsH$Codon_Percent
listpqsH <- as.numeric(as.character(listpqsH))
mean(listpqsH)

#pqsR
listpqsR <- pqsR$Codon_Percent
listpqsR <- as.numeric(as.character(listpqsR))
mean(listpqsR)

#qscR
listqscR <- qscR$Codon_Percent
listqscR <- as.numeric(as.character(listqscR))
mean(listqscR)

#psdR
listpsdR <- psdR$Codon_Percent
listpsdR <- as.numeric(as.character(listpsdR))
mean(listpsdR)

#rsaL
listrsaL <- rsaL$Codon_Percent
listrsaL <- as.numeric(as.character(listrsaL))
mean(listrsaL)

#vqsR
listvqsR <- vqsR$Codon_Percent
listvqsR <- as.numeric(as.character(listvqsR))
mean(listvqsR)

#rpsL
listrpsL <- rpsL$Codon_Percent
listrpsL <- as.numeric(as.character(listrpsL))
mean(listrpsL)

#mucA
listmucA <- mucA$Codon_Percent
listmucA <- as.numeric(as.character(listmucA))
mean(listmucA)