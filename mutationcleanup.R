rm(list=ls())
library(tidyverse)
library(ggvenn)
library(ggplot2)
library(gridExtra)


df <- read.csv2("gammaalldata.csv", sep = ",", header = TRUE)

df$Name <- str_extract_all(df$Contig, "[A-Z]+")

df1 <- df[!grepl("Native", df$Match_Type),]
df2 <- df1[!grepl("Contig Edge", df1$Match_Type),]
df3 <- df2[!grepl("contig edge", df2$Match_Type),]


#df3 <- apply(df3,2,as.character) #need to run this line to make csv, however re-run previous line to turn df3 back into df for future coding
#write.csv(df3, "mutants.csv")

#Create lists of isolates with mutations in each genes and save as a list to use in Venndiagrams            
lasRdf <- df3[grep("PAO1_lasR", df3$Gene), ]
listlasR <- as.list(lasRdf$Name)
print(listlasR)

lasIdf <- df3[grep("PAO1_lasI", df3$Gene), ]
listlasI <- as.list(lasIdf$Name)
print(listlasI)

rhlRdf <- df3[grep("PAO1_rhlR", df3$Gene), ]
listrhlR <- as.list(rhlRdf$Name)
print(listrhlR)

rhlIdf <- df3[grep("PAK_rhlI", df3$Gene), ]
listrhlI <- as.list(rhlIdf$Name)
print(listrhlI)

pqsAdf <- df3[grep("PAO1_pqsA", df3$Gene), ]
listpqsA <- as.list(pqsAdf$Name)
print(listpqsA)

pqsBdf <- df3[grep("PAO1_pqsB", df3$Gene), ]
listpqsB <- as.list(pqsBdf$Name)
print(listpqsB)

pqsCdf <- df3[grep("PAO1_pqsC", df3$Gene), ]
listpqsC <- as.list(pqsCdf$Name)
print(listpqsC)

pqsDdf <- df3[grep("PAK_pqsD", df3$Gene), ]
listpqsD <- as.list(pqsDdf$Name)
print(listpqsD)

pqsEdf <- df3[grep("PAO1_pqsE", df3$Gene), ]
listpqsE <- as.list(pqsEdf$Name)
print(listpqsE)

pqsHdf <- df3[grep("PAO1_pqsH", df3$Gene), ]
listpqsH <- as.list(pqsHdf$Name)
print(listpqsH)

pqsLdf <- df3[grep("PAO1_pqsL", df3$Gene), ]
listpqsL <- as.list(pqsLdf$Name)
print(listpqsL)

pqsRdf <- df3[grep("PAO1_pqsR", df3$Gene), ]
listpqsR <- as.list(pqsRdf$Name)
print(listpqsR)

#PAO1_vqsR
vqsRdf <- df3[grep("PAO1_vqsR", df3$Gene), ]
listvqsR <- as.list(vqsRdf$Name)
print(listvqsR)

#PAO1_rsaL
rsaLdf <- df3[grep("PAO1_rsaL", df3$Gene), ]
listrsaL <- as.list(rsaLdf$Name)
print(listrsaL)

#PAO1_psdR
psdRdf <- df3[grep("PAO1_psdR", df3$Gene), ]
listpsdR <- as.list(psdRdf$Name)
print(listpsdR)

#PAO1_qscR
qscRdf <- df3[grep("PAO1_qscR", df3$Gene), ]
listqscR <- as.list(qscRdf$Name)
print(listqscR)

#mexT
MPAO1mexTdf <- df3[grep("MPAO1_mexT", df3$Gene), ]
listMPAO1mexT <- as.list(MPAO1mexTdf$Name)
print(listMPAO1mexT)

mexTdf <- df3[grep("U0330A_mexT", df3$Gene), ]
listmexT <- as.list(mexTdf$Name)
print(listmexT)

#input gene lists
set.seed(20190708)
#here is where you make venn diagrams - just enter in what you want to compare, these are examples
x.mexT <- list(
  LasR = sample(listlasR),
  MexT = sample(listmexT))

x.psdR <- list(
  LasR = sample(listlasR),
  PsdR = sample(listpsdR))

x.lasI <- list(
  LasR = sample(listlasR),
  LasI = sample(listlasI))

x.rhlR <- list(
  LasR = sample(listlasR),
  RhlR = sample(listrhlR))

x.rhlI <- list(
  LasR = sample(listlasR),
  RhlI = sample(listrhlI))

x.pqsA <- list(
  LasR = sample(listlasR),
  PqsA = sample(listpqsA))

x.pqsB <- list(
  LasR = sample(listlasR),
  PqsB = sample(listpqsB))

x.pqsC <- list(
  LasR = sample(listlasR),
  PqsC = sample(listpqsC))

x.pqsD <- list(
  LasR = sample(listlasR),
  PqsD = sample(listpqsD))

x.pqsE <- list(
  LasR = sample(listlasR),
  PqsE = sample(listpqsE))

x.pqsL <- list(
  LasR = sample(listlasR),
  PqsL = sample(listpqsL))

x.pqsH <- list(
  LasR = sample(listlasR),
  PqsH = sample(listpqsH))

x.pqsR <- list(
  LasR = sample(listlasR),
  PqsR = sample(listpqsR))

x.rsaL <- list(
  LasR = sample(listlasR),
  RsaL = sample(listrsaL))

x.vqsR <- list(
  LasR = sample(listlasR),
  VqsR = sample(listvqsR))

x.qscR <- list(
  LasR = sample(listlasR),
  QscR = sample(listqscR))

#create Venn diagram with nice colors and labels
p.mexT <- ggvenn(
  x.mexT,
  show_percentage = FALSE,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

p.psdR <- ggvenn(
  x.psdR,
  show_percentage = FALSE,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

p.lasI <- ggvenn(
  x.lasI,
  show_percentage = FALSE,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

p.rhlR <- ggvenn(
  x.rhlR,
  show_percentage = FALSE,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

p.rhlI <- ggvenn(
  x.rhlI,
  show_percentage = FALSE,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

p.pqsA <- ggvenn(
  x.pqsA,
  show_percentage = FALSE,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

p.pqsB <- ggvenn(
  x.pqsB,
  show_percentage = FALSE,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

p.pqsC <- ggvenn(
  x.pqsC,
  show_percentage = FALSE,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

p.pqsD <- ggvenn(
  x.pqsD,
  show_percentage = FALSE,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

p.pqsE <- ggvenn(
  x.pqsE,
  show_percentage = FALSE,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

p.pqsL <- ggvenn(
  x.pqsL,
  show_percentage = FALSE,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

p.pqsH <- ggvenn(
  x.pqsH,
  show_percentage = FALSE,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

p.pqsR <- ggvenn(
  x.pqsR,
  show_percentage = FALSE,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

p.rsaL <- ggvenn(
  x.rsaL,
  show_percentage = FALSE,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

p.vqsR <- ggvenn(
  x.vqsR,
  show_percentage = FALSE,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

p.qscR <- ggvenn(
  x.qscR,
  show_percentage = FALSE,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)


#create multi panel figure with all of our graphs
plot <- grid.arrange(p.mexT, p.psdR, p.lasI, p.rhlR, p.rhlI, p.pqsA, p.pqsB, p.pqsC,
             p.pqsD, p.pqsE, p.pqsL, p.pqsH, p.pqsR, p.rsaL, p.vqsR, p.qscR, nrow = 4)

#save plot as figure
ggsave('proteinvenn-full.png', plot, device = 'png')

