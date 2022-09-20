library(tidyverse)
library(drc)
library(gtools)

standard <-read_csv("Proteasecombopaper/proteasestandard20min_4.13.22.csv") #our Trypsin standard
protease_data <- read_csv("Proteasecombopaper/proteasereadingscombo.csv") #our data


df <- standard
df1 <- protease_data


#use mod2
mod2 = lm(Trypsin ~ Fluorescence, data=df)
summary(mod2)
plot2 <- plot(df$`Fluorescence`, (df$Trypsin))
abline(mod2)

#process data using model for mod2
Protease <- predict(mod2, data.frame(Fluorescence=df1$Fluorescence), se.fit=FALSE)
Protease <- as.data.frame(Protease)
Protease <- cbind(Protease,df1$Strain,df1$Pa_OD)

#now nromalize protease readings by OD
Protease$Normalized_Protease <- Protease$Protease/Protease$`df1$Pa_OD`

#rename columns
colnames(Protease) <-c('Protease_raw','Strain','Pa_OD','Normalized_Protease')

#calculate mean and sd for all replicates, after normalizing by OD growth
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

proteasedata <- data_summary(Protease, varname="Normalized_Protease", 
                    groupnames=c("Strain")) 

proteasedata <- proteasedata[mixedorder(as.character(proteasedata$Strain)),]
rownames(proteasedata) = NULL
proteasedata <- proteasedata[c(15,3,4,5,6,2,7:14,16:18),]

proteasedata$Percent <- ((proteasedata$Normalized_Protease)/(proteasedata$Normalized_Protease[c(1)])*100)

#plot data, first make strain names a factor so it stays in order
proteasedata$Strain <- factor(proteasedata$Strain,levels = proteasedata$Strain)
p<- ggplot(proteasedata, aes(x=Strain, y=Normalized_Protease)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Normalized_Protease-sd, ymax=Normalized_Protease+sd), width=.2,
                position=position_dodge(.9)) 
print(p)

# Add labels and make x labels vertical so you can read them
p+labs(x="Strain", y = "Protease(µg/mL)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c('#999999','#E69F00'))

ggsave("Proteasecombopaper/paperfigure_05152022.png")

#if you want to output data you can do that here
write.table(proteasedata,
            "Proteasecombopaper/proteasedatacombo05152022.csv",
            sep=",",
            col.names=T,
            quote=F,
            row.names=F)
# 
# #determine significance using Welch t test
# PAO1 <- Protease$Normalized_Protease[Protease$Strain=="PAO1"]
# PAO1
# lasR_insert <- Protease$Normalized_Protease[Protease$Strain=="∆lasR insert"]
# lasR_insert
# lasR_clean <- Protease$Normalized_Protease[Protease$Strain=="∆lasR clean"]
# lasR_clean
# lasI <- Protease$Normalized_Protease[Protease$Strain=="∆lasI"]
# lasI
# lasRrhlR <- Protease$Normalized_Protease[Protease$Strain=="∆lasR∆rhlR"]
# lasRrhlR
# rhlR <- Protease$Normalized_Protease[Protease$Strain=="∆rhlR"]
# rhlR
# A17 <- Protease$Normalized_Protease[Protease$Strain=="A17"]
# A17
# SWPA15J <- Protease$Normalized_Protease[Protease$Strain=="SWPA15J = NSWPA15a"]
# SWPA15J
# CND03 <- Protease$Normalized_Protease[Protease$Strain=="CND03"]
# CND03
# e5BR2 <- Protease$Normalized_Protease[Protease$Strain=="5BR2"]
# e5BR2
# Jp238 <- Protease$Normalized_Protease[Protease$Strain=="Jp238"]
# Jp238
# Jp1155 <- Protease$Normalized_Protease[Protease$Strain=="Jp1155"]
# Jp1155
# W15Dec14 <- Protease$Normalized_Protease[Protease$Strain=="W15Dec14"]
# W15Dec14
# Aa249 <- Protease$Normalized_Protease[Protease$Strain=="Aa249"]
# Aa249
# PT31M <- Protease$Normalized_Protease[Protease$Strain=="PT31M"]
# PT31M
# JD303 <- Protease$Normalized_Protease[Protease$Strain=="JD303"]
# JD303
# 
# t.test(lasR_clean,lasI)

ls.Protease_byStrain = split(Protease, Protease$Strain)
ls.Protease_tTests = sapply(names(ls.Protease_byStrain), function(x){
  PAO1_in = ls.Protease_byStrain[['PAO1']]$Normalized_Protease
  Strain_in = ls.Protease_byStrain[[x]]$Normalized_Protease
  tTest_result = t.test(PAO1_in, Strain_in)
  tTest_pval = tTest_result$p.value
})

rownames(proteasedata) = as.character(proteasedata$Strain)
proteasedata$pval = ls.Protease_tTests[rownames(proteasedata)]
proteasedata$label = ifelse(proteasedata$pval<0.05, '*', '')

p_withstats<- ggplot(proteasedata, aes(x=Strain, y=Normalized_Protease)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Normalized_Protease-sd, ymax=Normalized_Protease+sd), width=.2,
                position=position_dodge(.9)) +
  labs(x="Strain", y = "Protease(µg/mL)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c('#999999','#E69F00')) +
  geom_text(aes(label=label), nudge_y = proteasedata$sd + 0.1, size=10)

print(p_withstats)
p_withstats
ggsave("Proteasecombopaper/paperfigure_09182022_withStats.png")
