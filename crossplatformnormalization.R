##1.SERIES_MATRIX_DATA_INPUT#################################################################################

#Load series matrix dataset in R.

GPL96 <- read.table(file="GSE9006_96series.txt",comment.char='!',sep='\t',row.names=NULL,header=T)

##2.STRING_SPLIT#############################################################################################

library(tidyr)
library(dplyr)

a <- GPL96 %>% 
  mutate(SYMBOL = strsplit(as.character(SYMBOL), "///")) %>% 
  unnest(SYMBOL)

write.table(a, "seriesgenesplit_96.txt", sep="\t")

##3.FILTER_AND_AVEREPS#############################################################################################

GPL96 <- read.table(file="seriesgenesplit_96.txt",comment.char='!',sep='\t',row.names=NULL,header=T)

GPL96$ID_REF <-  NULL

library(dplyr)
GPL96 <- filter(GPL96, SYMBOL != '0')
GPL96 <- filter(GPL96, SYMBOL != '#N/A')
GPL96 <- filter(GPL96, SYMBOL != ' ')


#summarize probes and remove replicates
library(limma)

GPL96 <- avereps(GPL96,ID=GPL96$SYMBOL)

write.table(GPL96, "GPL96avg.txt", sep="\t")

GPL96 <- as.data.frame(GPL96)

#load mergeavg in R

GPL96 <- read.table(file="GPL96avg.txt",comment.char='!',sep='\t',row.names=NULL,header=T)

rownames(GPL96) <- GPL96$SYMBOL
GPL96$SYMBOL <- NULL

boxplot(GPL96)

##DENSITYPLOT######################################################################################

colramp = colorRampPalette(c(3,"white",2))(20)
plot(density(GPL96[,1]),col=colramp[1],lwd=3,ylim=c(0,.50))
for(i in 1:36){lines(density(GPL96[,i]),lwd=3,col=colramp[i])}

##4.BOXCOX###########################################################################################
#transforms data so that it closely resembles a normal distribution

library(bestNormalize)

GPL96 <- GPL96 +1 # added 1 to avoid negative value -0.2013121
GPL962 <- GPL96  # making other dataframe to store new transformed values

for(i in 1:ncol(GPL96)){
  print(i)
  transformed <- boxcox(as.vector(GPL96[,i]))
  GPL962[,i]<-transformed$x.t
}

write.csv(GPL962, file = "GPL962_boxcox.csv") 

##DENSITYPLOT_after_boxcox##############################################################################

colramp = colorRampPalette(c(3,"white",2))(20)
plot(density(GPL962[,1]),col=colramp[1],lwd=3,ylim=c(0,.50))
for(i in 1:36){lines(density(GPL962[,i]),lwd=3,col=colramp[i])}

boxplot(GPL962)

##PCA############################################################################################

#pca analysis
library(devtools)
library(Biobase)
library(preprocessCore)
library(sva)
library(bladderbatch)
library(snpStats)

dim(GPL962)

datapc <- GPL962

pc.data <- princomp(datapc,cor=TRUE)

#information output
names(pc.data)

summary(pc.data)

#eigenvalues/eigenvectors
eigenvectors <- pc.data$loadings

#note:- these value are scaled so the ss=1.
eigenvalues <- pc.data$sdev*pc.data$sdev

round(cor(GPL962,pc.data$scores),3)

library("factoextra")
library("FactoMineR")

pca.data <- PCA(GPL962, scale.unit = TRUE, graph = FALSE)

fviz_eig(pca.data, addlabels = TRUE, ylim = c(0, 70))

fviz_pca_var(pca.data, col.var = "cos2",
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE) 

pca.data <- PCA(t(GPL962), scale.unit = TRUE, graph = FALSE)

fviz_pca_ind(pca.data, col.ind = "cos2", 
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
             repel = TRUE)

##5.DESIGN_MATRIX##################################################################################

ss <- GPL962
#design matrix for case and control

library(limma)

design <- matrix(data=0,nrow=36,ncol=2)
design

rownames(design) <- colnames(ss)
design

colnames(design) <- c('control','case')
design

design[1:24,'control']<-1
design[25:36,'case']<-1

design

##6.LMFIT_AND_ EBAYES################################################################################

#extracting DEGs

fit <-lmFit(ss,design)

contrast.matrix <- makeContrasts(case-control,levels=design)

fit2 <- contrasts.fit(fit,contrast.matrix)

fit2 <- eBayes(fit2)

tt <- topTable(fit2,coef=1,adjust="BH",number=nrow(ss))

write.table(tt, "ttseries96.txt", sep="\t")

