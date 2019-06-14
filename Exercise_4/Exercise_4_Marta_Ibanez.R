## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ------------------------------------------------------------------------
expression <- as.matrix(read.delim("~/Desktop/Universitat/2_ANY/3_TRIMESTRE/Omics Techniques/Ibanez_Marta-OmicsTechniques/Exercise_4/Data/expression.txt", row.names=1, comment.char="#"))
targets <- read.delim("~/Desktop/Universitat/2_ANY/3_TRIMESTRE/Omics Techniques/Ibanez_Marta-OmicsTechniques/Exercise_4/Data/targets.txt", header=TRUE)
require(Biobase)


## ------------------------------------------------------------------------
show(targets)
dim(expression)
summary(expression)
class(targets)
class(expression)
mypatient <- paste0("Patient", 1:11)


## ------------------------------------------------------------------------
myInfo = list(myName= "Marta Ibáñez", myLab = "Bioinformatics Lab", myContact = "marta.ibanez@alum.esci.upf.edu", myTitle = "Practical Exercise on ExpressionSets")
show(myInfo)


## ------------------------------------------------------------------------
myEset <- ExpressionSet(expression)
class(myEset)
show(myEset)


## ------------------------------------------------------------------------
columnDesc <- data.frame(labelDescription = c("Sample Names", "Patients", "Groups"))
myAnnotDF <- new("AnnotatedDataFrame", data = targets, varMetadata= columnDesc)
show(myAnnotDF)


## ------------------------------------------------------------------------
phenoData(myEset) <- myAnnotDF
rownames(pData(myAnnotDF)) <- pData(myAnnotDF)$SampleName
myEset <- ExpressionSet(assayData = expression, phenoData = myAnnotDF)
show(myEset)


## ------------------------------------------------------------------------
myDesc <- new("MIAME", name= myInfo[["myName"]], lab= myInfo[["myLab"]], contact = myInfo[["myContact"]], title = myInfo[["myTitle"]])
print(myDesc)


## ------------------------------------------------------------------------
myEset <- ExpressionSet(assayData = expression, fetureNames = mypatient, phenoData = myAnnotDF, experimentData = myDesc)
show(myEset)


## ------------------------------------------------------------------------
workingDir <-getwd()
dataDir <- file.path(workingDir, "Data")
resultsDir <- file.path(workingDir, "Results")


## ------------------------------------------------------------------------
targets <-read.csv(file=file.path(dataDir,"targets.csv"), header = TRUE, sep=";")
targets


## ------------------------------------------------------------------------
installifnot <- function (pkg){
  if (!require(pkg, character.only=T)){
    BiocManager::install(pkg)
}else{
  require(pkg, character.only=T)
  }
}
installifnot("oligo")
installifnot("limma")
installifnot("Biobase")
installifnot("arrayQualityMetrics")
installifnot("genefilter")
installifnot("multtest")
installifnot("annotate")
installifnot("xtable")
installifnot("gplots")
installifnot("scatterplot3d")
#install.packages("gridSVG")
CELfiles <- list.celfiles(file.path(dataDir)) #different cel files in data
CELfiles
rawData <- read.celfiles(file.path(dataDir,CELfiles)) #cel files in another way


## ------------------------------------------------------------------------
sampleNames<- as.character(targets$ShortName)
sampleNames
sampleColor <- as.character(targets$Color)
sampleColor


## ------------------------------------------------------------------------
boxplot(rawData, which="all",las=2, main="Intensity distribution of RAW data", 
        cex.axis=0.6, names = sampleNames, col=sampleColor)


## ------------------------------------------------------------------------
clust.euclid.average <- hclust(dist(t(exprs(rawData))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of RawData", 
     cex=0.7,  hang=-1)


## ------------------------------------------------------------------------
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE, formapunts=NULL, myCex=0.8,...)
{
  pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, pch=formapunts, 
       xlim=c(min(pcX$x[,1])-100000, max(pcX$x[,1])+100000),ylim=c(min(pcX$x[,2])-100000, max(pcX$x[,2])+100000))
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=myCex)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}
plotPCA(exprs(rawData), labels=sampleNames, dataDesc="raw data", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)


## ------------------------------------------------------------------------
pdf(file.path(resultsDir, "QCPlots_Raw.pdf"))
boxplot(rawData, which="all",las=2, main="Intensity distribution of RAW data", 
        cex.axis=0.6, col=sampleColor, names=sampleNames)
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of samples of RawData", 
     cex=0.7,  hang=-1)
plotPCA(exprs(rawData), labels=sampleNames, dataDesc="raw data", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)
dev.off()


## ------------------------------------------------------------------------
eset<-rma(rawData)

write.exprs(eset, file.path(resultsDir, "NormData.txt"))


## ------------------------------------------------------------------------
boxplot(eset, las=2, main="Intensity distribution of Normalized data", cex.axis=0.6, 
        col=sampleColor, names=sampleNames)


## ------------------------------------------------------------------------
clust.euclid.average <- hclust(dist(t(exprs(eset))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of NormData", 
     cex=0.7,  hang=-1)


## ------------------------------------------------------------------------
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE, formapunts=NULL, myCex=0.8,...)
{
  pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, pch=formapunts, 
       xlim=c(min(pcX$x[,1])-10, max(pcX$x[,1])+10),ylim=c(min(pcX$x[,2])-10, max(pcX$x[,2])+10))
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=myCex)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}
plotPCA(exprs(eset), labels=sampleNames, dataDesc="NormData", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)


## ------------------------------------------------------------------------
pdf(file.path(resultsDir, "QCPlots_Norm.pdf"))
boxplot(eset, las=2, main="Intensity distribution of Normalized data", cex.axis=0.6, 
        col=sampleColor, names=sampleNames)
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of NormData", 
     cex=0.7,  hang=-1)
plotPCA(exprs(eset), labels=sampleNames, dataDesc="selected samples", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)
dev.off()


## ------------------------------------------------------------------------

annotation(eset) <- "org.Hs.eg.db"
eset_filtered <- nsFilter(eset, var.func=IQR,
         var.cutoff=0.75, var.filter=TRUE,
         filterByQuantile=TRUE)
#NUMBER OF GENES OUT
print(eset_filtered$filter.log$numLowVar)

#NUMBER OF GENES IN
print(eset_filtered$eset)


## ------------------------------------------------------------------------
design <-model.matrix(~ 0+targets$Group)
colnames(design)<-c("ASC", "BM_BSC")
rownames(design)<- targets$ShortName
print(design); 


## ------------------------------------------------------------------------
cont.matrix <- makeContrasts (
  ASC_BMBSC = ASC-BM_BSC,
  levels=design)
comparison1 <- "Effect of cells"


## ------------------------------------------------------------------------
fit1 <- lmFit(eset_filtered$eset, design) # Fits a linear model, solving the equations to obtain the estimates. Estimates parameters using least sqaures
fit.main1 <- contrasts.fit(fit1, cont.matrix)
fit.main1 <- eBayes(fit.main1)
fit.main1


## ------------------------------------------------------------------------
topTab <-  topTable(fit.main1, number=nrow(fit.main1), coef=NULL, adjust="fdr")
write.csv2(topTab, file= file.path(resultsDir,paste("Selected.Genes.in.comparison.",
                                                    comparison1, ".csv", sep = "")))
print(xtable(topTab,align="lllllll"),type="html",html.table.attributes="",
      file=paste("Selected.Genes.in.comparison.",comparison1,".html", sep=""))


## ------------------------------------------------------------------------
volcanoplot(fit.main1, highlight=10, names=fit.main1$ID, 
            main = paste("Differentially expressed genes", colnames(cont.matrix), sep="\n"))



pdf(file.path(resultsDir,"Volcanos.pdf"))
volcanoplot(fit.main1, highlight = 10, names = fit.main1$ID, 
            main = paste("Differentially expressed genes", colnames(cont.matrix), sep = "\n"))



## ------------------------------------------------------------------------
all_anota<-data.frame(exprs(eset))
library(hugene21sttranscriptcluster.db)
Annot <- data.frame(SYMBOL=sapply(contents(hugene21sttranscriptclusterSYMBOL), paste, collapse=", "),
                    DESC=sapply(contents(hugene21sttranscriptclusterGENENAME), paste, collapse=", "))
Annot<-Annot[!Annot$SYMBOL=="NA",]
Annot<-Annot[!Annot$DESC=="NA",]
head(Annot)

anotaGenes <- merge(Annot,all_anota, by.x=0,by.y=0)
head(anotaGenes)
write.table(anotaGenes, file ="data.ann.txt",sep="\t")

rownames(anotaGenes) <- anotaGenes[,1]
anotaGenes <- anotaGenes[,-1]
anotaGenes.end <- merge(anotaGenes, topTab, by.x=0,by.y=0)

# topTab.end = anotaGenes.end
topTab.end <- anotaGenes.end[order(-anotaGenes.end$B),]

rownames(topTab.end) <- topTab.end[,1]
topTab.end <- topTab.end[, -1]
write.csv(topTab.end, file = file.path(resultsDir,"TopTable.end.csv"))

