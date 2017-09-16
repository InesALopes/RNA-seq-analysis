source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("pheatmap")
biocLite("sva")
biocLite("genefilter")


library("DESeq2")
library("pheatmap")
library("sva")
library("genefilter")


#First define the directory where the HtSeq output is located
directory<-"C:/Users/ines_/Desktop/Ana/mRNA/Results"

sampleName<-c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12")
sampleFiles<-c("Human1.txt", "Human2.txt", "Human3.txt", "Human4.txt", "Human5.txt", "Human6.txt", "Human7.txt", "Human8.txt", "Human9.txt", "Human10.txt", "Human11.txt", "Human12.txt")
sampleIndividual<-c("Pool3", "Pool3", "Pool3", "Pool3", "Pool1", "Pool1", "Pool1", "Pool1", "Pool2", "Pool2", "Pool2", "Pool2")
sampleTreatment<-c("CON", "Ab", "M2", "M1", "CON", "Ab", "M2", "M1", "CON", "Ab", "M2", "M1")

sampleTable<-data.frame(sampleName=sampleName, fileName=sampleFiles, treatment=sampleTreatment, individual=sampleIndividual)
sampleTable$treatment<-relevel(sampleTable$treatment, ref="CON")

dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~treatment)


#Choosing the DAVID background gene list
#Genes must have mapped at least 5 times for one condition
background<-counts(dds)
subGenes <- apply(background, 1, function(x) { any(x > 5)})
background<-background[subGenes==TRUE,]
writeLines(rownames(background), con="DAVIDBackground.txt")


#Pre-filtering to remove genes with 0 or 1 counts
dds<-dds[rowSums(counts(dds))>1,]

#DE analysis
dds<-DESeq(dds)

#PCA graph
rld <- rlog(dds, blind=FALSE)

plotPCA(rld, intgroup="treatment")


plotPCA(rld, intgroup="individual")


#Batch effect
dat <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~ treatment, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv=2)


ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + treatment

ddssva <- DESeq(ddssva)

#Apply any contrast here
res<-results(ddssva, contrast=c("treatment", "M1", "M2"), alpha = 0.05)


#results summary
summary(res)

#Making an heatmap

rld2 <- rlog(ddssva, blind=FALSE)
topVarGenes <- head(order(rowVars(assay(rld2)), decreasing = TRUE), 20)
mat  <- assay(rld2)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld2)[, c("individual","treatment")])
pheatmap(mat, annotation_col = anno)


#Significant genes
resSig=subset(res, padj<0.05)


#down-regulated
predown <- subset(resSig, log2FoldChange<= -1)
down<-predown[order(predown$log2FoldChange),]

write.table(down, file="down_regulatedeverything_M2Ab.txt", sep="\t")
writeLines(rownames(down), con="down_regulatedgenes_M2Ab.txt")

#up-regulated
preup <- subset(resSig, log2FoldChange>= 1)
up<-preup[order(preup$log2FoldChange, decreasing = TRUE),]

write.table(up, file="up_regulatedeverything_M2Ab.txt", sep="\t")
writeLines(rownames(up), con="up_regulatedgenes_M2Ab.txt")
