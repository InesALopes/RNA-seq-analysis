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
directory<-"C:/Users/ines_/Desktop/Ana/miRNA/FilteredResults"

sampleFiles<-c("P1_1_Mature.txt", "P1_2_Mature.txt", "P2_1_Mature.txt", "P2_2_Mature.txt", 
               "P3_1_Mature.txt", "P3_2_Mature.txt", "P4_1_Mature.txt", "P4_2_Mature.txt", 
               "P5_1_Mature.txt", "P5_2_Mature.txt", "P6_1_Mature.txt", "P6_2_Mature.txt", 
               "P7_1_Mature.txt", "P7_2_Mature.txt", "P8_1_Mature.txt", "P8_2_Mature.txt", 
               "P9_1_Mature.txt", "P9_2_Mature.txt", "P10_1_Mature.txt", "P10_2_Mature.txt", 
               "P11_1_Mature.txt", "P11_2_Mature.txt", "P12_1_Mature.txt", "P12_2_Mature.txt")
sampleIndividual<-c("Pool1", "Pool1", "Pool1", "Pool1", "Pool1", "Pool1", "Pool1", "Pool1", 
                    "Pool2", "Pool2", "Pool2", "Pool2", "Pool2", "Pool2", "Pool2", "Pool2", 
                    "Pool3", "Pool3", "Pool3", "Pool3", "Pool3", "Pool3", "Pool3", "Pool3")
sampleTreatment<-c("CON", "CON", "Ab", "Ab", "M2", "M2", "M1", "M1", "CON", "CON", "Ab", "Ab", 
                   "M2", "M2", "M1", "M1", "CON", "CON", "Ab", "Ab", "M2", "M2", "M1", "M1")
sampleReplicate<-c("P1", "P1", "P2", "P2", "P3", "P3", "P4", "P4", "P5", "P5", "P6", "P6", "P7", 
                   "P7", "P8", "P8", "P9", "P9", "P10", "P10", "P11", "P11", "P12", "P12")

sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, 
                        treatment=sampleTreatment, individual=sampleIndividual, 
                        replicate=sampleReplicate)
sampleTable$treatment<-relevel(sampleTable$treatment, ref="CON")

dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~treatment)

#Collapse the replicates
dds<-collapseReplicates(dds, dds$replicate)



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
res<-results(ddssva, contrast=c("treatment", "CON", "Ab"), alpha = 0.05)


#results summary
summary(res)


#Makin an heatmap
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
