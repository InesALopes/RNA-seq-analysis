source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
install.packages("readxl")
biocLite("plyr")
biocLite("DESeq2")

library("biomaRt")
library(readxl)
library("plyr")
library("DESeq2")

mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")



#get normalized data for miRNA using DESeq2
directory<-"C:/Users/ines_/Desktop/Ana/miRNA/FilteredResults"

sampleFiles<-c("P1_1_Mature.txt", "P1_2_Mature.txt", "P2_1_Mature.txt", "P2_2_Mature.txt", 
               "P3_1_Mature.txt", "P3_2_Mature.txt", "P4_1_Mature.txt", "P4_2_Mature.txt", 
               "P5_1_Mature.txt", "P5_2_Mature.txt", "P6_1_Mature.txt", "P6_2_Mature.txt", 
               "P7_1_Mature.txt", "P7_2_Mature.txt", "P8_1_Mature.txt", "P8_2_Mature.txt", 
               "P9_1_Mature.txt", "P9_2_Mature.txt", "P10_1_Mature.txt", "P10_2_Mature.txt", 
               "P11_1_Mature.txt", "P11_2_Mature.txt", "P12_1_Mature.txt", "P12_2_Mature.txt")
sampleIndividual<-c("E_J_RC", "E_J_RC", "E_J_RC", "E_J_RC", "E_J_RC", "E_J_RC", "E_J_RC", "E_J_RC", 
                    "I_H_GI", "I_H_GI", "I_H_GI", "I_H_GI", "I_H_GI", "I_H_GI", "I_H_GI", "I_H_GI", 
                    "Is_C_R", "Is_C_R", "Is_C_R", "Is_C_R", "Is_C_R", "Is_C_R", "Is_C_R", "Is_C_R")
sampleTreatment<-c("CON", "CON", "Ab", "Ab", "M2", "M2", "M1", "M1", "CON", "CON", "Ab", "Ab", 
                   "M2", "M2", "M1", "M1", "CON", "CON", "Ab", "Ab", "M2", "M2", "M1", "M1")
sampleReplicate<-c("P1", "P1", "P2", "P2", "P3", "P3", "P4", "P4", "P5", "P5", "P6", "P6", "P7", 
                   "P7", "P8", "P8", "P9", "P9", "P10", "P10", "P11", "P11", "P12", "P12")

sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, 
                        treatment=sampleTreatment, individual=sampleIndividual, 
                        replicate=sampleReplicate)

dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~treatment)

dds<-collapseReplicates(dds, dds$replicate)


dds<- estimateSizeFactors(dds)
normmirna<-counts(dds, normalized=TRUE)



#get normalized data for mRNA using DESeq2
directory<-"C:/Users/ines_/Desktop/Ana/mRNA/Results"

sampleName<-c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12")
sampleFiles<-c("Human1.txt", "Human2.txt", "Human3.txt", "Human4.txt", "Human5.txt", "Human6.txt", "Human7.txt", "Human8.txt", "Human9.txt", "Human10.txt", "Human11.txt", "Human12.txt")
sampleIndividual<-c("I", "I", "I", "I", "E", "E", "E", "E", "J", "J", "J", "J")
sampleTreatment<-c("CON", "Ab", "M2", "M1", "CON", "Ab", "M2", "M1", "CON", "Ab", "M2", "M1")

sampleTable<-data.frame(sampleName=sampleName, fileName=sampleFiles, treatment=sampleTreatment, individual=sampleIndividual)

dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~treatment)

dds<- estimateSizeFactors(dds)
normmrna<-counts(dds, normalized=TRUE)


#get the list of DE miRNA
CON_M1_up=scan("C:/Users/ines_/Desktop/Ana/miRNA/FilteredResults/up_regulatedgenes_CONM1.txt", character(), quote="")
CON_M1_down=scan("C:/Users/ines_/Desktop/Ana/miRNA/FilteredResults/down_regulatedgenes_CONM1.txt", character(), quote="")
CON_M2_up=scan("C:/Users/ines_/Desktop/Ana/miRNA/FilteredResults/up_regulatedgenes_CONM2.txt", character(), quote="")
CON_M2_down=scan("C:/Users/ines_/Desktop/Ana/miRNA/FilteredResults/down_regulatedgenes_CONM2.txt", character(), quote="")
M1_M2_up=scan("C:/Users/ines_/Desktop/Ana/miRNA/FilteredResults/up_regulatedgenes_M1M2.txt", character(), quote="")
M1_M2_down=scan("C:/Users/ines_/Desktop/Ana/miRNA/FilteredResults/down_regulatedgenes_M1M2.txt", character(), quote="")
M1_Ab_up=scan("C:/Users/ines_/Desktop/Ana/miRNA/FilteredResults/up_regulatedgenes_M1Ab.txt", character(), quote="")
M1_Ab_down=scan("C:/Users/ines_/Desktop/Ana/miRNA/FilteredResults/down_regulatedgenes_M1Ab.txt", character(), quote="")
M2_Ab_up=scan("C:/Users/ines_/Desktop/Ana/miRNA/FilteredResults/up_regulatedgenes_M2Ab.txt", character(), quote="")
M2_Ab_down=scan("C:/Users/ines_/Desktop/Ana/miRNA/FilteredResults/down_regulatedgenes_M2Ab.txt", character(), quote="")

demirna<-unique(c(CON_M1_up, CON_M1_down, CON_M2_up, CON_M2_down, M1_M2_up, 
                  M1_M2_down, M1_Ab_up, M1_Ab_down, M2_Ab_up, M2_Ab_down))

mirnadata<-data.frame(normmirna)
mirnames<-rownames(mirnadata)
mirnadata$miRNA<-mirnames
mirnadata<-mirnadata[,c(13,1,5,6,7,8,9,10,11,12,2,3,4)]

mirnadata<-subset(mirnadata, mirnadata$miRNA %in% demirna)

#Check for targets
db<-read_excel("hsa_MTI.xlsx")
db2<-subset(db, db$`Support Type`=="Functional MTI")


targets<-subset(db, db$miRNA %in% demirna)
targets<-targets[, c("miRNA", "Target Gene")]

targets<-unique(targets)



#transform the mrna normalized list
mrnadata<-data.frame(normmrna)
mrnaIDs<-rownames(mrnadata)
mrnadata$"ensembl_gene_id"<-mrnaIDs

IDtransf<-as.matrix(mrnaIDs)
colnames(IDtransf)<-"ensembl_gene_id"
mrnanamespre<-getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                    filters = "ensembl_gene_id",
                    values = mrnaIDs,
                    mart = mart)
mrnanames<-merge(x = IDtransf, y = mrnanamespre, by="ensembl_gene_id",all.x=TRUE)

mrnadata<-merge(x=mrnadata, y=mrnanames, by="ensembl_gene_id",all.x=TRUE)

mrnadata<-mrnadata[,c(14,2,3,4,5,6,7,8,9,10,11,12,13)]
mrnadata<-rename(mrnadata, c("external_gene_name"="Target Gene"))

targetsinFile<-subset(mrnadata, mrnadata$`Target Gene` %in% targets$`Target Gene`)

#joining all files
df<-join_all(list(mirnadata, targets, mrnadata), type="inner")

#to check for DE genes
CON_M1_upgenes=scan("C:/Users/ines_/Desktop/Ana/mRNA/RResultsComp_FDR/up_regulatedgenes_CONM1.txt", character(), quote="")
Ab_M1_downgenes=scan("C:/Users/ines_/Desktop/Ana/mRNA/RResultsComp_FDR/down_regulatedgenes_M1Ab.txt", character(), quote="")
CON_M1_downgenes=scan("C:/Users/ines_/Desktop/Ana/mRNA/RResultsComp_FDR/down_regulatedgenes_CONM1.txt", character(), quote="")
Ab_M1_upgenes=scan("C:/Users/ines_/Desktop/Ana/mRNA/RResultsComp_FDR/up_regulatedgenes_M1Ab.txt", character(), quote="")
CON_M2_upgenes=scan("C:/Users/ines_/Desktop/Ana/mRNA/RResultsComp_FDR/up_regulatedgenes_CONM2.txt", character(), quote="")
Ab_M2_downgenes=scan("C:/Users/ines_/Desktop/Ana/mRNA/RResultsComp_FDR/down_regulatedgenes_M2Ab.txt", character(), quote="")
CON_M2_downgenes=scan("C:/Users/ines_/Desktop/Ana/mRNA/RResultsComp_FDR/down_regulatedgenes_CONM2.txt", character(), quote="")
Ab_M2_upgenes=scan("C:/Users/ines_/Desktop/Ana/mRNA/RResultsComp_FDR/up_regulatedgenes_M2Ab.txt", character(), quote="")
M1_M2_upgenes=scan("C:/Users/ines_/Desktop/Ana/mRNA/RResultsComp_FDR/up_regulatedgenes_M1M2.txt", character(), quote="")
M1_M2_downgenes=scan("C:/Users/ines_/Desktop/Ana/mRNA/RResultsComp_FDR/down_regulatedgenes_M1M2.txt", character(), quote="")

DEgenes<-unique(c(CON_M1_upgenes, CON_M1_downgenes, CON_M2_upgenes, 
                  CON_M2_downgenes, M1_M2_upgenes, M1_M2_downgenes, 
                  Ab_M1_upgenes, Ab_M1_downgenes, Ab_M2_upgenes, 
                  Ab_M2_downgenes))


DEgenesnames<-getBM(attributes = c("external_gene_name"),
                    filters = "ensembl_gene_id",
                    values = DEgenes,
                    mart = mart)
DEgenesnames<-DEgenesnames[,1]


#correlation test
output<-data.frame()

for (line in 1:dim(df)[1]){
  hap<-as.numeric(as.vector(df[line, 2:13]))
  hep<-as.numeric(as.vector(df[line, 15:26]))
  cor.result<-cor.test(hap, hep)
  IC<-cor.result$conf.int
  DEpresent<-df[line,14] %in% DEgenesnames
  temp<-data.frame(miRNA=as.vector(df[line,1]), gene=as.vector(df[line,14]), 
                   pvalue=cor.result$p.value, estimate_cor=cor.result$estimate, 
                   confInt_Inferior=IC[1], confInt_Superior=IC[2], 
                   isDE=DEpresent, statistic_t=cor.result$statistic, 
                   parameter=cor.result$parameter)
  output<-rbind(output, temp)
}


#Write the file with every data
write.table(output, file="TargetCorrelation_Everything.txt", sep="\t")
write.csv(output, file="TargetCorrelation_Everything.csv")

#Write the file for pvalue<0.05
outputSig<-subset(output, output$pvalue<0.05)
write.table(outputSig, file="TargetCorrelation_pvalue.txt", sep="\t")
write.csv(outputSig, file="TargetCorrelation_pvalue.csv")

#Write the file for pvalue<0.5 and cor<-0.5
outputCor<-subset(outputSig, outputSig$estimate_cor< -0.5)
outputCor<-outputCor[order(outputCor$estimate_cor),]
write.table(outputCor, file="TargetCorrelation_cor.txt", sep="\t")
write.csv(outputCor, file="TargetCorrelation_cor.csv")

#Write the file for targets that are DE
DEoutput<-subset(output, output$isDE==TRUE)
write.table(DEoutput, file="TargetCorrelationDE_Everything.txt", sep="\t")
write.csv(DEoutput, file="TargetCorrelationDE_Everything.csv")
outputSigDE<-subset(DEoutput, DEoutput$pvalue<0.05)
write.table(outputSigDE, file="TargetCorrelationDE_pvalue.txt", sep="\t")
write.csv(outputSigDE, file="TargetCorrelationDE_pvalue.csv")
outputCorDE<-subset(outputSigDE, outputSigDE$estimate_cor< -0.5)
outputCorDE<-outputCorDE[order(outputCorDE$estimate_cor),]
write.table(outputCorDE, file="TargetCorrelationDE_cor.txt", sep="\t")
write.csv(outputCorDE, file="TargetCorrelationDE_cor.csv")

#Check for each miRNA the total targets, those who are DE and their perc
mirnainoutput<-unique(output$miRNA)

counttable<-data.frame()

for (mirna in mirnainoutput) {
  totalTargets<-sum(output$miRNA==mirna)
  DETargets<-sum(output$miRNA==mirna & output$isDE==TRUE)
  PercDE<-(DETargets*100)/totalTargets
  temp<-data.frame(miRNA=mirna, TotalTargets=totalTargets, DETargets=DETargets, 
                   PercentageDE=PercDE)
  counttable<-rbind(counttable, temp)
}

write.table(counttable, file="TargetCorrelation_CountandPerc.txt", sep="\t")
write.csv(counttable, file="TargetCorrelation_CountandPerc.csv")

###
#General statistics

#ler ficheiros
targettotal<-read.csv("TargetCorrelation_Everything.csv")
targetpvalue<-read.csv("TargetCorrelation_pvalue.csv")

mirnaintotal<-as.vector(unique(targettotal$miRNA))
mirnainpvalue<-as.vector(unique(targetpvalue$miRNA))


CorStat<-function(corvalues, mirnalist){
  tabletotal<-data.frame()
  
  for(mirna in mirnalist){
    totalpositivecor<-length(which(corvalues$miRNA==mirna & 
                                     corvalues$estimate_cor>=0))
    totalnegativecor<-length(which(corvalues$miRNA==mirna & 
                                     corvalues$estimate_cor<0))
    totalall<-totalpositivecor+totalnegativecor
    positiveper<-(totalpositivecor*100)/totalall
    negativeper<-100-positiveper
    temp<-data.frame(miRNA=mirna, PositiveCorrelations=totalpositivecor, 
                     NegativeCorrelations=totalnegativecor, TotalCorrelations=totalall, 
                     PercentageOfPositive=positiveper, PercentageOfNegative=negativeper)
    tabletotal<-rbind(tabletotal, temp)
  }
  topositive<-sum(tabletotal$PositiveCorrelations)
  tonegative<-sum(tabletotal$NegativeCorrelations)
  toall<-topositive+tonegative
  toposiper<-(topositive*100)/toall
  tonegaper<-100-toposiper
  temp<-data.frame(miRNA="Total", PositiveCorrelations=topositive, 
                   NegativeCorrelations=tonegative, TotalCorrelations=toall, 
                   PercentageOfPositive=toposiper, PercentageOfNegative=tonegaper)
  tabletotal<-rbind(tabletotal, temp)
  return(tabletotal)
}

alltargets<-CorStat(targettotal, mirnaintotal)
write.table(alltargets, file="TargetCorrelation_GeneralStatAll.txt", sep="\t")
write.csv(alltargets, file="TargetCorrelation_GeneralStatAll.csv")

pvaluetargets<-CorStat(targetpvalue, mirnainpvalue)
write.table(pvaluetargets, file="TargetCorrelation_GeneralStatPvalue.txt", sep="\t")
write.csv(pvaluetargets, file="TargetCorrelation_GeneralStatPvalue.csv")
