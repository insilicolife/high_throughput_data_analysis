#Task 1
sampleData<-read.table("gene_expr_microarrays/sample_data.txt", row.names=1,sep="\t")
sampleHeader<-read.table("gene_expr_microarrays/sample_header.txt", skip=2, header=FALSE,sep="\t")
sampleId<-unlist(sampleHeader[1,2:83])
colnames(sampleData)<-sampleId
probIds<-row.names(sampleData)
source("http://bioconductor.org/biocLite.R")
biocLite("hgu133plus2.db")
library("hgu133plus2.db")
geneName<-select(hgu133plus2.db, probIds, "ALIAS" , "PROBEID")
notDuplicatedByProbId<-geneName[!duplicated(geneName[,1]),]
#duplicatedGeneName<-geneName[duplicated(geneName[,2]),]
notDuplicatedByGeneNameProbId<-notDuplicatedByProbId[!duplicated(notDuplicatedByProbId[,2]),]
newDataMatrix<-matrix(nrow=dim(notDuplicatedByGeneNameProbId)[1],ncol=length(colnames(sampleData)))
for(i in 1:dim(notDuplicatedByGeneNameProbId)[1]){
	#search for all duplicated genes associated with the nonduplicated gene
	dupProbId<-notDuplicatedByProbId[grep(notDuplicatedByGeneNameProbId[i,2],notDuplicatedByProbId[,2]),1]
	#add the non duplicated gene with duplicated one
	#dupProbId<-c(dupProbId,notDuplicatedByGeneNameProbId[i,1])
	#replace the expression value by the avarage expression values of those duplicated genes for each sample
	newDataMatrix[i,]<-apply(sampleData[dupProbId,],2,mean)
}
summerizedSampleData<-newDataMatrix
rownames(summerizedSampleData)<-notDuplicatedByGeneNameProbId[,2]
colnames(summerizedSampleData)<-colnames(sampleData)
#print(head(sampleData))
#dim(sampleData)
#task 3
#before log transform
boxplot(summerizedSampleData)
plot(summerizedSampleData)
hist(as.matrix(summerizedSampleData))
#log transforming
logSummerizedSampleData<-log(summerizedSampleData)
boxplot(logSummerizedSampleData)# after log transformation
plot(logSummerizedSampleData)
hist(logSummerizedSampleData)
#quartile normalization
library(preprocessCore)
normalizedLogSummerizedSampleData<-normalize.quantiles(as.matrix(logSummerizedSampleData),copy=TRUE)
rownames(normalizedLogSummerizedSampleData)<-rownames(logSummerizedSampleData)
colnames(normalizedLogSummerizedSampleData)<-colnames(logSummerizedSampleData)
boxplot(normalizedLogSummerizedSampleData)
plot(normalizedLogSummerizedSampleData)
hist(as.matrix(normalizedLogSummerizedSampleData))
#Task 4 deferentially expression
#ETS rearrangement-positive sample group

ETV1ProbId<- intersect(geneName[which(geneName$ALIAS=="ETV1"),1],notDuplicatedByGeneNameProbId[,1])
ETV1GeneinSummerizedGeneSet<-notDuplicatedByGeneNameProbId[notDuplicatedByGeneNameProbId$PROBEID==ETV1ProbId,2]
ETV4ProbId<- intersect(geneName[which(geneName$ALIAS=="ETV4"),1],notDuplicatedByGeneNameProbId[,1])
ETV4GeneinSummerizedGeneSet<-notDuplicatedByGeneNameProbId[notDuplicatedByGeneNameProbId$PROBEID==ETV4ProbId,2]
#ETV1ProbId<-notDuplicatedGeneNames[notDuplicatedGeneNames$ALIAS=="ETV1",1]
#ETV1ProbId<-ETV1ProbId[!is.na(ETV1ProbId)]
#ETV4ProbId<-notDuplicatedGeneNames[notDuplicatedGeneNames$ALIAS=="ETV4",1]
#ETV4ProbId<-ETV4ProbId[!is.na(ETV4ProbId)]
EST_Rearrangene<-c(ETV1GeneinSummerizedGeneSet,ETV4GeneinSummerizedGeneSet)

highlyExpressdETV4<-names(which(normalizedLogSummerizedSampleData[EST_RearrangeProbIds[2],]>6))
highlyExpressdETV1<-names(which(normalizedLogSummerizedSampleData[EST_RearrangeProbIds[1],]>7.1))
highlyExpressedEST_RearrangeSamples<-c(highlyExpressdETV1,highlyExpressdETV4)

#ETV1ProbId<- intersect(geneName[which(geneName$ALIAS=="ETV1"),1],notDuplicatedByGeneNameProbId[,1])
#ETV1GeneinSummerizedGeneSet<-notDuplicatedByGeneNameProbId[notDuplicatedByGeneNameProbId$PROBEID==ETV1ProbId,2]

#plot(normalizedLogSummerizedSampleData)
#SPINK1-overexpression sample group
#SPINK1ProbId<-notDuplicatedGeneNames[notDuplicatedGeneNames$ALIAS=="SPINK1",1]
#SPINK1ProbId<-SPINK1ProbId[!is.na(SPINK1ProbId)]

SPINK1ProbId<- intersect(geneName[which(geneName$ALIAS=="SPINK1"),1],notDuplicatedByGeneNameProbId[,1])
SPINK1GeneinSummerizedGeneSet<-notDuplicatedByGeneNameProbId[notDuplicatedByGeneNameProbId$PROBEID==SPINK1ProbId,2]

highlyExpressedSPINK1Sample<-names(which(normalizedLogSummerizedSampleData[SPINK1GeneinSummerizedGeneSet,]>7.8))
#remaining sample group
allSamples<-colnames(normalizedLogSummerizedSampleData)
knownSamples<-c(highlyExpressedSPINK1Sample,highlyExpressedEST_RearrangeSamples)
unknownSamples<-setdiff(allSamples,knownSamples)
#Defferentialy exprestion calculation using Lima package
biocLite("limma")
library("limma")
samples<-c(rep(1,82))
EST_RearrangeSamples<-match(highlyExpressedEST_RearrangeSamples,colnames(normalizedLogSummerizedSampleData))
SPINK1Sample<-match(highlyExpressedSPINK1Sample,colnames(normalizedLogSummerizedSampleData))
samples[EST_RearrangeSamples]<-2
samples[SPINK1Sample]<-3
design <- model.matrix(~ 0+factor(samples))
colnames(design)<-c("unknownSample","EST_Samples","SPINK1Sample")
fit <- lmFit(normalizedLogSummerizedSampleData, design)

contrast.matrix <- makeContrasts(unknownSample-EST_Samples, SPINK1Sample-EST_Samples, SPINK1Sample-unknownSample, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
top_DE_gene_unknown_EST<-topTable(fit2, coef=1, adjust="BH")
top_DE_gene_SPINK1_EST<-topTable(fit2, coef=2, adjust="BH")
top_DE_gene_SPINK1_unknown<-topTable(fit2, coef=3, adjust="BH")
results <- decideTests(fit2)
#vendiagram to show the relation between thãe sample groups
venDiagramComparision<-vennDiagram(results)

#cdna_microarray

CRPC_278<-read.csv("cdna_microarrays/CRPC_278.csv", skip=1, sep="\t")
CRPC_543<-read.csv("cdna_microarrays/CRPC_543.csv", skip=1, sep="\t")
PCaN_5934<-read.csv("cdna_microarrays/PCaN_5934.csv", skip=1, sep="\t")
PCaN_6102<-read.csv("cdna_microarrays/PCaN_6102.csv", skip=1, sep="\t")
PCaP_17163<-read.csv("cdna_microarrays/PCãAaP_17163.csv", skip=1, sep="\t")
PCaP_470<-read.csv("cdna_microarrays/PCaP_470.csv", skip=1, sep="\t")
hist(CRPC_278)

petients<-data.frame(gene=CRPC_278[,1],
			petients1=CRPC_278[,2],
			petients2=CRPC_543[,2],
			petients3=PCaN_5934[,2],
			petients4=PCaN_6102[,2],
			petients5=PCaP_17163[,2],
			petients6=PCaP_470[,2])
install.packages("ggplot2")
install.packages("reshape")
library(reshape)
petients <- melt(petients,  id = 'gene', variable_name = 'genes')
histComparision<-ggplot(petients, aes(value,colour=genes)) + geom_density()
#quartile normalization
petientsNormalized<-normalize.quantiles(as.matrix(petients$value),copy=TRUE)
petients$value<-petientsNormalized
histComparissionNormalized<-ggplot(petients, aes(value,colour=genes)) + geom_density()

#expression level of the gene ERG
ERGexpressionLevele<-as.matrix(petients[(petients$gene=="ERG"),3])
colnames(ERGexpressionLevele)<-"ERGEspression"
rownames(ERGexpressionLevele)<-c("petient1","petient2","petient3","petient4","petient5","petient6")

ERGElivatedExpresedPeteints<- ERGexpressionLevele[ERGexpressionLevele[,1]>11,]
ERGfamilyGene<-petients[grep("ERG",petients$gene),1]
#ETS family of transcription factor
PEA3SubFamily<-c("ETV1","ETV4","ETV5")
PEA3SubFamilyExpressionLevele<-data.frame(data=NA)
for(i in 1:length(PEA3SubFamily)){

	PEA3SubFamilyExpressionLevele<-cbind(PEA3SubFamilyExpressionLevele,petients[petients$gene==PEA3SubFamily[i],3])
}
PEA3SubFamilyExpressionLevele<-subset(PEA3SubFamilyExpressionLevele, select=-(data))
colnames(PEA3SubFamilyExpressionLevele)<-PEA3SubFamily
rownames(PEA3SubFamilyExpressionLevele)<-c("petient1","petient2","petient3","petient4","petient5","petient6")
overEXpressedPEA3FamilySample<-c()
for(j in 1:ncol(PEA3SubFamilyExpressionLevele)){

		overEXpressedPEA3FamilySample<-c(overEXpressedPEA3Family,rownames(PEA3SubFamilyExpressionLevele[PEA3SubFamilyExpressionLevele[,j]>13,]))
}

#over expressed AR

#Expresion level of the gene ERG
ARexpressionLevele<-as.matrix(petients[(petients$gene=="ERG"),3])
colnames(ARexpressionLevele)<-"ARExpression"
rownames(ARexpressionLevele)<-c("petieãOAnt1","petient2","petient3","petient4","petient5","petient6")
