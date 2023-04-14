head(phe)
phe$Sex<- as.factor(phe$Sex)
phe$Status<- as.factor(phe$Status)
phe$AgeCat<- as.factor(phe$AgeCat)
phe$RINCat<- as.factor(phe$RINCat)
phe$PMDCat<- as.factor(phe$PMDCat)
phe$Site_Specimen_Collected<- as.factor(phe$Site_Specimen_Collected)
dds <- DESeqDataSetFromTximport(txiGenes, colData = phe, design = ~ Sex + AgeCat + RINCat + PMDCat + Site_Specimen_Collected +  Status)
dds
head(phe)
dds <- estimateSizeFactors(dds)
idx <- rowSums( counts(dds, normalized=T) >= 5 ) >= 10
dds <- dds[idx,]
dds

library("sva")
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ Sex + AgeCat + RINCat + PMDCat + Site_Specimen_Collected +  Status, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
nsv <- num.sv(dat,mod)
nsv
svseq <- svaseq(dat, mod, mod0, n.sv = nsv)
dds$SV1 <- svseq$sv[,1]
dds$SV2 <- svseq$sv[,2]
design(dds)
design(dds) <- ~Sex + AgeCat + RINCat + PMDCat + Site_Specimen_Collected + SV1 + SV2 + Status

write.table(colData(dds), file = "samples.design.updated.231.txt",quote=FALSE,sep="\t")
write.table(counts(dds,normalized=TRUE), file = "normalisedGeneMatrix.txt",quote=FALSE,sep="\t")
write.table(counts(dds,normalized=FALSE), file = "notnormalisedGeneMatrix.txt",quote=FALSE,sep="\t")
head(colData(dds))
dds<-DESeq(dds)
resIHW <- results(dds, filterFun=ihw)
library(IHW)
resIHW <- results(dds, filterFun=ihw)
resIHWOrdered <- resIHW[order(resIHW$pvalue),]
sum(resIHW$padj < 0.05, na.rm=TRUE)
write.table(resIHWOrdered, "geneLevel_DE_res.txt",sep="\t")
