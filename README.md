

## Prepare txdb
```
#install if needed
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("tximport")
BiocManager::install("GenomicFeatures")

library(tximport)
library(tidyverse)
library(GenomicFeatures)
library(readr)

txdb <- makeTxDbFromGFF("~/Documents/SalmonTarget/gencode.v38.chr_patch_hapl_scaff.annotation.gtf")
txdb.filename <- "gencode.v38.hr_patch_hapl_scaff.annotation.sqlite"
saveDb(txdb, txdb.filename)
#if you are starting from the db
#txdb<-loadDb("~/Documents/swimmongDownstream/gencode.v38.hr_patch_hapl_scaff.annotation.sqlite")

txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
tab
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]
head(txdf)
dim(txdf)
columns(txdb)
keytypes(txdb)
k <- keys(txdb, keytype = "GENEID")
tx2gene <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene<-tx2gene[,c(2,1)]
write.table(as.data.frame(tx2gene), file="tx2gene.txt", sep="\t")
```

## Get phenotype data for the motor cortex region - with new rules for cases/controls, see https://github.com/rkabiljo/editPhenoTargetALS

## Read with tximport
```
setwd("~/Documents/MNDquantNew")
#specify which files to read -for example whatever is in my motor cortex matrix
#this is to see what's in the directlry
fileList<-list.files(pattern = "quant.sf", recursive = TRUE)
#and this is from the phenotype table
files<-paste0(rownames(pheCortexAll),"/quant.sf")
all(files %in% fileList)
#[1] TRUE

txiGenes <- tximport(files, type="salmon", tx2gene=tx2gene, txOut = FALSE)
dim(txiGenes$counts)
#[1] 60230   272
head(txiGenes$counts)
colnames(txiGenes$counts)<-rownames(pheCortexAll)


#txiTranscripts <- tximport(files, type="salmon", txOut = TRUE)
txiTranscripts <- tximport(files, type="salmon", txOut =TRUE,countsFromAbundance = "scaledTPM",tx2gene=tx2gene)
head(txiTranscripts$counts)
dim(txiTranscripts$counts)
#[1] 236186    272
colnames(txiTranscripts$counts)<-rownames(pheCortexAll)
txiTranscripts$counts[1:5,1:5]
```
## another way to read - summarise transcripts to genes later
```
txi.sum <- summarizeToGene(txiTranscripts, tx2gene)
#check that it's the same as before
all.equal(txiGenes$counts, txi.sum$counts)
```
## or different flags to summarise transcripts, depends on the needs of post analysis
```
txiTranscripts <- tximport(files, type="salmon", txOut =TRUE,countsFromAbundance = "dtuScaledTPM",tx2gene=tx2gene )
```


## DESeq2 with the gene tximport object
see DESeq2_afterTximport.R     To compare these results with our previous DE results on the same dataset see compareDESeqResWithDEPaper.R


# DRIMSeq
## prepare for DRIMSeq
```
cts <- txiTranscripts$counts
dim(cts)
#[1] 236186    272
cts <- cts[rowSums(cts) > 0,]
dim(cts)
#[1] 216530    272
cts[1:3,1:3]
#if needed, set names
#colnames(cts)<-samples3$sample_name
#cts[1:3,1:3]
range(colSums(cts)/1e6)
all(rownames(cts) %in% txdf$TXNAME)
txdf <- txdf[match(rownames(cts),txdf$TXNAME),]
all(rownames(cts) == txdf$TXNAME)
counts <- data.frame(gene_id=txdf$GENEID,
feature_id=txdf$TXNAME,
cts)
#hyphens were replaced by ., and I am reverting to hyphens
 colnames(counts) <- gsub("\\.", "-", colnames(counts))
```
## actual DRIMSeq
```
BiocManager::install("DRIMSeq")
library(DRIMSeq)
head(phe)
table(phe$Status)
 samps<-phe
samps$sample_id <- rownames(samps)
samps<-samps[,c("sample_id","Status")]
colnames(samps)<-c("sample_id","condition")
samps$condition<-as.factor(samps$condition)
head(samps)

d <- dmDSdata(counts=counts, samples=samps)
d
#An object of class dmDSdata 
#with 50632 genes and 272 samples
#* data accessors: counts(), samples()
counts(d[1,])[,1:4]
#            gene_id        feature_id CGND.HRA.00596 CGND.HRA.00646
#1 ENSG00000223972.5 ENST00000456328.2              0       1.115553

n<-272
#n.small<-36
n.small<-37

d <- dmFilter(d,
               min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
               min_samps_gene_expr=n, min_gene_expr=10)
d
#An object of class dmDSdata 
#with 10266 genes and 272 samples
#* data accessors: counts(), samples()

table(table(counts(d)$gene_id))

#    2    3    4    5    6    7    8    9   10   11   12   14 
#2259 2411 2119 1596 1005  514  220   97   33    7    4    1 
 
design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
colnames(design_full)
#[1] "(Intercept)" "condition2" 
set.seed(1) 
system.time({
     d <- dmPrecision(d, design=design_full)
     d <- dmFit(d, design=design_full)
     d <- dmTest(d, coef="condition2")
 })
#! Using a subset of 0.1 genes to estimate common precision !
#! Using common_precision = 13.1556 as prec_init !
#! Using 0 as a shrinkage factor !

#    user   system  elapsed 
# 987.215   22.163 1012.409 
res <- DRIMSeq::results(d)
head(res)
```
Output:
```
#              gene_id         lr df      pvalue adj_pvalue
#1  ENSG00000241860.7 15.9176717  4 0.003131649 0.03386405
#2 ENSG00000237094.12  7.7036345  4 0.103057939 0.30914375
#3 ENSG00000237491.10 18.9855671  7 0.008232674 0.06409992
#4  ENSG00000230092.7  0.2030046  1 0.652306490 0.82104369
#5 ENSG00000228794.10 18.3125611  5 0.002579141 0.03007630
#6 ENSG00000188976.11  0.2613558  2 0.877500363 0.94361403

```
```
res.txp <- DRIMSeq::results(d, level="feature")
head(res.txp)
```
Output
```
#    gene_id        feature_id          lr df       pvalue adj_pvalue
#1  ENSG00000241860.7 ENST00000484859.1  0.10925423  1 0.7409947831 0.96718722
#2  ENSG00000241860.7 ENST00000466557.6  0.09413004  1 0.7589908896 0.97209538
#3  ENSG00000241860.7 ENST00000662089.1  4.43238170  1 0.0352632682 0.29019308
#4  ENSG00000241860.7 ENST00000491962.1 13.73742817  1 0.0002102231 0.01361563
#5  ENSG00000241860.7 ENST00000655252.1  1.35295757  1 0.2447617350 0.69140984
#6 ENSG00000237094.12 ENST00000455464.7  0.03503514  1 0.8515220399 1.00000000
```
Check how many na's - these can be useful later as they can represent complete isoform switches
```
table(is.na(res$pvalue))

#FALSE  TRUE 
#10262     4 
 nas<-res[is.na(res$pvalue),]
 nas
#                gene_id lr df pvalue adj_pvalue
#210  ENSG00000060656.20 NA  6     NA         NA
#6627 ENSG00000231607.13 NA  5     NA         NA
#6695 ENSG00000134871.19 NA  7     NA         NA
#7819 ENSG00000102984.15 NA  7     NA         NA
write.table(nas,"Target_NAs_drmseqRes.txt",sep="\t",quote=FALSE)
table(is.na(res.txp$pvalue))

#FALSE  TRUE 
#40908    29 
 nas<-res.txp[is.na(res.txp$pvalue),]
 write.table(nas,"Target_NAs_drmseqRes.Txp.txt",sep="\t",quote=FALSE)
 dim(res)
#[1] 10266     5
 dim(res.txp)
#[1] 40937     6
```
```
no.na <- function(x) ifelse(is.na(x), 1, x)
res$pvalue <- no.na(res$pvalue)
res.txp$pvalue <- no.na(res.txp$pvalue)
idx <- which(res$adj_pvalue < 0.05)
length(idx)
#[1] 1157

idx <- which(res$adj_pvalue < 0.05)[1]
res[idx,]
#             gene_id       lr df      pvalue adj_pvalue
#5 ENSG00000228794.10 20.77846  7 0.004112255 0.04317466
plotProportions(d, res$gene_id[idx], "condition")

plotProportions(d, "ENSG00000165280.18", "condition",plot_type = "boxplot1")
```
### Order and save results
```
res.txp.ordered <- res.txp[order(res.txp$adj_pvalue),]
write.table(res.txp.ordered,"~/Documents/swimmongDownstream/DrimseqRes_res.txp.ordered.txt",sep="\t",quote=FALSE)
res.ordered <- res[order(res$adj_pvalue),]
write.table(res.ordered,"~/Documents/swimmongDownstream/DrimseqRes_res.ordered.txt",sep="\t",quote=FALSE)

```

# stageR following DRIMSeq

```
#construct a vector of pValues
pScreen <- res$pvalue
strp <- function(x) substr(x,1,15)
names(pScreen) <- strp(res$gene_id)

#We construct a one column matrix of the confirmation p-values:

pConfirmation <- matrix(res.txp$pvalue, ncol=1)
rownames(pConfirmation) <- strp(res.txp$feature_id)
#We arrange a two column data.frame with the transcript and gene identifiers.

tx2gene <- res.txp[,c("feature_id", "gene_id")]
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])
```
## load and run StageR
```
library(stageR)
stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=FALSE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
  drim.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                  onlySignificantGenes=TRUE)
})

#The returned adjusted p-values are based on a stage-wise testing approach and are only valid for the provided target OFDR level of 5%. If a different #target OFDR level is of interest,the entire adjustment should be re-run. 

head(drim.padj)
```
output:
```
           geneID            txID       gene transcript
1 ENSG00000241860 ENST00000484859 0.03387725 1.000000000
2 ENSG00000241860 ENST00000466557 0.03387725 1.000000000
3 ENSG00000241860 ENST00000662089 0.03387725 0.939479356
4 ENSG00000241860 ENST00000491962 0.03387725 0.005600737
5 ENSG00000241860 ENST00000655252 0.03387725 1.000000000
6 ENSG00000228794 ENST00000445118 0.03008803 1.000000000
```
```
write.table(drim.padj,"drim_padj.txt",sep="\t",quote=FALSE)
```
### extract significant genes and transcripts
```
sigGenes<-getSignificantGenes(stageRObj)
sigTranscripts<-getSignificantTx(stageRObj)
write.table(sigGenes,"~/Documents/swimmongDownstream/significantGenesStageR.txt",sep="\t",quote=FALSE)
write.table(sigTranscripts,"~/Documents/swimmongDownstream/significantTxStageR.txt",sep="\t",quote=FALSE)
dim(sigGenes)
#[1] 1156    1
 dim(sigTranscripts)
#[1] 1394    1

```

## Post-hoc filtering on the standard deviation in proportions
```
res.txp.filt <- DRIMSeq::results(d, level="feature")
smallProportionSD <- function(d, filter=0.1) {
  cts <- as.matrix(subset(counts(d), select=-c(gene_id, feature_id)))
  gene.cts <- rowsum(cts, counts(d)$gene_id)
  total.cts <- gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
  props <- cts/total.cts
  propSD <- sqrt(rowVars(props))
  propSD < filter
}
filt <- smallProportionSD(d)
res.txp.filt$pvalue[filt] <- 1 
res.txp.filt$adj_pvalue[filt] <- 1

idx <- which(res.txp.filt$adj_pvalue < 0.05)
length(idx)
#[1] 487

 head(res.txp.filt)
 ```
 Output:
 ```
          gene_id        feature_id          lr df       pvalue adj_pvalue
1  ENSG00000241860.7 ENST00000484859.1  0.10925423  1 1.0000000000 1.00000000
2  ENSG00000241860.7 ENST00000466557.6  0.09413004  1 1.0000000000 1.00000000
3  ENSG00000241860.7 ENST00000662089.1  4.43238170  1 0.0352632682 0.29019308
4  ENSG00000241860.7 ENST00000491962.1 13.73742817  1 0.0002102231 0.01361563
5  ENSG00000241860.7 ENST00000655252.1  1.35295757  1 0.2447617350 0.69140984
6 ENSG00000237094.12 ENST00000455464.7  0.03503514  1 1.0000000000 1.00000000
```
```
write.table(res.txp.filt,"res.txp.filt.txt",sep="\t",quote=FALSE)
```

## DEXSeq - an alternative, starting from DRIM
```
library(DEXSeq)
sample.data <- DRIMSeq::samples(d)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
head(sample.data)

sample.data$condition<-as.factor(sample.data$condition)
dxd <- DEXSeqDataSet(countData=count.data,
sampleData=sample.data,
design=~sample + exon + condition:exon,
featureID=counts(d)$feature_id,
groupID=counts(d)$gene_id)

system.time({
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, quiet=TRUE)
save.image(file="afterEstimateDispersionTarget.RData")
dxd <- testForDEU(dxd, reducedModel=~sample + exon)
save.image(file="afterTestForDEUTarget.RData")
})

```
### Look at DEXSeq results
```
#if you are starting from saved results
#load("afterTestForDEU.RData")
dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
qval <- perGeneQValue(dxr)
dxr.g <- data.frame(gene=names(qval),qval)
columns <- c("featureID","groupID","pvalue")
dxr <- as.data.frame(dxr[,columns])
head(dxr)
```
### Output:
```
#                                             featureID            groupID
#ENSG00000241860.7:ENST00000484859.1  ENST00000484859.1  ENSG00000241860.7
#ENSG00000241860.7:ENST00000466557.6  ENST00000466557.6  ENSG00000241860.7
#ENSG00000241860.7:ENST00000662089.1  ENST00000662089.1  ENSG00000241860.7
#ENSG00000241860.7:ENST00000491962.1  ENST00000491962.1  ENSG00000241860.7
#ENSG00000241860.7:ENST00000655252.1  ENST00000655252.1  ENSG00000241860.7
#ENSG00000237094.12:ENST00000455464.7 ENST00000455464.7 ENSG00000237094.12
#                                         pvalue
#ENSG00000241860.7:ENST00000484859.1  0.78245764
#ENSG00000241860.7:ENST00000466557.6  0.01831049
#ENSG00000241860.7:ENST00000662089.1  0.20590015
#ENSG00000241860.7:ENST00000491962.1  0.94145980
#ENSG00000241860.7:ENST00000655252.1  0.49742273
#ENSG00000237094.12:ENST00000455464.7 0.04361211
```
```
write.table(dxr,"DEXseq_result_dxr",sep="\t",quote=FALSE)
library(stageR)
strp <- function(x) substr(x,1,15)
pConfirmation <- matrix(dxr$pvalue,ncol=1)
dimnames(pConfirmation) <- list(strp(dxr$featureID),"transcript")
pScreen <- qval
names(pScreen) <- strp(names(pScreen))
tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=TRUE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
  dex.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                 onlySignificantGenes=TRUE)
})

head(dex.padj)
```
### Output:
```
#           geneID            txID        gene transcript
#1 ENSG00000228794 ENST00000445118 0.005537052 0.00150867
#2 ENSG00000228794 ENST00000669922 0.005537052 0.35596946
#3 ENSG00000228794 ENST00000667414 0.005537052 0.11874327
#4 ENSG00000228794 ENST00000658846 0.005537052 1.00000000
#5 ENSG00000228794 ENST00000657837 0.005537052 1.00000000
#6 ENSG00000188976 ENST00000327044 0.005323858 0.67044126
```
```
dim(dex.padj)
#[1] 12989     4
 write.table(dex.padj,"DEXseq_Rstage_result_dex.padj",sep="\t",quote=FALSE)
```
