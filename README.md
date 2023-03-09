# postSalmon
analyse Salmon results

Sulev's annotation had gene ID appended to transcript ID, like this <br>

```
Name	Length	EffectiveLength	TPM	NumReads
ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|processed_transcript|	1657	1307.000	0.000000	0.000
```
and I needed it to be like this:
```
Name	Length	EffectiveLength	TPM	NumReads
ENST00000456328.2 1657 1356.000 0.000000 0.000
```
Run this to change it, just change the path:
```
 find /Users/renatakabiljo/Documents/MNDquant/pawsey0360/MNDsalmon/MNDquant/ -type f -name '*.sf' -exec sh -c 'awk -F "\t" "{gsub(/\|.*/, \"\", \$1)}1" OFS="\t" "$1" > tmp && mv tmp "$1"' _ {} \;  
```

or on Create

```
find /scratch/prj/alstargets/HG38Salmon/quant/ -type f -name '*.sf' -exec sh -c 'awk -F "\t" "{gsub(/\|.*/, \"\", \$1)}1" OFS="\t" "$1" > tmp && mv tmp "$1"' _ {} \;
```

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
columns(txdb)
keytypes(txdb)
k <- keys(txdb, keytype = "GENEID")
tx2gene <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
write.table(as.data.frame(tx2gene), file="tx2gene.txt", sep="\t")
```
## Get phenotype data for the motor cortex region
```
cortexTarget<-read.table("~/Documents/SalmonTarget/samples.design.updated.site.sv1.sizefactors234.txt",header=TRUE,row.names = 1,sep="\t")
dim(cortexTarget)
samples2 <- read.table("~/Documents/SalmonTarget/TargetALSRNAMetadata09092021_1270.txt",sep="\t",header=TRUE,row.names=1)
#coldata2 <- data.frame(files, names=names(files), samples2)
pheno<-samples2[samples$sample_name,]



## Read with tximport
setwd("~/Documents/MNDquantNew")
#specify which files to read -for example whatever is in my motor cortex matrix
files<-list.files(pattern = "quant.sf", recursive = TRUE)
samples3 <- data.frame(sample_name = gsub("\\.sf", "", files), file_name = files)
samples3$sample_name <- sapply(strsplit(samples3$file_name, "/"), head, n = 1)
notHere<-rownames(cortexTarget)[!rownames(cortexTarget) %in% samples3$sample_name]
here<-rownames(cortexTarget)[rownames(cortexTarget) %in% samples3$sample_name]
cortexTarget<-cortexTarget[rownames(cortexTarget) %in% here,]
phe<-cortexTarget[,c(1,2,3,4,5,6)]
samples3<-samples3[samples3$sample_name %in% here,]
files<-samples3$file_name
txiGenes <- tximport(files, type="salmon", tx2gene=tx2gene, txOut = FALSE)
dim(txiGenes$counts)
head(txiGenes$counts)
colnames(txiGenes$counts)<-samples3$sample_name


txiTranscripts <- tximport(files, type="salmon", txOut = TRUE)
head(txiTranscripts$counts)
colnames(txiTranscripts$counts)<-samples3$sample_name
head(txiTranscripts$counts)
```

