setwd("~/Documents/R_Objects_DTU")

variantInfo="CFAP410_25selected.txt"
geneID="ENSG00000160226.16"

load("txiTranscriptsBB.RData")
source("plotExpressionLabFunc.R")
library(tximport)
library(tidyverse)
library(GenomicFeatures)
library(ggbeeswarm)
library(DRIMSeq)
#this includes transcripts read, using scaledTPM, for 171 BB subjects
#read the accompanying phenotype table
pheno<-read.table("BB_pheno.171subjects.txt",sep="\t",header=TRUE)
pheno<-pheno[colnames(txiTranscriptsBB$counts),]
#text file linking genes and transcripts
tx2gene<-read.table("tx2gene.txt",header=TRUE,sep="\t")
#head(tx2gene)

##summarise to genes, and save for later
#library(tximport)
#txiGenesBB <- summarizeToGene(txiTranscriptsBB, tx2gene[,c(2,1)])
#colnames(txiGenesBB$counts)<-colnames(txiTranscriptsBB$counts)
#save(txiGenesBB,file="txiGenesBB.RData")

txdb<-loadDb("~/Documents/swimmongDownstream/gencode.v38.hr_patch_hapl_scaff.annotation.sqlite")
txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
cts <- txiTranscriptsBB$counts
#dim(cts)
cts <- cts[rowSums(cts) > 0,]
#dim(cts)
#cts[1:3,1:3]
#range(colSums(cts)/1e6)
all(rownames(cts) %in% txdf$TXNAME)

txdf <- txdf[match(rownames(cts),txdf$TXNAME),]
#all(rownames(cts) == txdf$TXNAME)

#this is the place to select subjects/genes. the following line is for the complete set
#counts <- data.frame(gene_id=txdf$GENEID,
#                     feature_id=txdf$TXNAME,
#                     cts)

#head(counts[1:5,1:5])
#select all transcripts for a gene of interest - look up the exact version

tx2gene[tx2gene$GENEID==geneID,"TXNAME"]
#tx2gene[tx2gene$GENEID=="ENSG00000160226.16","TXNAME"]
#[1] "ENST00000339818.9" "ENST00000470196.5" "ENST00000496321.5" "ENST00000462742.1"
#[5] "ENST00000397956.7" "ENST00000325223.7" "ENST00000478674.1"

#filter cts, only for genes of interest
cts<-cts[tx2gene[tx2gene$GENEID=="ENSG00000160226.16","TXNAME"],]
dim(cts)
#[1]   7 171

#to filter also subjects, load a separate phenotype data with only these subjects you want to keep

#selected<-read.table("CFAP410_25selected.txt",sep="\t",header=TRUE)
selected<-read.table(variantInfo,sep="\t",header=TRUE)

phe<-selected
#table(phe$carrier)
#1  2 
#20  5

#make sure the cts matrix contains only this subset of selected subjects
cts<-cts[,rownames(selected)]

samps<-phe
samps$sample_id <- rownames(samps)
samps<-samps[,c("sample_id","carrier","Status")]
colnames(samps)<-c("sample_id","condition","Status")
samps$condition<-as.factor(samps$condition)
head(samps)

all(rownames(cts) %in% txdf$TXNAME)
txdf <- txdf[match(rownames(cts),txdf$TXNAME),]
all(rownames(cts) == txdf$TXNAME)
counts <- data.frame(gene_id=txdf$GENEID,
                     feature_id=txdf$TXNAME,
                     cts)


d <- dmDSdata(counts=counts, samples=samps)
d

#optional filtering
n = nrow(samps)
n.small = min(table(samps$condition))

d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
              min_samps_gene_expr=n, min_gene_expr=10)
d

design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
colnames(design_full)
system.time({
  d <- dmPrecision(d, design=design_full)
  d <- dmFit(d, design=design_full)
  d <- dmTest(d, coef="condition2")
})



#the following is needed for the more sophisticated plotting function, calculating proportions
drim.prop = reshape2::melt(counts[counts$feature_id %in% proportions(d)$feature_id,], id = c("gene_id", "feature_id"))
drim.prop = drim.prop[order(drim.prop$gene_id, drim.prop$variable, drim.prop$feature_id),]

# Calculate proportions from counts
system.time({
  drim.prop = drim.prop %>%
    group_by(gene_id, variable) %>%
    mutate(total = sum(value)) %>%
    group_by(variable, add=TRUE) %>%
    mutate(prop = value/total)
})
drim.prop = reshape2::dcast(drim.prop[,c(1,2,3,6)], gene_id + feature_id ~ variable)

plotExpressionLab(drim.prop, "ENSG00000160226.16", samps, isProportion = TRUE, labelVar= "Status")

