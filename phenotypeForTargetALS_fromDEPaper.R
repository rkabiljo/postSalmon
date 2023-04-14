cortexTarget<-read.table("~/Documents/SalmonTarget/samples.design.updated.site.sv1.sizefactors234.txt",header=TRUE,row.names = 1,sep="\t")
dim(cortexTarget)
samples2 <- read.table("~/Documents/SalmonTarget/TargetALSRNAMetadata09092021_1270.txt",sep="\t",header=TRUE,row.names=1)
#coldata2 <- data.frame(files, names=names(files), samples2)
pheno<-samples2[samples$sample_name,]
