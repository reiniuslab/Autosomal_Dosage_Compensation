library(data.table)
library(magrittr)

##

scale.by.allele <- function(x, ref, alt, labs=c("c57","cast")){
  total <- ref+alt
  refratio <- ref/total
  x.ref <- x*refratio
  x.alt <- x*(1-refratio)
  x.allele <- list(x.ref, x.alt)
  names(x.allele) <- labs
  return(x.allele)
}

fpkm2tpm <- function(fpkm){
  (fpkm / sum(fpkm,na.rm = T)) * 1e6
}

add.info <- function(x, rowdata, coldata, rowjoin = "Var1==name2", coljoin = "Var2==LibraryName"){
  require(data.table)
  require(magrittr)
  
  x %<>% as.data.table()
  rowdata %<>% as.data.table()
  coldata %<>% as.data.table()
  
  x[rowdata, chr := chrom, on = rowjoin]
  x[chr != "chrY",chrx := chr == "chrX"]
  x %<>% .[coldata, on = coljoin]
  
  idx <- match(c(gsub("==.*","",rowjoin),gsub("==.*","",coljoin)),colnames(x))
  
  colnames(x)[idx] <- c("gene", "sample")
  colnames(x) %<>% tolower()
  return(x[!is.na(gene)])
}

##

## smart-seq2 data from Reinius et al. doi:10.1038/ng.3678 (GSE75659)
# load data, files from from https://github.com/RickardSandberg/Reinius_et_al_Nature_Genetics_2016
counts.c57.fibro <- readRDS("GSE75659/mouse.c57.counts.rds")
counts.cast.fibro <- readRDS("GSE75659/mouse.cast.counts.rds")
rpkm.all.fibro <- readRDS("GSE75659/mouse.rpkm.3digits.rds")

refseq.gene.fibro <- readRDS("GSE75659/mouse.gene.annotation.rds")

# Table S1a
meta.fibro <- fread("GSE75659/41588_2016_BFng3678_MOESM55_ESM_TableS1a.txt")
colnames(meta.fibro) <- gsub(" ", "_", colnames(meta.fibro))

# calculate reference ratios
refratio.fibro <- counts.c57.fibro / (counts.c57.fibro + counts.cast.fibro)
refratio.melt.fibro <- reshape2::melt(refratio.fibro)
refratio.melt.fibro %<>% add.info(refseq.gene.fibro, meta.fibro, coljoin = "Var2==Sample_id")
refratio.melt.fibro[, c("det", "var") := list(mean(!is.na(value)),var(value,na.rm=T)), by="gene"]

# scale rpkm matrix by allelic counts and convert to tpm
tpm.all.fibro <- apply(rpkm.all.fibro,2,fpkm2tpm)
tpm.allele.fibro <- scale.by.allele(tpm.all.fibro, counts.c57.fibro, counts.cast.fibro)
tpm.allele.fibro %<>% lapply(as.matrix)

# melt data
tpm.melt.all.fibro <- reshape2::melt(tpm.all.fibro)
tpm.melt.all.fibro %<>% add.info(refseq.gene.fibro, meta.fibro, coljoin = "Var2==Sample_id")

tpm.melt.fibro <- reshape2::melt(tpm.allele.fibro)
tpm.melt.fibro %<>% add.info(refseq.gene.fibro, meta.fibro, coljoin = "Var2==Sample_id")

# add expression annotations
tpm.melt.all.fibro[,expressed := mean(value,na.rm = T)>0, by = c("gene","sample_type")]
tpm.melt.fibro[,expressed := tpm.melt.all.fibro$expressed[match(paste0(gene,sample_type),with(tpm.melt.all.fibro,paste0(gene,sample_type)))]]
refratio.melt.fibro[,expressed := tpm.melt.all.fibro$expressed[match(paste0(gene,sample_type),with(tpm.melt.all.fibro,paste0(gene,sample_type)))]]

##
