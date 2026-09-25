suppressMessages(source("R/preprocess_data_fibroblast_ss2.R"))

##

correct.values <- function(x,all){
  # corrects allelic values depending on overall detection
  # if the gene is detected => 0, else missing
  x[which(is.na(x),arr.ind = T)] <- 0
  x[which(all == 0,arr.ind = T)] <- ""
  x <- x[apply(all,1,mean) != 0,]
  return(x)
}

##

# get subgroup
tpm.melt.avg.fibro <- tpm.melt.fibro[expressed == T & !is.na(chr), base::mean(value,na.rm=T),by=c("sample","clone_id", "sample_type","l1","chr", "sex")]
tpm.melt.avg.fibro[,refratio := V1[l1 == "c57"] / (V1[l1 == "c57"] + V1[l1 == "cast"]), by=c("sample", "chr")]
tpm.melt.avg.fibro[,subgroup := refratio %between% c(0.25,0.75)]

# x status
tpm.melt.avg.fibro[sample_type == "single.fibroblast", x.status := ifelse(refratio[chr == "chrX"] > 0.9, "XaXi", ifelse(refratio[chr == "chrX"] < 0.1, "XiXa", NA))]

# correct values and export split data
idx.row <- tpm.melt.all.fibro[!is.na(chr) & !duplicated(gene), as.character(gene)]
idx.col <- tpm.melt.all.fibro[sample_type == "single.fibroblast" & !duplicated(sample), as.character(sample)]

tpm.allele.corr.fibro <-  lapply(tpm.allele.fibro, function(x) correct.values(x[idx.row, idx.col], tpm.all.fibro[idx.row, idx.col]) )

spl.idx.fibro <- factor(tpm.melt.avg.fibro[chr == "chr3"][match(idx.col, sample), paste(ifelse(clone_id == "non.clonal","nonclonal","clonal"),ifelse(subgroup,"normal","loss"),sep=".")])

sapply(levels(spl.idx.fibro)[table(spl.idx.fibro)>20],function(x){
  write.table(tpm.allele.corr.fibro$c57[,spl.idx.fibro == x],paste0("kinetics/input/tpm.C57.",x,".csv"),quote = F,sep = ",",col.names = NA,row.names = T)
  write.table(tpm.allele.corr.fibro$cast[,spl.idx.fibro == x],paste0("kinetics/input/tpm.CAST.",x,".csv"),quote = F,sep = ",",col.names = NA,row.names = T)
})

##                                 
