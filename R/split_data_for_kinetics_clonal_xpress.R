library(data.table)

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

dat.normumi.xpress <- readRDS("data/dat_normumi_xpress.rds")

# add wl read fraction and allelic tpm
dat.normumi.xpress[melt(fread("out_mRNA/clonal_xpress/zUMIs_output/allelic/clonal_xpress.fract_CAST_reads.txt"), id.vars = "GeneID"), readfrac := value, on = c("GeneID == GeneID", "Var2 == variable") ]
dat.normumi.xpress[, c("tpm.c57", "tpm.cast") := list((tpm/sizeFactor)*(1-readfrac), (tpm/sizeFactor)*readfrac)]

# correct values
ls.normumi.allelic.corr.xpress <- list(
  "c57" = correct.values(data.frame(dcast(dat.normumi.xpress, GeneID~Var2, value.var = "c57"), row.names=1), data.frame(dcast(dat.normumi.xpress, GeneID~Var2, value.var = "normumi"), row.names=1)), 
  "cast" = correct.values(data.frame(dcast(dat.normumi.xpress, GeneID~Var2, value.var = "cast"), row.names = 1), data.frame(dcast(dat.normumi.xpress, GeneID~Var2, value.var = "normumi"), row.names=1))
)

ls.tpm.allelic.corr.xpress <- list(
  "c57" = correct.values(data.frame(dcast(dat.normumi.xpress, GeneID~Var2, value.var = "tpm.c57"), row.names=1), data.frame(dcast(dat.normumi.xpress, GeneID~Var2, value.var = "tpm"), row.names=1)), 
  "cast" = correct.values(data.frame(dcast(dat.normumi.xpress, GeneID~Var2, value.var = "tpm.cast"), row.names = 1), data.frame(dcast(dat.normumi.xpress, GeneID~Var2, value.var = "tpm"), row.names=1))
)

spl.idx.xpress <- factor(sce.xpress.filt$clone)

# split data and export
sapply(levels(spl.idx.xpress),function(x){
  write.table(ls.normumi.allelic.corr.xpress$c57[,spl.idx.xpress == x],paste0("kinetics/input_xpress/normumi/xpress.normumi.C57.",sub("_","",x),".csv"),quote = F,sep = ",",col.names = NA,row.names = T)
  write.table(ls.normumi.allelic.corr.xpress$cast[,spl.idx.xpress == x],paste0("kinetics/input_xpress/normumi/xpress.normumi.CAST.",sub("_","",x),".csv"),quote = F,sep = ",",col.names = NA,row.names = T)
})

sapply(levels(spl.idx.xpress),function(x){
  write.table(ls.tpm.allelic.corr.xpress$c57[,spl.idx.xpress == x],paste0("kinetics/input_xpress/tpm/xpress.tpm.C57.",sub("_","",x),".csv"),quote = F,sep = ",",col.names = NA,row.names = T)
  write.table(ls.tpm.allelic.corr.xpress$cast[,spl.idx.xpress == x],paste0("kinetics/input_xpress/tpm/xpress.tpm.CAST.",sub("_","",x),".csv"),quote = F,sep = ",",col.names = NA,row.names = T)
})
