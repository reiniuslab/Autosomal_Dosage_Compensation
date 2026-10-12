library(data.table)
library(magrittr)

##

rpkm2tpm <- function(x){
  t(t(x) / colSums(x, na.rm=T) * 1e6)
}

scale.allelic <- function(x, mat_frac, alleles = c("C57", "CAST")){
  require(data.table)
  mat_frac2 <- mat_frac[apply(mat_frac,1,var, na.rm=T)>0,] # filter invariant genes
  
  idx.col <- intersect(colnames(x), colnames(mat_frac2))
  idx.row <- intersect(row.names(x), row.names(mat_frac2))
  
  out <- list(
    as.data.table(x[idx.row,idx.col] * (1 - mat_frac2[idx.row, idx.col]), keep.rownames=T),
    as.data.table(x[idx.row,idx.col] * mat_frac2[idx.row, idx.col], keep.rownames=T)
  )
  names(out) <- alleles
  return(out)
}

##

## Load data
# annotations
library(rtracklayer)
anno.xpress <- list(
  "barcode" = fread("out_mRNA/clonal_xpress_1/barcode_annotation_1.tsv.gz"),
  "barcode2" = fread("out_mRNA/clonal_xpress_2/barcode_annotation_2.tsv.gz"),
  "gene" = import("out_mRNA/clonal_xpress_1/clonal_xpress_1.final_annot.gtf.gz", feature.type = "gene")
)

# expression
dge.xpress <- readRDS("out_mRNA/clonal_xpress_1/zUMIs_output/expression/clonal_xpress.dgecounts.rds")
ls.dge.xpress <- list(
  "umi_c57" = fread("out_mRNA/clonal_xpress_1/zUMIs_output/allelic/clonal_xpress_1.BL6_direct_UMIs.txt"),
  "umi_cast" = fread("out_mRNA/clonal_xpress_1/zUMIs_output/allelic/clonal_xpress_1.CAST_direct_UMIs.txt"),
  "read_c57" = fread("out_mRNA/clonal_xpress_1/zUMIs_output/allelic/clonal_xpress_1.BL6_reads.txt"),
  "read_cast" = fread("out_mRNA/clonal_xpress_1/zUMIs_output/allelic/clonal_xpress_1.CAST_reads.txt"),
  "molspike" = fread("out_mRNA/clonal_xpress_1/zUMIs_output/molspikes/clonal_xpress_1.spikes.filtered.txt.gz")
)

dge.xpress2 <- readRDS("out_mRNA/clonal_xpress_2/zUMIs_output/expression/clonal_xpress_2.dgecounts.rds")
ls.dge.xpress2 <- list(
  "umi_c57" = fread("out_mRNA/clonal_xpress_2/zUMIs_output/allelic/clonal_xpress_2.BL6_direct_UMIs.txt"),
  "umi_cast" = fread("out_mRNA/clonal_xpress_2/zUMIs_output/allelic/clonal_xpress_2.CAST_direct_UMIs.txt"),
  "read_c57" = fread("out_mRNA/clonal_xpress_2/zUMIs_output/allelic/clonal_xpress_2.BL6_reads.txt"),
  "read_cast" = fread("out_mRNA/clonal_xpress_2/zUMIs_output/allelic/clonal_xpress_2.CAST_reads.txt"),
  "molspike" = fread("out_mRNA/clonal_xpress_2/zUMIs_output/molspikes/clonal_xpress_2.spikes.filtered.txt.gz")
)

## Spike-in normalisation
library(scater)
library(scran)
library(magrittr)
# batch1
sce.xpress <- SingleCellExperiment(assays = list(counts = dge.xpress$umicount$exon$all[!grepl("molspike",row.names(dge.xpress$umicount$exon$all)),], rpkm = dge.xpress$umicount$exon$all[!grepl("molspike",row.names(dge.xpress$rpkm$exon$all)),]), colData=anno.xpress$barcode[match(colnames(dge.xpress$umicount$exon$all), barcode),])
sce.xpress %<>% quickPerCellQC()
altExp(sce.xpress, "spumi") <- SummarizedExperiment(list(counts=as.matrix(ls.dge.xpress$molspike,rownames=1)[, colnames(sce.xpress)]))
sce.xpress$spike.pct <- colSums(assay(altExp(sce.xpress))) / (colSums(assay(altExp(sce.xpress)))+colSums(assay(sce.xpress)))
sce.xpress <- subset(sce.xpress,,spike.pct<0.1)
sce.xpress$outlier <- isOutlier(sce.xpress$detected, type="low", nmads = 3, log=T)
sce.xpress.filt <- subset(sce.xpress,,!outlier)
sce.xpress.filt %<>% computeSpikeFactors("spumi")
sce.xpress.filt %<>% logNormCounts()
sce.xpress.filt %<>% logNormCounts(transform = "none")

#batch2
sce.xpress2 <- SingleCellExperiment(assays = list(counts = dge.xpress2$umicount$exon$all[!grepl("molspike",row.names(dge.xpress2$umicount$exon$all)),], rpkm = dge.xpress2$umicount$exon$all[!grepl("molspike",row.names(dge.xpress2$rpkm$exon$all)),]), colData=anno.xpress$barcode2[match(colnames(dge.xpress2$umicount$exon$all), barcode),])
sce.xpress2 %<>% quickPerCellQC()
altExp(sce.xpress2, "spumi") <- SummarizedExperiment(list(counts=as.matrix(ls.dge.xpress2$molspike,rownames=1)[, colnames(sce.xpress2)]))
sce.xpress2$spike.pct <- colSums(assay(altExp(sce.xpress2))) / (colSums(assay(altExp(sce.xpress2)))+colSums(assay(sce.xpress2)))
sce.xpress2 <- subset(sce.xpress2,,spike.pct<0.1)
sce.xpress2$outlier <- isOutlier(sce.xpress2$detected, type="low", nmads = 3, log=T)
sce.xpress2.filt <- subset(sce.xpress2,,!outlier)
sce.xpress2.filt %<>% computeSpikeFactors("spumi")
sce.xpress2.filt %<>% logNormCounts()
sce.xpress2.filt %<>% logNormCounts(transform = "none")

# merge
sce.xpress.merge <- cbind(removeAltExps(sce.xpress.filt[intersect(rownames(sce.xpress.filt),rownames(sce.xpress2.filt))]),removeAltExps(sce.xpress2.filt[intersect(rownames(sce.xpress.filt),rownames(sce.xpress2.filt))]))
sce.xpress.merge$dataset <- as.character(as.numeric(colnames(sce.xpress.merge) %in% colnames(sce.xpress2.filt))+1)

## melt and combine UMI data
dat.normumi.xpress <- as.data.table(reshape2::melt(as.matrix(assay(sce.xpress.merge,"normcounts")[rowSums(assay(sce.xpress.merge, "normcounts"))>0,]), value.name = "normumi"))
colnames(dat.normumi.xpress)[1] <- "GeneID"
# add expression
dat.normumi.xpress[as.data.table(reshape2::melt(as.matrix(assay(sce.xpress.merge, "rpkm")[rowSums(assay(sce.xpress.merge, "rpkm"))>0,]))), rpkm := value, on = c("GeneID == Var1", "Var2 == Var2")]
dat.normumi.xpress[, tpm := rpkm / sum(rpkm,na.rm=T) * 10e6, by="Var2"]
dat.normumi.xpress[, c("tpm.avg", "gene.det") := list(mean(tpm, na.rm=T), mean(normumi>0)), by="GeneID"]
# add annotations
dat.normumi.xpress[as.data.table(anno.xpress$gene), c("chrom", "pos") := list(seqnames, start), on = c("GeneID == gene_id")]
dat.normumi.xpress[as.data.table(colData(sce.xpress.merge)), c("clone", "spike.pct", "sizeFactor", "dataset") := list(clone, spike.pct, sizeFactor, dataset), on = "Var2 == barcode"]
# add allelic expression
dat.normumi.xpress[rbind(melt(ls.dge.xpress$umi_c57, id.vars = "GeneID"),melt(ls.dge.xpress2$umi_c57, id.vars = "GeneID")), c57 := value/sizeFactor, on = c("GeneID == GeneID", "Var2 == variable")]
dat.normumi.xpress[rbind(melt(ls.dge.xpress$umi_cast, id.vars = "GeneID"),melt(ls.dge.xpress2$umi_cast, id.vars = "GeneID")), cast := value/sizeFactor, on = c("GeneID == GeneID", "Var2 == variable")]

## allelic tpm
dat.tpm.allelic.xpress <- as.data.table(reshape2::melt(scale.allelic(rpkm2tpm(data.frame(ls.dge.xpress$rpkm, row.names = 1)), data.frame(ls.dge.xpress$read_cast, row.names=1)/(data.frame(ls.dge.xpress$read_c57,row.names = 1)+data.frame(ls.dge.xpress$read_cast,row.names=1))), id.vars="rn"))
colnames(dat.tpm.allelic.xpress)[1] <- "GeneID"
dat.tpm.allelic.xpress[as.data.table(colData(sce.xpress.filt)), c("clone", "sizeFactor") := list(clone, sizeFactor), on = "variable == XC_DNBPE"]
dat.tpm.allelic.xpress[as.data.table(anno.xpress$gene), c("chrom", "gene_name") := list(seqnames, gene_name), on = c("GeneID == gene_id")]
dat.tpm.allelic.xpress[, tpm.avg := rowMeans(rpkm2tpm(assay(sce.xpress.filt, "rpkm")))[GeneID]]
dat.tpm.allelic.xpress %<>% .[!is.na(clone)]

## add DNA-seq segments
dat.clonal.cnv <- fread("data/dat_clonal_cnv_fixed.tsv")
dat.clonal.cnv[, variable := gsub("_repeat","",variable)]

# find overlaps using data.table
setkey(dat.clonal.cnv, variable, chr, start, end)

dat.normumi.xpress[, temp := pos] # dummy variable
dat.normumi.xpress[, c("ploidy", "segment", "copies", "ratio", "note") := dat.clonal.cnv[foverlaps(dat.normumi.xpress, dat.clonal.cnv, by.x = c("clone","chrom","pos", "temp"), mult="first", which = T), list(ploidy, segment, copies, ratio, note)]]
dat.normumi.xpress[, temp := NULL]

## EXPORT
saveRDS(sce.xpress, "data/sce_xpress1.rds")
saveRDS(sce.xpress.filt, "data/sce_xpress1_filt.rds")
saveRDS(sce.xpress2, "data/sce_xpress2.rds")
saveRDS(sce.xpress2.filt, "data/sce_xpress2_filt.rds")
saveRDS(sce.xpress.merge, "data/sce_xpress_merge.rds")

saveRDS(dat.normumi.xpress, "data/dat_normumi_xpress.rds")

write.table(assay(aggregateAcrossCells(sce.xpress.merge, ids = sce.xpress.merge$clone, statistics = "mean", use.assay.type = "normcounts")), "data/dat_normumi_pseudobulk_xpress.tsv.gz", sep = "\t", quote = F, col.names = NA)

## DONE
