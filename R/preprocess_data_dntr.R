library(data.table)

##

correctReadcount_nomap <- function(x, mappability = 0.9, verbose = F, samplesize=50000) {
  # from HMMcopy but mappability normalisation removed
  if (length(x$reads) == 0 | length(x$gc) == 0 | length(x$map) == 0) {
    stop("Missing one of required columns: reads, gc, map")
  }
  
  if(verbose) { message("Applying filter on data...") }
  x$valid <- TRUE
  x$valid[x$reads <= 0 | x$gc < 0] <- FALSE
  x$ideal <- TRUE
  routlier <- 0.01
  range <- quantile(x$reads[x$valid], prob = c(0, 1 - routlier), na.rm = TRUE)
  doutlier <- 0.001
  domain <- quantile(x$gc[x$valid], prob = c(doutlier, 1 - doutlier),
                     na.rm = TRUE)
  x$ideal[!x$valid | x$map < mappability | x$reads <= range[1] |
          x$reads > range[2] | x$gc < domain[1] | x$gc > domain[2]] <- FALSE
  
  if (verbose) { message("Correcting for GC bias...") }
  set <- which(x$ideal)
  select <- sample(set, min(length(set), samplesize))
  rough = loess(x$reads[select] ~ x$gc[select], span = 0.03)
  i <- seq(0, 1, by = 0.001)
  final = loess(predict(rough, i) ~ i, span = 0.3)
  x$cor.gc <- x$reads / predict(final, x$gc)
  
  x$copy <- x$cor.gc
  x$copy[x$copy <= 0] = NA
  x$copy <- log(x$copy, 2)
  return(x)
}

counts2tpm <- function(x, len){
  rpk <- x/(len/1e3)
  sf <- sum(rpk, na.rm = T)/1e6
  rpk/sf
}

counts2rpgc <- function(x, len = 74, eff = 2747877777){
  sf <- (sum(x,na.rm = T) * len) / eff
  x / sf
}

##

## DNTR-seq data from Zachariadis et al. doi:10.1016/j.molcel.2020.09.025 (GSE144296)
# dna data
fls.dna <- list.files("GSE144296", "VZA.*bincounts", full.names = T)
names(fls.dna) <- gsub(".*GSM[[:digit:]]{7}_|-.*","", fls.dna)

fls.rna <- list.files("GSE144296", "VZA.*htseq", full.names = T)
names(fls.rna) <- gsub(".*GSM[[:digit:]]{7}_|.htseq.*","",gsub("-", ".", fls.rna, fixed = T))

meta <- fread("GSE144296/metadata.tsv")
meta[, sample := gsub(".", "_", cell_id, fixed = T)]
meta.filt <- meta[
  !is.na(rna_library_id) & 
  !is.na(dna_library_id) & 
  rna_library_id %in% names(fls.rna) &
  dna_library_id %in% names(fls.dna) & 
  dna_dups_pct < 0.5 &
  ercc_frac < 0.1 &
  sort_gate == "G1"
  ]

# rna data
gene.anno <- fread("GSE144296/geneExtents_noERCC_sorted.bed")

dat.rna <- rbindlist(lapply(fls.rna[meta.filt[,rna_library_id]], fread), idcol = T, use.names=F)
colnames(dat.rna) <- c("id", "gene", "counts")

library(scater)
library(scran)
library(magrittr)
mat.rna <- as.matrix(data.frame(dcast(dat.rna[!grepl("^__", gene)], gene~id, value.var = "counts"), row.names = 1))

sce.rna <- SingleCellExperiment(assay = list(counts = mat.rna[!grepl("^ERCC-", row.names(mat.rna)),]))
altExp(sce.rna, "ercc") <- SummarizedExperiment(list(counts=mat.rna[grepl("^ERCC-", row.names(mat.rna)),] ))
sce.rna %<>% quickPerCellQC()
sce.rna$outlier <- isOutlier(sce.rna$detected, type="low", nmads = 3, log=T)
sce.rna.filt <- subset(sce.rna,,!outlier)
sce.rna.filt %<>% computeSpikeFactors("ercc")
sce.rna.filt %<>% logNormCounts()
sce.rna.filt %<>% logNormCounts(transform = "none")

dat.rna <- data.table(reshape2::melt(assay(sce.rna.filt, "counts")))
colnames(dat.rna) <- c("gene", "id", "counts")
dat.rna[,id := gsub(".","_",id,fixed = T)]
dat.rna[,sample := gsub("R_","_",id, fixed = T)]
dat.rna[gene.anno, c("chr", "start", "end") := list(V1,V2,V3), on = "gene == V5"]
dat.rna[!is.na(chr), tpm := counts2tpm(counts, end-start), by="sample"]
dat.rna$normcounts <- reshape2::melt(assay(sce.rna.filt, "normcounts"))$value
dat.rna$logcounts <- reshape2::melt(assay(sce.rna.filt, "logcounts"))$value

dat.rna[meta.filt, c("treatment", "timepoint") := list(source, timepoint), on = c("sample == sample")]

# dna data

dna.anno <- fread("GSE144296/variable_500000_37_mappability_k50.bedgraph.gz")
colnames(dna.anno) <- c("chr", "start", "end", "map")
dna.anno$map <- as.numeric(dna.anno$map)
dna.anno$gc <- as.numeric(unlist(fread("GSE144296/variable_500000_37_nucfreq.tsv.gz")[,"5_pct_gc"]))

dat.dna <- rbindlist(lapply(meta.filt[,dna_library_id], function(x){
  library(HMMcopy)
  input <- fread(fls.dna[x])
  dat <- cbind(id = x,dna.anno, "reads" = as.numeric(unlist(input)) )
  
  # GC-correct readcounts and segment
  dat.cor <- correctReadcount_nomap(dat[gc > 0.05], mappability = 0.55, samplesize = 700)
  
  # modify parameters
  param <- HMMsegment(dat.cor, verbose = F, getparam = T)
  param$e <- 0.9999999999999999
  param$strength <- 1e30

  dat.seg <- as.data.table(HMMsegment(dat.cor, param, verbose = F)$segs)
  dat.seg[, seg := seq_along(start), by="chr"]
  
  # add segment data
  setkey(dat.cor, chr, start, end)
  setkey(dat.seg, chr, start, end)
  
  out <- foverlaps(dat.seg, dat.cor)[, list(chr, id, start, end, map, gc, reads, valid, ideal, cor.gc, copy, "seg_median" = median, seg)]
  return(out)
}), use.names = F)
dat.dna[, sample := gsub("D_", "_", id, fixed=T)]
dat.dna[, cov1x := counts2rpgc(reads) / (end-start), by="sample"]

dat.seg <- dat.dna[, list(start=min(start), end=max(end)), by=c("id", "sample", "chr", "seg", "seg_median")]

# combine
setkey(dat.rna, sample, chr, start, end)
setkey(dat.seg, sample, chr, start, end)

dat.rna <- dat.rna[foverlaps(dat.seg, dat.rna), c("seg", "seg_median") := list(seg, seg_median), on=c("sample == sample", "gene == gene")]

## EXPORT
saveRDS(dat.rna, "data/dntr_rna.rds")
saveRDS(dat.dna, "data/dntr_dna.rds")

##
