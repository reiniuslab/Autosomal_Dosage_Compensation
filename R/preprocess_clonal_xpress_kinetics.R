library(data.table)

## Load data
# expression
dat.normumi.xpress <- readRDS("data/dat_normumi_xpress.rds")
dat.normumi.cloneavg.xpress <- dat.normumi.xpress[, lapply(.SD, mean, na.rm=T), by=c("GeneID", "clone"), .SDcols = c("normumi", "tpm", "c57", "cast")]
dat.normumi.cloneavg.xpress.melt <- melt(dat.normumi.cloneavg.xpress, measure.vars = c("c57", "cast"))
dat.normumi.cloneavg.xpress.melt[, variable := toupper(variable)]

# annotations
library(rtracklayer)
anno.xpress <- list(
  "barcode" = fread("out_mRNA/clonal_xpress/barcode_annotation.tsv.gz"),
  "gene" = import("out_mRNA/clonal_xpress/clonal_xpress.final_annot.gtf.gz", feature.type = "gene")
)

# burst kinetics
ls.burst.xpress <- list(
  "bf_tpm" = fread("kinetics/output_xpress/tpm/split/clonal_xpress_tpm_bf.tsv"),
  "bs_tpm" = fread("kinetics/output_xpress/tpm/split/clonal_xpress_tpm_bs.tsv"),
  "bf_normumi" = fread("kinetics/output_xpress/normumi/split/clonal_xpress_normumi_bf.tsv"),
  "bs_normumi" = fread("kinetics/output_xpress/normumi/split/clonal_xpress_normumi_bs.tsv")
)

dat.burst.xpress <- as.data.table(reshape2::melt(ls.burst.xpress, id.vars="gene"))
dat.burst.xpress[, metric := gsub("_.*", "", L1)]
dat.burst.xpress[, c("unit","allele", "library") := tstrsplit(variable, "_")]
dat.burst.xpress[, clone := sub("^(.{2})", "\\1_", library)]

dat.burst.xpress[as.data.table(anno.xpress$gene), c("name", "chrom", "pos", "width") := list(gene_name, seqnames, start, width), on = "gene == gene_id"]
dat.burst.xpress[unique(dat.normumi.xpress[,list(GeneID, tpm.avg)]), tpm.avg := tpm.avg, on = c("gene == GeneID")]
dat.burst.xpress[dat.normumi.cloneavg.xpress, normumi := normumi, on = c("gene == GeneID", "clone == clone")]
dat.burst.xpress[dat.normumi.cloneavg.xpress.melt, allelic.expression := i.value, on = c("gene == GeneID", "clone == clone", "allele == variable")]

# add DNA-seq segments
dat.clonal.cnv <- fread("data/dat_clonal_cnv_fixed.tsv")
dat.clonal.cnv[, variable := gsub("_repeat","",variable)]

setkey(dat.clonal.cnv, variable, chr, start, end)
dat.burst.xpress[, temp := pos] # dummy variable
dat.burst.xpress[, c("ploidy", "copies", "ratio", "segment", "note") := dat.clonal.cnv[foverlaps(dat.burst.xpress, dat.clonal.cnv, by.x = c("clone","chrom","pos", "temp"), mult="first", which = T), list(ploidy, copies, ratio, segment, note)]]
dat.burst.xpress[, temp := NULL]

# annotate alleles
dat.burst.xpress[, nalleles := ifelse(allele=="C57", round(copies*ratio), round(copies*(1-ratio)) )]
dat.burst.xpress[, typeallele := ifelse(nalleles >= (ploidy*0.5), "maj", "min")]
dat.burst.xpress[, alteration := c("loss", "neutral", "gain")[sign(copies-ploidy)+2]]

# relative to median
dat.burst.xpress[!is.na(segment), relative.burst := value/median(value,na.rm=T), by=c("variable","L1")]
dat.burst.xpress[!is.na(segment), relative.expression := allelic.expression/median(allelic.expression,na.rm=T), by=c("variable","L1")]

saveRDS(dat.burst.xpress, "data/dat_burst_xpress.rds")

##
