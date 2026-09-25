library(data.table)

##

outlier.mad <- function(x, nmad=2){
  x > median(x, na.rm=T)+(nmad*mad(x,na.rm=T))
}

##

## Load genome metrics
genome.nuc <- fread("ucsc/nucfreq_100kb.tsv.gz")

genome.map <- fread("ucsc/mappability_100kb.bedgraph.gz")
genome.map[, V4 := as.numeric(V4)][is.na(V4), V4 := 0]

genome.gap <- fread("ucsc/gap_100kb.bedgraph.gz")

genome.rmsk <- fread("ucsc/rmsk_100kb.bedgraph.gz")

# genome.ideo <- fread("ucsc/cytoBand.txt.gz")
# genome.ideo[, V5 := factor(V5, levels=c("gneg","gpos33","gpos66", "gpos75","gpos100"))]

## Load data
# binned read depth
cnt.depth <- melt(fread("out_dna/cnv/clonal_DNAseq/read.counts.100kb.bin.tsv.gz"), id.vars = c("CHROM", "START", "END"))

colnames(cnt.depth) <- c("chr", "start", "end", "variable", "reads")
cnt.depth[, variable := gsub(".*/|\\..*", "", variable)]
cnt.depth[, chr := factor(chr, levels = c(1:19, "X", "Y"))]

cnt.depth[genome.nuc, c("gc","n") := list(`5_pct_gc`, `10_num_N`/`12_seq_len`), on = c("chr == `#1_usercol`", "start == 2_usercol")]
cnt.depth[genome.map, map := as.numeric(V4), on = c("chr == V1", "start == V2")]
cnt.depth[genome.gap, gap := V4, on = c("chr == V1", "start == V2")]
cnt.depth[genome.rmsk, rmsk := V4, on = c("chr == V1", "start == V2")]
cnt.depth[, exclude := n>0.025|map<0.55|gap>0|rmsk>0.7|is.na(chr)]
cnt.depth[frollapply(exclude, 5, median, align="center") == 1, exclude := TRUE, by=c("variable", "chr")] # smoothen excluded regions

#fwrite(cnt.depth, "out_dna/cnv/clonal_DNAseq/dat_clonal_depth.tsv", quote = F, sep = "\t")

## identify genomic copy segments with HMMcopy
library(HMMcopy)
ls.depth.cor <- lapply(split(cnt.depth[exclude == F], cnt.depth[exclude == F, variable]), correctReadcount, mappability = 0.55, samplesize = 5e3, verbose=F)
cnt.depth.cor <- rbindlist(ls.depth.cor)
# normalise to control
cnt.depth.cor[, copy.raw := copy]
cnt.depth.cor[, copy := copy - copy[variable == "mTTF_female"], by=c("chr", "start", "end")]

ls.depth.seg <- lapply(ls.depth.cor, function(x){
  param <- HMMsegment(x, getparam = T)
  param$e <- 0.9999999999999999
  param$strength <- 1e30
  HMMsegment(x, param, verbose=F)
})
cnt.depth.seg <- rbindlist(lapply(ls.depth.seg, function(x) as.data.table(x$segs) ), idcol = "variable")
# merge segments with a smaller absolute diff than 0.2
cnt.depth.seg[, seg.grp := cumsum(c(TRUE, abs(diff(median)) > 0.2)), by=c("variable", "chr")]

cnt.depth.seg.red <- cnt.depth.seg[, list(start=min(start), end=max(end), state=state[which.max(end-start)], median=weighted.mean(median, end-start)), by=c("variable", "chr", "seg.grp")]
cnt.depth.seg.red[, label := c("HOMD", "HETD", "NEUT", "GAIN", "AMPL", "HLAMP")[as.numeric(state)]]
cnt.depth.seg.red[, segment := seq_along(state), by=c("variable", "chr")]

#fwrite(cnt.depth.seg, "out_dna/cnv/clonal_DNAseq/genome_segments.tsv", quote = F, sep = "\t")

# add segment info
setkey(cnt.depth.seg.red, variable, chr, start, end)
cnt.depth.cor[, c("seg_end", "label", "seg_median") := cnt.depth.seg.red[foverlaps(cnt.depth.cor, cnt.depth.seg.red, by.x = c("variable","chr","start", "end"), mult="first", which = T), list(end, label, median)]]

## ploidy estimation using ACE

# sample ploidy cellularity
# "C1_B8" 2.00 1.00 # diploid, 3rd best fit
# "C1_D10" 2.00 1.00 # diploid
# "C2_B8" 3.96 1.00 # tetra
# "C2_F9" 4.00 1.00 # tetra
# "C2_H5" 3.96 1.00 # tetra
# "C3_C8" 3.96 1.00 # tetra
# "C3_D9" 2.00 1.00 # diploid
# "C3_F1_repeat" 4.08 0.99 # tetra
# "C3_F3" 2.00 1.00 # diploid
# "C3_F4" 3.96 1.00 # tetra
# "C3_G8" 2.00 1.00 # diploid
# "M3_G12" 3.96 1.00 # tetra 
# "M5_C3" 4.00 1.00 # tetra 2nd best fit
# "M8_C10" 4.00 1.00 # tetra
# "M8_D9" 4.00 1.00 # tetra 3rd best fit
# "M9_B10" 3.84 1.00 # tetra
# "M9_B3" 3.96 1.00 # tetra
# "M9_D5" 3.96 0.94 # tetra
# "M9_E3" 2.00 1.00 # diploid
# "MEF_Bm4" 2.00 1.00 # diploid
# "mTTF_female" 2.00 1.00 # diploid
# "new_C2_B8" 3.96 1.00 # tetra

##

cnt.depth.ace <- rbindlist(lapply(split(cnt.depth.cor, cnt.depth.cor$variable), function(x){
  require(ACE)
  dat <- x[ideal == T, list(bin = seq_along(chr), chr, start, end, copynumbers=2^(copy-median(copy[label=="NEUT"],na.rm=T)), segments=2^(seg_median-median(seg_median[label=="NEUT"], na.rm=T)) )] 
  fit <- squaremodel(dat, penalty = 1, penploidy = 0.5, method = "MAE")
  
  fit.min <- fit$minimadf[ifelse(unique(x$variable) %in% c("C1_B8", "M8_D9"),3,1),]
  message(unique(x$variable))
  #print(head(fit$minimadf))
  print(fit.min)
  
  dat[, c("ploidy","cellularity") := list(round(fit.min[,"ploidy"]), fit.min[,"cellularity"])]
  seg <- getadjustedsegments(dat, cellularity = fit.min[,'cellularity'], ploidy = fit.min[,'ploidy'])
  
  dat[!is.na(copynumbers), copies := with(seg, rep(Copies, Num_Bins))]
  dat[, adjustedcopynumbers := copynumbers*(ploidy+2/cellularity-2)/median(segments,na.rm=T) - ifelse(chr == "Y", 1,2)/cellularity + ifelse(chr == "Y", 1,2)]
  dat[, adjustedsegments := segments*(ploidy+2/cellularity-2)/median(segments,na.rm=T) - ifelse(chr == "Y", 1,2)/cellularity + ifelse(chr == "Y", 1,2)]
  return(dat)
}), idcol = "variable")

cnt.depth.seg.red[cnt.depth.ace, ploidy := round(ploidy), on = c("variable == variable")] # integer ploidy

# add segment info
setkey(cnt.depth.ace, variable, chr, start, end)
cnt.depth.ace[, segment := cnt.depth.seg.red[foverlaps(cnt.depth.ace, cnt.depth.seg.red, mult="first", which = T), segment]]

## allele frequencies
cnt.baf <- melt(fread("out_dna/cnv/clonal_DNAseq/variants_cast_matched_counts.tsv.gz", dec = "."), id.vars = c("CHROM", "POS", "REF", "ALT"))
cnt.baf[, variable := gsub(".*/|\\..*", "", variable)]
cnt.baf[, CHROM := factor(CHROM, levels = c(1:19, "X", "Y"))]
cnt.baf[, c("count_ref", "count_alt") := tstrsplit(value, ",")]
cnt.baf[, count_ref := as.numeric(count_ref)]
cnt.baf[, count_alt := as.numeric(count_alt)]
cnt.baf[, count_total := count_ref + count_alt]
cnt.baf <- cnt.baf[count_total>0]
cnt.baf[, c("POS_100kb") := list(round(POS/1e5)*1e5)] # NOTE: for whatever reason round(POS, -5) is incorrect at the 10th decimal for some values
cnt.baf[cnt.depth, exclude := exclude, on = c("CHROM == chr", "POS_100kb == start", "variable == variable")]

# summarise for 100kb regions
cnt.baf.sum <- cnt.baf[exclude == F & count_total >= 3, list("count_ref"=sum(count_ref, na.rm=T), "count_total" = sum(count_total, na.rm=T)), by=c("CHROM", "POS_100kb", "variable")]
cnt.baf.sum[, ratio := count_ref/count_total]

# add segment info
cnt.baf.sum[, temp := POS_100kb] # dummy variable
cnt.baf.sum[, c("seg_end", "segment", "label") := cnt.depth.seg.red[foverlaps(cnt.baf.sum, cnt.depth.seg.red, by.x = c("variable","CHROM","POS_100kb", "temp"), mult="first", which = T), list(end, segment, label )]]
cnt.baf.sum[cnt.depth.seg.red, ploidy := round(ploidy), on=c("variable == variable")]
#cnt.baf.sum[, ratio_bg := sum(count_ref[label=="NEUT"])/(sum(count_total[label=="NEUT"])), by="variable"]

cnt.baf.seg <- cnt.baf.sum[count_total >= 10, list("start" = min(POS_100kb), "ratio" = median(ratio), "ratio_bg" = median(ratio_bg)), by=c("variable", "ploidy", "CHROM", "segment")]

# summarise
cnt.depth.seg.red[cnt.baf.seg, ratio := ratio, on = c("variable == variable", "chr == CHROM", "segment == segment")]
cnt.depth.ace[cnt.baf.sum, ratio := ratio, on = c("variable == variable", "chr == CHROM", "start == POS_100kb")]

## export for manual checking
fwrite(cnt.depth.seg.red[, list(variable, chr, start, end, median, segment, ploidy, "ratio" = round(ratio,2) )], "data/dat_clonal_cnv.tsv", quote = F, sep = "\t")
fwrite(cnt.depth.ace[, list(variable, chr, start, end, segment, adjustedcopynumbers, adjustedsegments, ratio)], "data/dat_clonal_dnaseq.tsv", quote = F, sep = "\t")

##