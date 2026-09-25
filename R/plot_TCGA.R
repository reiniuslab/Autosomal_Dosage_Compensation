library(data.table)

##

median_cl_boot <- function(x, conf = 0.95) {
  # function from http://rstudio-pubs-static.s3.amazonaws.com/28101_41a7995107d94c8dbb07bbf7cd7e8291.html
  lconf <- (1 - conf)/2
  uconf <- 1 - lconf
  require(boot)
  bmedian <- function(x, ind) median(x[ind])
  bt <- boot(x, bmedian, 1000)
  bb <- boot.ci(bt, type = "perc")
  data.frame(y = median(x), ymin = quantile(bt$t, lconf), ymax = quantile(bt$t, uconf))
}

compare_slopes <- function(x,y,which=2){
  x.est <- coef(summary(x))[which, "Estimate"]
  y.est <- coef(summary(y))[which, "Estimate"]
  
  x.se <- coef(summary(x))[which, "Std. Error"]
  y.se <- coef(summary(y))[which, "Std. Error"]
  
  x.n <- nrow(x$model)
  y.n <- nrow(y$model)
  
  # Welch's T-test statistic
  t.stat <- (x.est - y.est) / sqrt(x.se^2 + y.se^2)
  
  # Degrees of freedom
  df <- ((x.se^2 + y.se^2)^2) / (((x.se^2)^2 / (x.n - 2)) + ((y.se^2)^2 / (y.n - 2)))
  
  # Calculate p-value
  p.value <- 2 * pt(-abs(t.stat), df)
  
  print(paste0(
    "x_slope=", round(x.est,3),
    " y_slope=", round(y.est,3),
    " p=", signif(p.value,3)
  ))
}

##

tcga.anno <- list(
  "sample" = fread("TCGA/gdc_sample_sheet.2021-01-28.tsv"),
  "gene" = fread("TCGA/gencode.gene.info.v22.tsv")
)
tcga.anno$sample[, sample_id := gsub("\\,.*","",`Sample ID`)]

# masked copy number
fls.tcga.masked <- list.files("TCGA/Masked_copy_number_segment", "*.seg.v2.txt", full.names = T, recursive = T)
names(fls.tcga.masked) <- sapply(strsplit(fls.tcga.masked,"/"), "[[", 3)

tcga.masked <- rbindlist(lapply(fls.tcga.masked, fread), idcol = T)
tcga.masked[tcga.anno$sample, c("sample_id", "sample_type") := list(sample_id, `Sample Type`), on = ".id == `File ID`"]

# gene copy number
fls.tcga.copy <- list.files("TCGA/Gene_level_copy_number", "*.gene_level_copy_number.tsv", full.names = T, recursive = T)
names(fls.tcga.copy) <- sapply(strsplit(fls.tcga.copy,"/"), "[[", 3)

tcga.copy <- rbindlist(lapply(fls.tcga.copy, fread), idcol = T)
tcga.copy[tcga.anno$sample, sample_id := sample_id, on = ".id == `File ID`"]

# expression
fls.tcga.exprs <- list.files("TCGA/Gene_expression_quantification", "*.FPKM.txt.gz", full.names = T, recursive = T)
names(fls.tcga.exprs) <- sapply(strsplit(fls.tcga.exprs,"/"), "[[", 3)

tcga.exprs <- rbindlist(lapply(fls.tcga.exprs, fread, col.names = c("gene_id", "fpkm")), idcol = T)
tcga.exprs[, tpm := fpkm/sum(fpkm,na.rm=T)*1e6, by=".id"]
tcga.exprs[, tpm.avg := mean(tpm,na.rm=T), by="gene_id"]
tcga.exprs[tcga.anno$sample, c("sample_id", "sample_type") := list(sample_id, `Sample Type`), on = ".id == `File ID`"]
tcga.exprs[tcga.anno$gene, c("chrom", "start") := list(seqname, start), on = c("gene_id == gene_id")]
tcga.exprs[tcga.copy, "copies" := list(copy_number), on = c("sample_id == sample_id", "gene_id == gene_id")]

tcga.exprs.avg <- tcga.exprs[tpm.avg > 1 & !is.na(copies), base::mean(tpm, trim=0.2, na.rm=T), by=c("copies", "sample_id", "sample_type")]
tcga.exprs.avg[sample_type == "Primary Tumor", exp := mean(V1[copies == 2], na.rm=T)*(copies/2)]

# merge segments with a smaller absolute diff than 0.2
tcga.copy.seg <- tcga.masked[Num_Probes > 1e3, list(seg_start = min(Start), seg_end = max(End)), by=c("sample_id", "Chromosome", "Segment_Mean")]

tcga.copy.seg[, seg.grp := cumsum(c(TRUE, abs(diff(Segment_Mean )) > 0.2)), by=c("sample_id", "Chromosome")]
tcga.copy.seg[, N := .N, by=c("sample_id", "Chromosome", "seg.grp")]

tcga.copy.seg.red <- tcga.copy.seg[N>1, list(start=min(seg_start), end=max(seg_end), mean=weighted.mean(Segment_Mean, seg_end-seg_start)), by=c("sample_id", "Chromosome", "seg.grp")]
tcga.copy.seg.red[, chrom := paste0("chr", Chromosome)]
tcga.copy.seg.red[, seg := seq_along(seg.grp), by=c("sample_id", "Chromosome")]

# add segments
tcga.exprs[, temp := start] # dummy variable for overlap
setkey(tcga.exprs, sample_id, chrom, start, temp)
setkey(tcga.copy.seg.red, sample_id, chrom, start, end)

tcga.exprs[, c("seg") := tcga.copy.seg.red[foverlaps(tcga.exprs, tcga.copy.seg.red, mult="first", which = T), list(seg)]]
tcga.exprs[, temp := NULL] # cleanup

## PLOT
library(ggplot2)
library(cowplot)
# expression v copies
p.tcga.exprs.box <-
ggplot(tcga.exprs.avg[copies %between% c(1,4)], aes(x=factor(copies), y=V1 )) + 
  geom_boxplot(outlier.stroke = NA, outlier.color = "grey50", show.legend = F, aes(fill=factor(copies-2) )) +
  geom_point(data=tcga.exprs.avg[copies != 2 & copies %between% c(1,4), median(exp,na.rm=T), by=c("copies", "sample_type")]) +
  facet_grid(~sample_type) +
  labs(x="Copies", y="RNA expression") +
  coord_cartesian(ylim=c(0,30)) +
  scale_fill_manual(values=c("#93C5DE", "#F8F7F7", "#F5A582", "#CA1523")) +
  theme_cowplot() +
  theme(strip.background = element_blank())

ggsave2("plots/tcga_expression_boxplot.pdf", height = 3, width = 3, p.tcga.exprs.box)

# linear slope
p.tcga.exprs.line <-
ggplot(tcga.exprs.avg[sample_type == "Primary Tumor" & copies > 0 & n > 2], aes(x=log2(copies/2), y=log2(V1+1) )) +
  geom_point(col="grey50", stroke=NA, alpha=0.33) +
  geom_smooth(method="lm", col="black") +
  geom_smooth(method="lm", se=F, lty=2, col="black", aes(y=log2(exp+1))) +
  annotate("text", x=-1, y=7, hjust=0, label=gsub(" ", "\n", tcga.exprs.avg[sample_type == "Primary Tumor" & copies > 0 & n > 2, compare_slopes(lm(log2(V1+1)~log2(copies/2)), lm(log2(exp+1)~log2(copies/2)))]) ) +
  labs(x="Relative copies (log2)", y="Relative expression (log2+1)") +
  coord_cartesian(ylim=c(0,8)) +
  theme_cowplot()

ggsave2("plots/tcga_expression_line.pdf", height = 3, width = 3, p.tcga.exprs.line)

# chromosome plot for one sample
# TCGA-AA-3973-01A
# TCGA-AA-A00L-01A
idx <- "TCGA-AA-3693-01A"

p.tcga.expr.01A <-
ggplot(tcga.exprs[sample_id == idx & tpm.avg > 1 & !is.na(seg) & chrom != "chrY"], aes(x=factor(seg), y=tpm)) +
  stat_summary(fun.data="median_cl_boot", stroke=NA) +
  geom_hline(data=tcga.exprs[sample_id == idx & tpm.avg > 1 & !is.na(seg) & chrom != "chrY" & copies==2, list(y=c(0.5,1,1.5,2)*median(tpm), col=LETTERS[1:4])],lty=2, show.legend = F, aes(yintercept=y, col=col)) +
  facet_grid(~factor(chrom,levels=paste0("chr",c(1:22,"X"))), scales="free_x") +
  labs(x=NULL, y="Expression", title=idx) +
  coord_cartesian(ylim=c(0,15)) +
  scale_color_manual(values=c("#A3CBE4", "#999999", "#F4AB86", "#C81618")) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), panel.spacing.x = unit(0,"mm"))

p.tcga.seg.01A <-
ggplot(tcga.masked[sample_type == "Primary Tumor" & Num_Probes > 1e3 & sample_id == idx], aes(x=Start, y=Segment_Mean)) + 
  geom_point(stroke=NA) +
  geom_segment(data=tcga.copy.seg.red[sample_id == idx], col="#FF7F00", aes(x=start, xend=end, y=mean)) +
  facet_grid(~factor(Chromosome,levels=c(1:22,"X")), scales="free_x") +
  coord_cartesian(ylim=c(-1,1)) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing.x = unit(0,"mm"))

library(patchwork)
ggsave2("plots/tcga_expression_singlesample.pdf",height = 4, width = 10, p.tcga.expr.01A/p.tcga.seg.01A)

# all individual samples
samps <- c("TCGA-AA-A00L-01A", "TCGA-CM-5864-01A", "TCGA-AA-3679-01A", "TCGA-G4-6322-01A", "TCGA-AA-3975-01A")

ls.p.tcga <- lapply(samps, function(i){
  p.expr <-
  ggplot(tcga.exprs[sample_id == i & tpm.avg > 1 & !is.na(seg) & chrom != "chrY"], aes(x=factor(seg), y=tpm)) +
    stat_summary(fun.data="median_cl_boot", stroke=NA) +
    geom_hline(data=tcga.exprs[sample_id == i & tpm.avg > 1 & !is.na(seg) & chrom != "chrY" & copies==2, list(y=c(0.5,1,1.5,2)*median(tpm), col=LETTERS[1:4])],lty=2, show.legend = F, aes(yintercept=y, col=col)) +
    facet_grid(~factor(chrom,levels=paste0("chr",c(1:23,"X"))), scales="free_x") +
    labs(x=NULL, y="Expression", title=i) +
    coord_cartesian(ylim=c(0,15)) +
    scale_color_manual(values=c("#A3CBE4", "#999999", "#F4AB86", "#C81618")) +
    theme_cowplot() +
    theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), panel.spacing.x = unit(0,"mm"))
  
  p.seg <-
  ggplot(tcga.masked[sample_type == "Primary Tumor" & Num_Probes > 1e3 & sample_id == i], aes(x=Start, y=Segment_Mean)) + 
    geom_point(stroke=NA) +
    geom_segment(data=tcga.copy.seg.red[sample_id == i], col="#FF7F00", aes(x=start, xend=end, y=mean)) +
    facet_grid(~factor(Chromosome,levels=c(1:22,"X")), scales="free_x") +
    coord_cartesian(ylim=c(-1,1)) +
    theme_cowplot() +
    theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing.x = unit(0,"mm"))

  library(patchwork)
  p.expr/p.seg
})

ggsave2("plots/tcga_expression_singlesample_all.pdf",height = length(samps)*4, width = 10, limitsize = F, plot_grid(plotlist=ls.p.tcga, ncol = 1))

##