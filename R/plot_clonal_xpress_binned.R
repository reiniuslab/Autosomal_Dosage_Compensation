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

##

dat.normumi.xpress <- readRDS("data/dat_normumi_xpress.rds")

# load protein class annotations
orth <- fread("humanproteinatlas/mouse_human_marmoset_macaque_orthologs_20231113.csv")
anno.class <- rbindlist(list(
  "enzymes" = fread("humanproteinatlas/protein_class_Enzymes_Any.tsv")[,c(1,3,5)],
  "transporters" = fread("humanproteinatlas/protein_class_Transporters_Any.tsv")[,c(1,3,5)],
  "metabolic" = fread("humanproteinatlas/protein_class_Metabolic.tsv")[,c(1,3,5)],
  "transcription" = fread("humanproteinatlas/protein_class_Transcription.tsv")[,c(1,3,5)]
), idcol = "class")

anno.class[orth, EnsemblMouse := mouse_EnsemblID, on = c("Ensembl == human_EnsemblID")]



##
dat.normumi.xpress[normumi > 0 & !is.na(segment) & chrom %in% 1:19 & note != "Exclude", exprs_bin := factor(as.numeric(cut(normumi, quantile(normumi, seq(0,1, length.out=5))))), by="ploidy"]

# summarise per clone and expression bin
dat.normumi.avggnbin.xpress <- dat.normumi.xpress[!is.na(segment) & chrom %in% 1:19 & note != "Exclude", list("normumi" = mean(normumi)), by=c("GeneID", "clone", "ploidy", "copies", "exprs_bin")]
dat.normumi.avggnbin.xpress[, relative := normumi / median(normumi, na.rm=T), by=c("ploidy", "clone", "exprs_bin")]

# summarise per gene and expand and summarise to protein class
dat.normumi.avggn.xpress <- dat.normumi.xpress[normumi > 0 & !is.na(segment) & chrom %in% 1:19 & note != "Exclude", list("normumi" = mean(normumi)), by=c("GeneID", "clone", "ploidy", "copies")]
dat.normumi.avggn.class.xpress <- dat.normumi.avggn.xpress[anno.class, on = c("GeneID == EnsemblMouse")][,c("GeneID", "ploidy", "copies", "clone", "normumi", "class")]
dat.normumi.avggn.class.xpress[, relative := normumi / median(normumi, na.rm=T), by=c("ploidy", "clone", "class")]


#dat.normumi.avggn.class.xpress[, obsexp.class := normumi / (mean(normumi[copies == ploidy], na.rm=T)*(copies/ploidy)), by=c("ploidy", "class")]

## PLOT
library(ggplot2)
library(cowplot)

# expression by copies binned by expression
p.exprs.bin <-
ggplot(dat.normumi.avggnbin.xpress[!is.na(exprs_bin)], aes(x=exprs_bin, y=log2(relative), fill=factor(copies-ploidy))) + 
  geom_boxplot(outlier.stroke = NA, outlier.colour = "#999999") + 
  geom_hline(data=dat.normumi.avggnbin.xpress[!is.na(exprs_bin), median(relative), by = c("copies", "ploidy", "exprs_bin")][,log2(median(V1)), by=c("copies","ploidy")], aes(yintercept=V1, col=factor(copies-ploidy))) +
  facet_grid(~ploidy, scales="free_x", space="free_x") + 
  labs(x="Expression bin", y="Expression (log2 relative)") +
  coord_cartesian(ylim=c(-0.75,0.75)) +
  scale_fill_manual(values=c("#0573B2", "#93C5DE", "#F8F7F7", "#F5A582"), name="DNA copies\n(relative)", aesthetics = c("col","fill")) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour = "black"))

ggsave2("plots/xpress_exprs_binned.pdf", width = 6, height = 3, p.exprs.bin)

# expression by copies binned by protein class
p.exprs.class <-
ggplot(dat.normumi.avggn.class.xpress[!is.na(GeneID) & !is.na(class) & !is.na(relative)], aes(x=class, y=log2(relative), fill=factor(copies-ploidy) )) +
  geom_boxplot(outlier.stroke = NA, outlier.colour = "#999999") +
  facet_grid(~ploidy, scales="free_x", space="free_x") +
  labs(x="Protein class", y="Expression (log2 relative)") +
  coord_cartesian(ylim=c(-4,4)) +
  scale_fill_manual(values=c("#0573B2", "#93C5DE", "#F8F7F7", "#F5A582"), name="DNA copies\n(relative)") +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour = "black"))

ggsave2("plots/xpress_exprs_class.pdf", width = 6, height = 3, p.exprs.class)

##
