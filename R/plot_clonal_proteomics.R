library(data.table)

##

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

#
dat.normumi.xpress <- readRDS("data/dat_normumi_xpress.rds")

dat.normumi.avgploidy.xpress <- dat.normumi.xpress[normumi > 0 & !is.na(segment) & chrom %in% 1:19 & note != "Exclude", list(base::mean(normumi, trim = 0.2, na.rm=T), .N), by=c("Var2", "clone", "ploidy", "copies")]
dat.normumi.avgploidy.xpress[, relative := V1 / median(V1[copies == ploidy], na.rm=T), by=c("clone")]

#

dat.prot <- rbindlist(list(
  melt(fread("proteomics/20250122_TMT18_Exp1_Frac_Ratio.txt"), id.vars=c("accession", "gene_id", "gene_name")),
  melt(fread("proteomics/20250122_TMT18_Exp2_Frac_Ratio.txt"), id.vars=c("accession", "gene_id", "gene_name")),
  melt(fread("proteomics/20240127_TMT18_Exp3_Frac_Ratio.txt"), id.vars=c("accession", "gene_id", "gene_name"))
), fill=T)
dat.prot[dat.normumi.xpress, c("ploidy", "copies", "segment", "chrom", "note") := list(ploidy, copies, segment, chrom, note), on = c("gene_id == GeneID", "variable == clone")]

dat.prot.avg <- dat.prot[chrom %in% 1:19 & !is.na(segment) & note != "Exclude", median(value, na.rm=T), by=c("variable", "copies", "ploidy")]
dat.prot.avg[, exp := mean(V1[copies == ploidy], na.rm=T)*(copies/ploidy), by=c("ploidy")]

## PLOT
library(ggplot2)
library(cowplot)

# linear regression
dat.prot.avg[, summary(lm(V1~log2(copies)+variable))]

compare_slopes(
  dat.prot.avg[, lm(V1~log2(copies)+variable)],
  dat.prot.avg[, lm(exp~log2(copies)+variable)]
) # v expected
compare_slopes(
  dat.prot.avg[, lm(V1~log2(copies)+variable)],
  dat.normumi.avgploidy.xpress[, lm(relative~log2(copies)+clone)]
) # v RNA

# rna vs protein compensation levels
dat.prot.comp <- cbind(dat.prot.avg[order(ploidy, copies), list("prot" = median(V1)), by=c("copies", "ploidy")], "rna" = dat.normumi.avgploidy.xpress[, median(relative, na.rm=T), by=c("ploidy", "copies")][order(ploidy,copies), V1])
dat.prot.comp[, abs(1-prot/rna), by=c("copies","ploidy")][copies != ploidy, range(V1)]

## protein levels vs copies
p.proteomics.copy <-
ggplot(dat.prot.avg, aes(x=factor(copies), y=V1 )) +
  geom_boxplot(outlier.shape = NA, show.legend = F, aes(fill=factor(copies-ploidy) )) +
  ggbeeswarm::geom_quasirandom(colour="grey50") +
  geom_point(data=dat.normumi.avgploidy.xpress[, median(relative, na.rm=T), by=c("ploidy", "copies")], show.legend = F, aes(y=V1, col=factor(copies-ploidy) )) +
  geom_point(data=unique(dat.prot.avg[copies != ploidy, list(ploidy, copies, exp)]), show.legend = F, aes(y=exp)) +
  geom_text(data=dat.prot.avg[copies != ploidy, t.test(V1, mu=exp[1])$p.value, by=c("ploidy", "copies")], nudge_x = -0.25, angle = 90, hjust = 1, aes(y=1.5, label=paste0("P = ",signif(V1,2) ) )) +
  labs(x="DNA copies", y="Relative protein abundance") +
  facet_grid(~ploidy, scales="free_x", space="free_x") +
  coord_cartesian(ylim=c(0.5,1.5)) +
  scale_fill_manual(values=c("#0573B2", "#93C5DE", "#F8F7F7", "#F5A582"), aesthetics = c("fill", "col")) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"))

ggsave2("plots/proteomics_ratio_copy.pdf", height = 4, width = 3, p.proteomics.copy)

## correlation vs copies
p.proteomics.corr <-
ggplot(dat.prot[chrom %in% 1:19 & note != "Exclude", median(value, na.rm=T), by=c("variable", "copies", "ploidy")], aes(x=log2(copies/ploidy), y=log2(V1/(copies/ploidy)) )) +
  geom_point(aes(col=paste(ploidy, variable))) +
  geom_line(aes(group=variable, col=paste(ploidy, variable))) +
  geom_smooth(method="lm", lty=2, fullrange=T) +
  geom_smooth(method="lm", se=F) +
  labs(x="DNA copies (log2 relative)", y="Relative protein abundance (log2 obs/exp)") +
  facet_grid(~ploidy, scales="free_x", space="free_x") +
  coord_cartesian(ylim=c(-1,1)) +
  scale_x_continuous(limits=c(-1,1)) +
  scale_color_brewer(palette="Set2", name=NULL) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"))

ggsave2("plots/proteomics_obsexp_correlation.pdf", height = 3, width = 5, p.proteomics.corr)

## correlation protein vs rna
dat.prot2 <- rbindlist(list(
  melt(fread("proteomics/20250122_TMT18_Exp1_Grouped_Abundance.txt"), id.vars=c("accession", "gene_id", "gene_name")),
  melt(fread("proteomics/20250122_TMT18_Exp2_Grouped_Abundance.txt"), id.vars=c("accession", "gene_id", "gene_name")),
  melt(fread("proteomics/20240127_TMT18_Exp3_Grouped_Abundance.txt"), id.vars=c("accession", "gene_id", "gene_name"))
), fill=T)
dat.prot2[dat.normumi.xpress, c("ploidy", "copies", "segment", "chrom", "note") := list(ploidy, copies, segment, chrom, note), on = c("gene_id == GeneID", "variable == clone")]

# add pseudobulk norm-umi counts
dat.prot2[dat.normumi.xpress[, mean(normumi), by=c("GeneID", "clone")], normumi := i.V1, on = c("gene_id == GeneID", "variable == clone")]

p.protrna.cor.2n <- 
ggplot(dat.prot2[value > 1 & normumi >1 & chrom %in% 1:19 & note != "Exclude" & ploidy == 2], aes(x=log10(normumi), y=log10(value) )) +
  ggrastr::rasterise(ggpointdensity::geom_pointdensity(stroke=NA, size=1, method.args = list(n = 100), aes(color = after_stat(ndensity)) ), dpi=600) +
  geom_smooth(method="lm", col="black", se = F) +
  geom_text(data=dat.prot2[value > 1 & normumi >1 & chrom %in% 1:19 & note != "Exclude" & ploidy == 2, round(cor.test(log10(value),log10(normumi), method="spearman")$estimate, 2), by=c("variable")], hjust = 0, vjust=1, aes(x=0, y=6, label=V1)) +
  facet_wrap(~variable, nrow=1) +
  labs(x="UMI (log10)", y="Protein abundance (log10)") +
  coord_cartesian(xlim=c(0,3), ylim=c(0.5,6)) +
  scale_color_viridis_c(name="Density", limits=c(0,1)) +
  theme_cowplot() +
  theme(aspect.ratio=1, strip.background = element_blank(), panel.border = element_rect(colour="black"))

p.protrna.cor.4n <- 
ggplot(dat.prot2[value > 1 & normumi >1 & chrom %in% 1:19 & note != "Exclude" & ploidy == 4], aes(x=log10(normumi), y=log10(value) )) +
  ggrastr::rasterise(ggpointdensity::geom_pointdensity(stroke=NA, size=1, method.args = list(n = 100), aes(color = after_stat(ndensity)) ), dpi=600) +
  geom_smooth(method="lm", col="black", se = F) +
  geom_text(data=dat.prot2[value > 1 & normumi >1 & chrom %in% 1:19 & note != "Exclude" & ploidy == 4, round(cor.test(log10(value),log10(normumi), method="spearman")$estimate, 2), by=c("variable")], hjust = 0, vjust=1, aes(x=0, y=6, label=V1)) +
  facet_wrap(~variable, nrow=1) +
  labs(x="UMI (log10)", y="Protein abundance (log10)") +
  coord_cartesian(xlim=c(0,3), ylim=c(0.5,6)) +
  scale_color_viridis_c(option = "C", name="Density", limits=c(0,1)) +
  theme_cowplot() +
  theme(aspect.ratio=1, strip.background = element_blank(), panel.border = element_rect(colour="black"))

library(patchwork)
ggsave2("plots/proteomics_protein_rna_correlation.pdf", height = 3, width = 10, p.protrna.cor.2n+p.protrna.cor.4n+plot_layout(widths=c(3,5.1)))

##
