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

## Load data
dat.burst.xpress <- readRDS("data/dat_burst_xpress.rds")


## PLOT
library(ggplot2)
library(cowplot)
library(patchwork)

# burst freq v size 2d comparison
dat.burst.xpresss.cast <- dcast(dat.burst.xpress[unit=="normumi" & !is.na(ploidy) & note != "Exclude", mean(value,na.rm=T), by=c("gene","metric","allele","ploidy", "tpm.avg")], gene+allele+ploidy+tpm.avg~metric)

p.kinetics.comp <-
ggplot(dat.burst.xpresss.cast[, lapply(.SD, median, na.rm=T), by=c("gene", "allele","ploidy"), .SDcols=c("tpm.avg", "bf", "bs")], aes(x=bs, y=bf, z=tpm.avg )) +
  ggrastr::rasterise(stat_summary_2d(fun = function(x) log2(median(x,na.rm=T)+1), bins = 100 )) +
  labs(x="Burst size", y="Burst frequency", fill="Expression\n(TPM)") +
  coord_cartesian(ylim=c(0.003,6), xlim=c(1,300)) +
  #facet_wrap(~ploidy) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  scale_fill_viridis_c(limits = c(0,15), oob=scales::squish) +
  theme_cowplot() +
  theme(aspect.ratio = 1)

ggsave2("plots/kinetics_xpress_comp.pdf", width = 6, height = 4, p.kinetics.comp)

# burst freq v size density comparison
# earths mover distance between density distributions
library(emdist)
emdist::emd2d(
  MASS::kde2d(dat.burst.xpresss.cast[ploidy=="2" & (!is.na(bf) & !is.na(bs)),bs], dat.burst.xpresss.cast[ploidy=="2" & (!is.na(bf) & !is.na(bs)),bf])$z,
  MASS::kde2d(dat.burst.xpresss.cast[ploidy=="4" & (!is.na(bf) & !is.na(bs)),bs], dat.burst.xpresss.cast[ploidy=="4" & (!is.na(bf) & !is.na(bs)),bf])$z
)

p.kinetics.compdens <-
ggplot(dat.burst.xpresss.cast, aes(x=bs, y=bf, col=factor(ploidy) )) +
  geom_density2d() +
  labs(x="Burst size", y="Burst frequency", col="Ploidy") +
  coord_cartesian(ylim=c(0.03,6), xlim=c(1,10)) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  scale_colour_brewer(palette="Set1", direction=1) +
  theme_cowplot() +
  theme(aspect.ratio = 1)

ggsave2("plots/kinetics_xpress_compdens.pdf", width = 6, height = 4, p.kinetics.compdens)

# burst kinetics per clone, segment and chromosome
p.kinetics.bf.tpm <-
ggplot(dat.burst.xpress[metric == "bf" & unit == "tpm" & !is.na(segment) & chrom %in% c(1:19,"X")],aes(x=factor(segment), y=relative.burst, col=allele)) +
  stat_summary(fun.data="median_cl_boot") +
  facet_wrap(~ploidy+clone+chrom, scales="free_x", ncol=20) +
  labs(x=NULL, y="Burst frequency (relative)") +
  coord_cartesian(ylim=c(0,3.5)) +
  scale_color_brewer(palette="Set1") +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour = "black"))

ggsave2("plots/kinetics_xpress_bf_tpm.pdf", width = 18, height = 22, p.kinetics.bf.tpm)

p.kinetics.bf.normumi <-
ggplot(dat.burst.xpress[metric == "bf" & unit == "normumi" & !is.na(segment) & chrom %in% c(1:19,"X")],aes(x=factor(segment), y=relative.burst, col=allele)) +
  stat_summary(fun.data="median_cl_boot") +
  facet_wrap(~ploidy+clone+chrom, scales="free_x", ncol=20) +
  labs(x=NULL, y="Burst frequency (relative)") +
  coord_cartesian(ylim=c(0,2)) +
  scale_color_brewer(palette="Set1") +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour = "black"))

ggsave2("plots/kinetics_xpress_bf_normumi.pdf", width = 18, height = 22, p.kinetics.bf.normumi)

p.kinetics.bs.normumi <-
ggplot(dat.burst.xpress[metric == "bs" & unit == "normumi" & !is.na(segment) & chrom %in% c(1:19,"X")],aes(x=factor(segment), y=relative.burst, col=allele)) +
  stat_summary(fun.data="median_cl_boot") +
  facet_wrap(~ploidy+clone+chrom, scales="free_x", ncol=20) +
  labs(x=NULL, y="Burst size (relative)") +
  coord_cartesian(ylim=c(0,2)) +
  scale_color_brewer(palette="Set1") +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour = "black"))

ggsave2("plots/kinetics_xpress_bs_normumi.pdf", width = 18, height = 22, p.kinetics.bs.normumi)

# Burst frequency per copy, only for the major allele (nalleles >= ploidy*0.5)
p.kinetics.bf.copy <-
ggplot(dat.burst.xpress[nalleles >= ploidy*0.5 & metric == "bf" & unit == "tpm" & !is.na(segment) & chrom %in% 1:19 & note != "Exclude"], aes(x=factor(copies), col=factor(copies-ploidy), y=relative.burst)) +
  stat_summary(fun.data="median_cl_boot", show.legend = F) +
  facet_grid(~ploidy, scales="free_x", space="free_x") +
  labs(x="DNA copies", y="Burst frequency (relative)") +
  coord_cartesian(ylim=c(0,3)) +
  scale_color_brewer(palette="RdBu", direction = -1) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"))

p.kinetics.bs.copy <-
ggplot(dat.burst.xpress[nalleles >= ploidy*0.5 & metric == "bs" & unit == "normumi" & !is.na(segment) & chrom %in% 1:19 & note != "Exclude"], aes(x=factor(copies), col=factor(copies-ploidy), y=relative.burst)) +
  stat_summary(fun.data="median_cl_boot", show.legend = F) +
  facet_grid(~ploidy, scales="free_x", space="free_x") +
  labs(x="DNA copies", y="Burst size (relative)") +
  coord_cartesian(ylim=c(0,3)) +
  scale_color_brewer(palette="RdBu", direction = -1) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"))

ggsave2("plots/kinetics_xpress_copy.pdf", width = 3, height = 4, p.kinetics.bf.copy / p.kinetics.bs.copy)

## density plot
# center data to neutral
dat.burst.xpress[!is.na(relative.burst) & !is.na(relative.expression) & chrom %in% 1:19 & note != "Exclude" & normumi > 1, c("relative.burst2", "relative.expression2") := list(value/median(value[nalleles == ploidy/2]), allelic.expression/median(allelic.expression[nalleles == ploidy/2])), by=c("L1","variable")]

p.kinetics.bf.density2d <-
ggplot(dat.burst.xpress[nalleles >= ploidy*0.5 & metric == "bf" & unit == "tpm" & !is.na(segment) & chrom %in% 1:19 & note != "Exclude" & normumi > 1], aes(x=log2(relative.expression2), y=log2(relative.burst2), col=factor(copies-ploidy) )) + 
  geom_density2d() +
  geom_hline(yintercept = 0, col= "#999999", lty=2, show.legend = F) +
  geom_vline(xintercept = 0, col= "#999999", lty=2, show.legend = F) +
  facet_wrap(~ploidy) +
  labs(x="Relative expression (log2)", y="Burst frequency (log2 relative)") +
  scale_color_manual(values = c("#0573B2", "#93C5DE", "#999999", "#F5A582"), name="DNA copies") +
  coord_cartesian(xlim=c(-3.5,3.5), ylim=c(-3.5,3.5)) +
  theme_cowplot() +
  theme(aspect.ratio = 1, strip.background = element_blank(), panel.border = element_rect(colour="black"))

p.kinetics.bs.density2d <-
ggplot(dat.burst.xpress[nalleles >= ploidy*0.5 & metric == "bs" & unit == "normumi" & !is.na(segment) & chrom %in% 1:19 & note != "Exclude" & normumi > 1], aes(x=log2(relative.expression2), y=log2(relative.burst2), col=factor(copies-ploidy) )) + 
  geom_density2d() +
  geom_hline(yintercept = 0, col= "#999999", lty=2, show.legend = F) +
  geom_vline(xintercept = 0, col= "#999999", lty=2, show.legend = F) +
  facet_wrap(~ploidy) +
  labs(x="Relative expression (log2)", y="Burst size (log2 relative)") +
  scale_color_manual(values = c("#0573B2", "#93C5DE", "#999999", "#F5A582"), name="DNA copies") +
  coord_cartesian(xlim=c(-3.5,3.5), ylim=c(-3.5,3.5)) +
  theme_cowplot() +
  theme(aspect.ratio = 1, strip.background = element_blank(), panel.border = element_rect(colour="black"))

ggsave2("plots/kinetics_xpress_exprs_density.pdf", width = 6, height = 6, p.kinetics.bf.density2d/p.kinetics.bs.density2d)

##
