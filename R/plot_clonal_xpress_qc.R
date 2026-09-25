library(data.table)

## Load data
library(scuttle)
dat.xpress.qc <- rbind(
  as.data.table(colData(readRDS("data/sce_xpress1.rds"))),
  as.data.table(colData(readRDS("data/sce_xpress2.rds")))
)


dat.normumi.xpress <- readRDS("data/dat_normumi_xpress.rds")

## PLOT
library(ggplot2)
library(cowplot)

p.qc.sat.xpress <-
ggplot(dat.xpress.qc[order(spike.pct, decreasing = F)], aes(x=sum, y=detected, col=spike.pct*100)) +
  ggrastr::geom_point_rast(scale=0.75) +
  geom_hline(yintercept= median(sce.xpress$detected) - 3*mad(sce.xpress$detected), lty=2, col="grey50") +
  annotate("text", x=4e5, y=9e3, label=paste0(table(sce.xpress$outlier),collapse="\n"), col="grey50") +
  labs(x="Reads", y="Detected genes", col="%Spike-in") +
  coord_cartesian(ylim=c(0,15e3)) +
  scale_color_viridis_c(limits=c(0,10)) +
  theme_cowplot() +
  theme(aspect.ratio = 1)

p.qc.det.xpress <-
ggplot(dat.normumi.xpress[normumi > 0, .N, by="Var2"], aes(x=N)) +
  geom_histogram() +
  labs(x="Detected genes", y="Count") +
  #expand_limits(x=0) +
  theme_cowplot()

p.qc.allele.xpress <-
ggplot(dat.normumi.xpress[normumi > 0, mean((c57+cast)>0,na.rm=T), by="Var2"], aes(x=V1)) +
  geom_histogram() +
  labs(x="Fraction allelic detection", y="Count") +
  #expand_limits(x=0) +
  theme_cowplot()

## EXPORT
library(patchwork)
ggsave2("plots/xpress_qc.pdf", width = 5, height = 6, p.qc.sat.xpress / p.qc.det.xpress / p.qc.allele.xpress)  
