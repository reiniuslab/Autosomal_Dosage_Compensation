library(data.table)

dat.normumi.xpress <- readRDS("data/dat_normumi_xpress.rds")
dat.normumi.xpress[, pos_100kb := list(round(pos/1e5)*1e5)]

dat.umifrac.sum <- dat.normumi.xpress[tpm.avg > 1 & gene.det > 0.25 & !is.na(c57), sum(c57, na.rm=T)/sum(c57+cast,na.rm=T), by=c("clone", "Var2", "chrom", "pos_100kb")]

library(ggplot2)
library(cowplot)
p.allelefrac.heat <-
ggplot(dat.umifrac.sum[!is.na(V1) & chrom %in% c(1:19,"X")], aes(x=factor(pos_100kb), y=Var2, fill=V1)) +
  ggrastr::geom_tile_rast() + 
  facet_grid(clone~chrom, scales="free", space="free") +
  scale_fill_distiller(palette="Spectral", name="C57 fraction") +
  labs(x="Gene position", y="Cells") +
  theme_cowplot() +
  theme(strip.background = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title.y = element_blank(), panel.border = element_rect(colour="black"), panel.spacing.x = unit(0,"mm"))

ggsave2("plots/xpress_allelic_heatmap.pdf", width = 12, height = 12, p.allelefrac.heat)
