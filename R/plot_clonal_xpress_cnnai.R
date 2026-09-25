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

# compare CNN-AI genes
dat.normumi.xpress[normumi > 0 & !is.na(segment) & note != "Exclude", rel := normumi/median(normumi[chrom %in% 1:19 & copies == ploidy]), by=c("Var2")]
dat.normumi.avgcnn.xpress <- dat.normumi.xpress[normumi > 0 & !is.na(segment) & note != "Exclude" & copies == ploidy & GeneID %in% dat.normumi.xpress[note == "CNN-AI", unique(GeneID)], mean(rel, na.rm=T), by=c("GeneID", "chrom", "note")]
dat.normumi.avgcnn.xpress[, n := .N, by="GeneID"]

## PLOT
library(ggplot2)
library(cowplot)

p.rel.cnn <-
ggplot(dat.normumi.avgcnn.xpress[n == 2], aes(x=note, y=log2(V1) )) + 
  ggrastr::rasterise(geom_line(alpha=0.1, col="grey50", aes(group=GeneID)), dpi=100, scale=0.5) +
  geom_point(stroke = NA, col="grey50") +
  geom_boxplot(outlier.shape = NA) +
  annotate("text", x=1.5, y=10, label=paste0("Paired Wilcoxon\nP = ", signif(dat.normumi.avgcnn.xpress[n == 2][order(GeneID), wilcox.test(V1[note == "CNN-AI"], V1[note == ""], paired=T)$p.value], 4))) +
  labs(x=NULL, y="Relative expression (log2)") +
  coord_cartesian(ylim=c(-5,10)) +
  theme_cowplot()

ggsave2("plots/xpress_exprs_cnn.pdf", width = 3, height = 3, p.rel.cnn)

##
