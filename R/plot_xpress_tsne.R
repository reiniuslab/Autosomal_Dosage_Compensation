library(data.table)

##

## load data
# cnv data
dat.clonal.cnv <- fread("data/dat_clonal_cnv_fixed.tsv")
dat.clonal.cnv[, variable := gsub("_repeat","",variable)]

# scrna-seq data
library(scater)
library(scran)
sce.xpress <- readRDS("data/sce_xpress_scvi.rds")

# add information
sce.xpress$ploidy <- factor(dat.clonal.cnv[!duplicated(variable)][match(sce.xpress$clone, variable), round(ploidy)])
sce.xpress$aneuploidy_score <- dat.clonal.cnv[!duplicated(variable)][match(sce.xpress$clone, variable), aneuploidy_score]
sce.xpress$treatment <- substr(sce.xpress$clone, 1,1)

## PLOT
library(ggplot2)
library(cowplot)
# umi counts
dat.umicnt <- data.table(barcode=sce.xpress$barcode ,clone=sce.xpress$clone,counts=colSums(assay(sce.xpress,"counts")), normcounts = colSums(assay(sce.xpress,"normcounts")), ploidy=sce.xpress$ploidy, aneuploidy_score=sce.xpress$aneuploidy_score, dataset=sce.xpress$dataset)

p.xpress.umicounts <- 
ggplot(dat.umicnt, aes(x=ploidy, y=normcounts/1e3)) +
  geom_violin() +
  geom_boxplot(width = 0.25, outlier.shape = NA) +
  labs(y="Spike-normalised UMI counts (x1,000)", x="Ploidy") +
  coord_cartesian(ylim=c(0,600)) +
  theme_cowplot()

# t-sne
p.tsne.batch <- plotTSNE(sce.xpress, colour_by="dataset", order_by = "barcode") + scale_color_brewer(palette="Set2", name="dataset") + theme(aspect.ratio = 1)
p.tsne.clone <- plotTSNE(sce.xpress, colour_by="clone", order_by = "barcode") + theme(aspect.ratio = 1)
p.tsne.ploidy <- plotTSNE(sce.xpress, colour_by="ploidy", order_by = "barcode") + scale_color_brewer(palette="Set1", name="ploidy", direction = -1) + theme(aspect.ratio = 1)
p.tsne.aneusc <- plotTSNE(sce.xpress, colour_by="aneuploidy_score", order_by = "barcode") + theme(aspect.ratio = 1)
p.tsne.treat <- plotTSNE(sce.xpress, colour_by="treatment", order_by = "barcode") + theme(aspect.ratio = 1)

## EXPORT
library(patchwork)
ggsave2("plots/xpress_umicounts.pdf", height = 4, width = 3, p.xpress.umicounts)
ggsave2("plots/xpress_tsne.pdf", width = 9, height = 5, p.tsne.batch+p.tsne.clone+p.tsne.ploidy+p.tsne.aneusc+p.tsne.treat)

##