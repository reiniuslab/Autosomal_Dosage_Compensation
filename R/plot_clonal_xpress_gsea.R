library(data.table)

##

dat.normumi.xpress <- readRDS("data/dat_normumi_xpress.rds")

v.fc.xpress <- dat.normumi.xpress[normumi > 0, log2(mean(normumi[ploidy == 4],na.rm=T)/mean(normumi[ploidy == 2],na.rm=T)), by = c("GeneID")][!is.na(V1), sort(setNames(V1,GeneID), decreasing=T)]

# GSEA
library(clusterProfiler)
library(org.Mm.eg.db)
dat.gsea.ploidy <- as.data.table(simplify(gseGO(v.fc.xpress, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL"))@result)

fwrite(dat.gsea.ploidy, "data/dat_gsea_ploidy.tsv", quote = F, sep = "\t")

## PLOT
library(ggplot2)
library(cowplot)

p.gsea.ploidy <-
ggplot(dat.gsea.ploidy[qvalue < 1e-4], aes(x=-log10(qvalue), y=reorder(Description,-log10(qvalue), mean), col= NES, size= -log10(qvalue) )) +
  geom_point() +
  labs(x="-log10 FDR (4N vs. 2N)", y=NULL) +
  scale_color_distiller(palette="Reds", direction = 1, limits=c(1,2.25)) +
  theme_cowplot()

ggsave2("plots/xpress_gsea_ploidy.pdf", width = 9, height = 6)
