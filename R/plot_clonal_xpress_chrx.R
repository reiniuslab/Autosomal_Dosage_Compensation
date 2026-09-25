library(data.table)

dat.normumi.xpress <- readRDS("data/dat_normumi_xpress.rds")

# XCI escapees from doi:10.1038/s41467-022-29414-1 
x.escapees <- c("1810030O07Rik","2010000l03Rik","2010308F09Rik","2610029G23Rik","5530601H04Rik","6720401G13Rik","Abcd1","Araf","Atp6ap2","BC022960","Bgn","Car5b","D330035K16Rik","D930009K15Rik","Ddx3x","Eif1ax","Eif2s3x","Fam50a","Flna","Ftsj1","Fundc1","Gdi1","Gemin8","Gpkow","Huwe1","Idh3g","Igbp1","Ikbkg","Kdm5c","Kdm6a","Lamp2","Maged1","Mbtps2","Med14","Mid1","Mmgt1","Mpp1","Msl3","Ndufb11","Nkap","Ogt","Pbdc1","Pdha1","Prdx4","Rbm3","Renbp","Sh3bgrl","Shroom4","Sms","Suv39h1","Syap1","Tbc1d25","Timp1","Trap1a","Uba1","Usp9x","Utp14a","Uxt","Xist","Yipf6")
dat.normumi.xpress[, x.escapee := GeneID %in% as.data.table(rtracklayer::import("out_mRNA/clonal_xpress/clonal_xpress.final_annot.gtf.gz", feature.type = "gene"))[match(x.escapees, gene_name), na.omit(gene_id)]]

#

dat.normumi.xpress[tpm > 1 & !is.na(segment) & x.escapee == F & note != "Exclude", chraratio := tpm/median(tpm[chrom %in% 1:19 & copies == ploidy],na.rm=T), by=c("Var2")]
dat.normumi.ratio.xpress <- dat.normumi.xpress[!is.na(chraratio) & note != "Exclude", list(V1 = median(chraratio), w = .N), by=c("Var2", "chrom", "clone", "ploidy", "copies", "segment")]
dat.normumi.ratio.xpress[chrom %in% c(1:19,"X"), chra := ifelse(chrom == "X", "chrX", ifelse(copies != ploidy, paste0("chr", chrom), "chrA"))]

dat.normumi.xpress[chrom %in% c(1:19,"X") & note != "Exclude", chra := ifelse(chrom == "X", "chrX", ifelse(copies != ploidy, paste0("chr", chrom), "chrA"))]
dat.normumi.xpress[chrom %in% c(1:19,"X"), chrx := chrom == "X"]
dat.normumi.ratiomerge.xpress <- dat.normumi.xpress[!is.na(chraratio), median(chraratio), by=c("Var2", "chra", "clone", "ploidy", "copies")]

dat.normumi.avgsegx.xpress <- melt(dat.normumi.xpress[normumi > 0 & !is.na(segment) & x.escapee == F, list("total"=base::mean(normumi, trim=0.2, na.rm=T), "c57"=base::mean(c57, trim=0.2, na.rm=T), "cast"=base::mean(cast, trim=0.2, na.rm=T)), by=c("chrom", "segment", "Var2", "clone", "ploidy", "copies")], measure.vars=c("total", "c57", "cast"))

## PLOT
library(ggplot2)
library(cowplot)

# degree of compensation
dat.normumi.ratiomerge.xpress[ploidy == 2 & copies <= ploidy, median(V1)/0.5, by=c("chra")]

# compare chrX
dat.normumi.ratiomerge.xpress[ploidy == 2 & copies <= ploidy & chra == "chrX", pairwise.wilcox.test(V1/0.5, clone, p.adjust.method = "fdr")]

# merged
p.ratio.merge.xpress <-
ggplot(dat.normumi.ratio.xpress[ploidy == 2 & !is.na(chra) & copies <= ploidy, median(rep(V1,w)), by=c("Var2", "clone","copies", "chra")], aes(x=factor(chra, levels=c("chrA", "chr4", "chr7", "chr17", "chrX")), y=V1 )) +
  geom_hline(yintercept=1, col="grey50", lty=2) +
  geom_violin(scale="width", show.legend = F,aes(fill=factor(copies))) +
  geom_boxplot(outlier.stroke = NA, width=0.25, aes(fill=factor(copies))) +
  geom_text(data=dat.normumi.ratio.xpress[ploidy == 2 & !is.na(chra) & copies <= ploidy, median(rep(V1,w)), by=c("Var2", "clone","copies", "chra")][, reshape2::melt(as.list(pairwise.wilcox.test(V1, chra, p.adjust.method="fdr")[[3]]["chrX",])), by="clone"], aes(x=L1, y=1.5, label=signif(value,2))) +
  facet_grid(~clone, scales="free_x", space="free_x") +
  coord_cartesian(ylim=c(0,2)) +
  labs(x=NULL, y="Autosomal expression ratio") +
  scale_fill_manual(values = c("#93C5DE", "#F8F7F7"), name="Copies") +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"))

ggsave2("plots/xpress_chraratio_merge.pdf", width = 5, height = 3, p.ratio.merge.xpress)

##
