library(data.table)

##

# load data
dntr.rna <- readRDS("data/dntr_rna.rds")
dntr.dna <- readRDS("data/dntr_dna.rds")

# segments
dntr.dna.seg <- dntr.dna[!is.na(seg), list("start"=min(start), "end"=max(end), "copy"=median(seg_median,na.rm=T) ), by=c("chr","sample","clust", "seg") ]

# cluster dna-seq
mat.dna <- as.matrix(data.frame(dcast(unique(dntr.dna[chr %in% paste0("chr", 1:22), list(sample, chr, start, end, seg, seg_median)]), paste(chr,start,end,seg, sep=":")~sample, value.var="seg_median"), row.names = 1))
clus <- hclust(dist(t(mat.dna), method = "euclidean"), method= "complete")
grp <- cutree(clus, h=60)


dntr.dna[, clust := grp[sample]]
dntr.rna[, clust := grp[sample]]

#
dntr.rna[, rel.copy := round((2*2^seg_median)-2)]
dntr.rna[normcounts > 0 & chr %in% paste0("chr",1:22), relative := normcounts / median(normcounts[rel.copy == 0], na.rm=T), by=c("sample")]

dntr.rna.avg <- dntr.rna[normcounts > 0 & chr %in% paste0("chr",1:22), base::mean(relative, na.rm=T, trim=0.2), by=c("sample", "treatment", "timepoint", "rel.copy")]
dntr.rna.avg[, N := .N, by=c("treatment", "rel.copy")]
dntr.rna.avg[!is.na(rel.copy) & N > 50, exp := mean(V1[rel.copy == 0], na.rm=T) * (rel.copy+2)/2, by=c("treatment")]

dntr.rna.avgseg <- dntr.rna[normcounts > 0 & chr %in% paste0("chr",1:22) & !is.na(seg), base::mean(normcounts, na.rm=T, trim=0.2), by=c("sample", "clust", "chr", "seg", "rel.copy")]
dntr.rna.avgseg[, N := .N, by=c("chr", "clust", "seg")]
dntr.rna.avgseg[N > 10, exp := mean(V1[rel.copy == 0], na.rm=T) * (rel.copy+2)/2]

## PLOT
library(ggplot2)
library(cowplot)

p.dntr.heat <- 
ggplot(dntr.dna[chr %in% paste0("chr",1:22)], aes(x=factor(start), y=reorder(sample,clust), fill=(2*2^seg_median)-2 )) +
  ggrastr::geom_tile_rast() + 
  facet_grid(factor(clust)~factor(gsub("chr","",chr), levels=1:22), scales="free", space="free") +
  scale_fill_fermenter(palette="RdBu", breaks=seq(-1.5,1.5,1), limits=c(-2,2), oob=scales::squish, name="Copies") +
  labs(x="Position along chromosome", y="Cells") +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), panel.spacing = unit(0, "mm"), axis.text = element_blank(), axis.ticks = element_blank())

ggsave2("plots/dntr_dna_heatmap.pdf", width = 6, height = 4, p.dntr.heat)

p.dntr.exprs.copy <-
ggplot(dntr.rna.avg[!is.na(rel.copy) & N > 50], aes(x=factor(rel.copy), y=V1, fill=factor(rel.copy))) +
  geom_boxplot(outlier.stroke = NA, outlier.colour = "#999999", show.legend = F) +
  geom_point(data=dntr.rna.avg[!is.na(rel.copy) & N > 50, median(exp, na.rm=T), by=c("treatment", "rel.copy")], show.legend = F) +
  facet_grid(~gsub("","",treatment), space="free_x", scales="free_x") +
  labs(x="DNA copies", y="Relative expression") +
  scale_fill_manual(values=c("#93C5DE", "#F8F7F7", "#F5A582")) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"))

ggsave2("plots/dntr_rna_copies.pdf", width = 6, height = 4, p.dntr.exprs.copy)

#

p.dntr.rna <-
ggplot(dntr.rna.avgseg[(chr %in% paste0("chr", c(11,12,14)) & clust %in% 13 & N > 10) | (chr %in% paste0("chr", c(1:3)) & clust %in% 10 & N > 10) | (chr %in% paste0("chr", c(5:7)) & clust %in% 16 & N > 10) | (chr %in% paste0("chr", c(1:3)) & clust %in% 30 & N > 10)], aes(x=factor(seg), y=V1)) +
  ggbeeswarm::geom_quasirandom(stroke=NA, col="#999999") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(data=dntr.rna.avgseg[(chr %in% paste0("chr", c(11,12,14)) & clust %in% 13 & N > 10) | (chr %in% paste0("chr", c(1:3)) & clust %in% 10 & N > 10) | (chr %in% paste0("chr", c(5:7)) & clust %in% 16 & N > 10) | (chr %in% paste0("chr", c(1:3)) & clust %in% 30 & N > 10), median(exp), by=c("chr", "clust", "seg")]) +
  facet_grid(~clust+factor(chr,levels=paste0("chr",1:22)), scales="free_x", space="free_x") +
  labs(x=NULL, y="Expression (ERCC-normalised counts)") +
  coord_cartesian(ylim=c(0,30)) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"))

p.dntr.dna <-
ggplot(dntr.dna[(chr %in% paste0("chr",c(11,12,14)) & clust == 13) | (chr %in% paste0("chr",c(1:3)) & clust == 10) | (chr %in% paste0("chr",c(5:7)) & clust == 16) | (chr %in% paste0("chr",c(1:3)) & clust == 30), median((2*2^copy)-2, na.rm=T), by=c("chr", "start", "clust")], aes(x=start, y=V1 )) +
  geom_point(stroke=NA) +
  geom_segment(data=dntr.dna[(chr %in% paste0("chr",c(11,12,14)) & clust == 13) | (chr %in% paste0("chr",c(1:3)) & clust == 10) | (chr %in% paste0("chr",c(5:7)) & clust == 16) | (chr %in% paste0("chr",c(1:3)) & clust == 30), list("start" = min(start), "end" = max(end), "copy" = median((2*2^seg_median)-2, na.rm=T) ), by=c("chr", "clust", "sample", "seg")][, list("start" = median(as.numeric(start)), "end" = median(as.numeric(end)), "copy" = median(copy), "N" = length(unique(sample))), by=c("chr", "clust", "seg")][N>10], col="#FF7F00", aes(x=start, xend=end, y=copy, yend=copy)) +
  facet_grid(~clust+factor(chr,levels=paste0("chr",1:22)), space="free_x", scales="free_x") +
  labs(x="Position", y="DNA copies (relative)") +
  coord_cartesian(ylim=c(-1.5,1.5)) +
  theme_cowplot() +
  theme(strip.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.border = element_rect(colour="black"))
  
p.dntr.rna.all <-
ggplot(dntr.rna.avgseg[chr %in% paste0("chr", c(8,10,20)) & N > 10], aes(x=factor(seg), y=V1)) +
  geom_jitter(stroke=NA, col="#999999") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(data=dntr.rna.avgseg[chr %in% paste0("chr", c(8,10,20)) & N > 10, median(exp), by=c("chr", "seg")]) +
  facet_grid(~"All"+factor(chr,levels=paste0("chr",1:22)), scales="free_x", space="free_x") +
  labs(x=NULL, y="Expression (ERCC-normalised counts)") +
  coord_cartesian(ylim=c(0,30)) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"))

p.dntr.dna.all <-
ggplot(dntr.dna[chr %in% paste0("chr",c(8,10,20)), median((2*2^copy)-2, na.rm=T), by=c("chr", "start")], aes(x=start, y=V1 )) +
  geom_point(stroke=NA) +
  geom_segment(data=dntr.dna[chr %in% paste0("chr",c(8,10,20)), list("start" = min(start), "end" = max(end), "copy" = median((2*2^seg_median)-2, na.rm=T) ), by=c("chr", "sample", "seg")][, list("start" = median(as.numeric(start)), "end" = median(as.numeric(end)), "copy" = median(copy), "N" = length(unique(sample))), by=c("chr", "seg")][N>100], col="#FF7F00", aes(x=start, xend=end, y=copy, yend=copy)) +
  facet_grid(~"All"+factor(chr,levels=paste0("chr",1:22)), space="free_x", scales="free_x") +
  labs(x="Position", y="DNA copies (relative)") +
  coord_cartesian(ylim=c(-1.5,1.5)) +
  theme_cowplot() +
  theme(strip.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.border = element_rect(colour="black"))


## EXPORT
library(patchwork)
ggsave2("plots/dntr_chr.pdf", width = 8, height = 5, p.dntr.rna / p.dntr.dna)
ggsave2("plots/dntr_chr_all.pdf", width = 4, height = 5, p.dntr.rna.all / p.dntr.dna.all)

##