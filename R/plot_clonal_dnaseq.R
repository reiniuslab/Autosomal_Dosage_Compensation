library(data.table)

## Load data
cnt.depth <- fread("data/dat_clonal_dnaseq_fixed.tsv")
cnt.depth.seg <- fread("data/dat_clonal_cnv_fixed.tsv")

# median per segment 
cnt.depth[, segment_median := median(copynumbers, na.rm=T), by=c("variable", "chr", "segment")]
cnt.depth[, segment_ratio := median(ratio, na.rm=T), by=c("variable", "chr", "segment")]

# aneuploidy score per ploidy
cnt.depth.seg[!duplicated(variable), median(aneuploidy_score), by="ploidy"]

## PLOT
library(ggplot2)
library(cowplot)
library(patchwork)
## Heatmap
p.dna.heat <-
ggplot(cnt.depth[chr != "Y" & !grepl("new_", variable)], aes(x=factor(start), y=reorder(variable,-aneuploidy_score), fill= pmax(copynumbers,0)-round(ploidy) )) + 
  #geom_tile() +
  ggrastr::geom_tile_rast() +
  facet_grid(ploidy~factor(chr, levels=c(1:19,"X")), scales="free", space="free") +
  scale_fill_distiller(palette="RdBu", limits=c(-2,2), oob=scales::squish, name="Copies", na.value = "grey50") +
  scale_y_discrete(expand=c(0,0)) +
  theme_cowplot() +
  theme(strip.background = element_blank(), legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_blank(), panel.border = element_rect(colour="black"), panel.spacing = unit(0, "lines") )

p.dna.heat2 <-
ggplot(cnt.depth[chr != "Y" & !grepl("new_", variable)], aes(x=factor(start), y=reorder(variable,-aneuploidy_score), fill= pmax(copynumbers,0) )) + 
  #geom_tile() +
  ggrastr::geom_tile_rast() +
  facet_grid(ploidy~factor(chr, levels=c(1:19,"X")), scales="free", space="free") +
  scale_fill_viridis_c(oob=scales::squish, limits=c(1,6), name="Copies", na.value = "grey50") +
  scale_y_discrete(expand=c(0,0)) +
  theme_cowplot() +
  theme(strip.background = element_blank(), legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_blank(), panel.border = element_rect(colour="black"), panel.spacing = unit(0, "lines") )

p.dna.aneu <-
ggplot(cnt.depth.seg[!duplicated(variable) & !grepl("new_", variable)], aes(y=reorder(variable, -aneuploidy_score), x=aneuploidy_score)) +
  geom_col() +
  facet_grid(ploidy~., scales="free", space="free") +
  labs(x="Aneuploidy score", y=NULL) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), panel.spacing = unit(0, "lines") )

# chromosome plot
p.dna.cnv <-
ggplot(cnt.depth[chr != "Y"], aes(x=start, y=copynumbers)) + 
  ggrastr::geom_point_rast(stroke=NA, size=0.75) + 
  geom_segment(data=cnt.depth[chr != "Y", list("start"=min(start), "end"=max(end)), by=c("variable", "chr", "segment_median")], col="#FF7F00", aes(x=start, xend=end, y=segment_median)) +
  geom_text(data=cnt.depth.seg[note != ""], col="red", hjust=0, aes(label=note, x=start, y=copies)) +
  facet_grid(variable~factor(chr, levels=c(1:19,"X")), scales="free_x") +
  labs(y="DNA copies") +
  coord_cartesian(ylim=c(0,6.5)) +
  scale_y_continuous(breaks = 0:6) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.y = element_line(colour="grey50", linetype=2), )

p.dna.baf <- 
ggplot(cnt.depth[chr != "Y"], aes(x=start, y=ratio)) +
  ggrastr::geom_point_rast(stroke=NA, size=0.75) +
  geom_segment(data=cnt.depth[chr != "Y", list("start"=min(start), "end"=max(end)), by=c("variable", "chr", "segment_ratio")], col="#FF7F00", aes(x=start, xend=end, y=segment_ratio)) +
  labs(x="Chromosome Position (Mbp)", y="Allele ratio") +
  coord_cartesian(ylim=c(0,1)) + 
  facet_grid(variable~factor(chr, levels=c(1:19,"X")), scales="free_x") + 
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank() )

## EXPORT
ggsave2("plots/clonal_dnaseq_heatmap.pdf", width = 12, height = 6, p.dna.heat+p.dna.aneu+plot_layout(widths = c(1,0.3)))
ggsave2("plots/clonal_dnaseq_heatmap_absolute.pdf", width = 12, height = 6, p.dna.heat2+p.dna.aneu+plot_layout(widths = c(1,0.3)))
ggsave2("plots/clonal_dnaseq_cnv.pdf", width = 14, height = 28, p.dna.cnv/p.dna.baf+plot_layout(heights = c(1,0.75)))

