suppressMessages(source("R/preprocess_data_fibroblast_ss2.R"))

##

reformat.zlm <- function(fit, hypothesis){
  require(MAST)
  res <- summary(fit, doLRT=hypothesis)
  res.dt <- res$datatable
  hurdle <- merge(res.dt[contrast==hypothesis & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  res.dt[contrast==hypothesis & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  hurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  setnames(hurdle, c("primerid","coef","Pr(>Chisq)"), c("gene","logFC","p-value"))
  setorder(hurdle,fdr)
  return(hurdle)
}

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

# add gene position
tpm.melt.fibro[as.data.table(refseq.gene.fibro), pos := txStart, on= "gene == name2"]
tpm.melt.all.fibro[as.data.table(refseq.gene.fibro), pos := txStart, on= "gene == name2"]

# average data
tpm.melt.avg.fibro <- tpm.melt.fibro[expressed == T & !is.na(chr), base::mean(value, trim=0.2,na.rm=T),by=c("sample","clone_id", "sample_type","l1","chr")]
tpm.melt.avg.fibro[,refratio := V1[l1 == "c57"] / (V1[l1 == "c57"] + V1[l1 == "cast"]), by=c("sample", "chr")]
tpm.melt.avg.fibro[,subgroup := unique(refratio[chr == "chr3"]) < 0.25, by="sample"]

tpm.melt.all.avg.fibro <- tpm.melt.all.fibro[expressed == T & !is.na(chr), mean(value,na.rm=T),by=c("sample","clone_id","sample_type","sex","cross","chr")]
tpm.melt.all.avg.fibro[unique(tpm.melt.avg.fibro[,list(sample,chr,subgroup)]), subgroup := subgroup, on = c("sample == sample", "chr == chr")]

tpm.melt.fibro[unique(tpm.melt.avg.fibro[,list(sample,subgroup)]), subgroup := subgroup, on = c("sample == sample")]

## PLOT
library(ggplot2)
library(cowplot)
library(ggbeeswarm)

p.refratio.fibro <-
ggplot(refratio.melt.fibro[expressed == T & det > 0.25 & chr %in% c("chr3","chr6") & sample_type == "single.fibroblast" & !is.na(value)], aes(y=reorder(sample, value ), x=gene, fill=value)) +
  ggrastr::geom_tile_rast() +
  facet_grid(clone_id~chr, scales="free", space="free") +
  scale_fill_distiller(palette="Spectral", name="C57") +
  theme_void() +
  theme(panel.border = element_rect(colour="black", fill=NA))

p.allelic.fibro <-
ggplot(tpm.melt.avg.fibro[chr %in% paste0("chr",c(3,6)) & sample_type == "single.fibroblast"], aes(x=reorder(sample,-refratio), y=V1, col=l1 )) +
  geom_point(stroke=NA) +
  geom_hline(data=tpm.melt.avg.fibro[chr %in% paste0("chr",c(3,6)) & sample_type == "single.fibroblast" & subgroup == T, mean_sdl(V1, mult=2), by="chr"], aes(yintercept=ymin)) +
  geom_hline(data=tpm.melt.avg.fibro[chr %in% paste0("chr",c(3,6)) & sample_type == "single.fibroblast" & subgroup == T, mean_sdl(V1, mult=2), by="chr"], aes(yintercept=ymax)) +
  facet_grid(chr~clone_id, scales = "free_x") +
  labs(x="Cells", y="Expression (TPM)") +
  coord_cartesian(ylim = c(0,40) ) +
  scale_color_brewer(palette="Set1", name="Allele") +
  theme_cowplot() +
  theme(strip.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.border = element_rect(colour="black") )

p.genomic.fibro <-
ggplot(tpm.melt.fibro[clone_id == "fibroblast.cl7" & subgroup == T & !is.na(value) & chr %in% c("chr3", "chr6", "chr2", "chr7"), median(value, na.rm=T), by=c("l1", "chr", "pos")], aes(x=pos/1e6, y=log2(V1+1), col=l1, fill=l1)) + 
  ggrastr::geom_point_rast(alpha=0.1, stroke=NA) +
  geom_smooth(method="lm") +
  facet_grid(~chr, scales="free_x", space="free_x") + 
  labs(x="Mbp", y="TPM (log2+1)") +
  coord_cartesian(ylim=c(0,5)) +
  scale_color_brewer(palette="Set1", name="Allele", aesthetics = c("col", "fill")) +
  theme_cowplot() +
  theme(strip.background = element_blank())


## EXPORT
ggsave2("plots/fibroblast_refratio_heatmap.pdf", width = 4, height = 4, p.refratio.fibro)
ggsave2("plots/fibroblast_allelic_point.pdf", width = 8, height = 4, p.allelic.fibro)
ggsave2("plots/fibroblast_genomic_point.pdf", width = 5, height = 4, p.genomic.fibro)
