suppressMessages(source("R/preprocess_data_fibroblast_ss2.R"))

##

load.kinetics.fibro <- function(file){
  require(data.table)
  dat <- fread(file)
  names(dat)[1] <- "gene"
  dat.melt <- melt(dat, id.vars="gene")
  dat.melt[as.data.table(refseq.gene.fibro), chr := chrom, on = "gene == name2"] # add chromosome annotations
  dat.melt[, c("allele","clone", "allele_loss") := tstrsplit(variable, "_", fixed=TRUE)] # split sample annotations
  return(dat.melt)
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

dat.bf.fibro <- load.kinetics.fibro("kinetics/output/split/fibroblast_bf.tsv")
dat.bs.fibro <- load.kinetics.fibro("kinetics/output/split/fibroblast_bs.tsv")

#
library(ggplot2)
library(cowplot)
p.bf <- 
ggplot(dat.bf.fibro[chr %in% c("chr3", "chr6")], aes(x=allele_loss, y=value, col=allele)) +
  stat_summary(fun.data="median_cl_boot", position=position_dodge(0.9)) +
  facet_grid(chr~clone, scale="free_x", space="free_x") +
  labs(x=NULL, y="Burst frequency") +
  scale_color_brewer(palette="Set1",name="Allele") +
  expand_limits(y=0) +
  theme_cowplot() +
  theme(strip.background = element_blank(), legend.position = "top")

p.bs <- 
ggplot(dat.bs.fibro[chr %in% c("chr3", "chr6")], aes(x=allele_loss, y=value, col=allele)) +
  stat_summary(fun.data="median_cl_boot", position=position_dodge(0.9)) +
  facet_grid(chr~clone, scale="free_x", space="free_x") +
  labs(x=NULL, y="Burst size") +
  scale_color_brewer(palette="Set1",name="Allele") +
  expand_limits(y=0) +
  theme_cowplot() +
  theme(strip.background = element_blank(), legend.position = "top")

ggsave2("plots/fibroblast_kinetics.pdf",width = 4,height = 4,
  plot_grid(p.bf, p.bs)
)

