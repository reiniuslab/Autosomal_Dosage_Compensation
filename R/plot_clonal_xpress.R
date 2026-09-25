library(data.table)

dat.normumi.xpress <- readRDS("data/dat_normumi_xpress.rds")

##

compare_slopes <- function(x,y,which=2){
  x.est <- coef(summary(x))[which, "Estimate"]
  y.est <- coef(summary(y))[which, "Estimate"]
  
  x.se <- coef(summary(x))[which, "Std. Error"]
  y.se <- coef(summary(y))[which, "Std. Error"]
  
  x.n <- nrow(x$model)
  y.n <- nrow(y$model)
  
  # Welch's T-test statistic
  t.stat <- (x.est - y.est) / sqrt(x.se^2 + y.se^2)
  
  # Degrees of freedom
  df <- ((x.se^2 + y.se^2)^2) / (((x.se^2)^2 / (x.n - 2)) + ((y.se^2)^2 / (y.n - 2)))
  
  # Calculate p-value
  p.value <- 2 * pt(-abs(t.stat), df)
  
  print(paste0(
    "x_slope=", round(x.est,3),
    " y_slope=", round(y.est,3),
    " p=", signif(p.value,3)
  ))
}

tost_slopes <- function(x,y,which=2,eqb){
  require(TOSTER)
  x.est <- coef(summary(x))[which, "Estimate"]
  y.est <- coef(summary(y))[which, "Estimate"]
  
  x.se <- coef(summary(x))[which, "Std. Error"]
  y.se <- coef(summary(y))[which, "Std. Error"]

  x.n <- nrow(x$model)
  y.n <- nrow(y$model)

  tsum_TOST(x.est, x.se^2, x.n, y.est, y.se^2, y.n, eqb=eqb)
}  

# summarise per region

dat.normumi.avgseg.xpress <- melt(dat.normumi.xpress[normumi > 0 & !is.na(segment), list("total"=base::mean(normumi, trim=0.2, na.rm=T), "c57"=base::mean(c57, trim=0.2, na.rm=T), "cast"=base::mean(cast, trim=0.2, na.rm=T)), by=c("chrom", "segment", "Var2", "clone", "ploidy", "copies")], measure.vars=c("total", "c57", "cast"))

dat.normumi.avgploidy.xpress <- dat.normumi.xpress[normumi > 0 & !is.na(segment) & chrom %in% 1:19 & note != "Exclude", list(base::mean(normumi, trim = 0.2, na.rm=T), .N), by=c("Var2", "clone", "ploidy", "copies")]
dat.normumi.avgploidy.xpress[, exp := mean(V1[copies == ploidy])*(copies/ploidy), by=c("ploidy")]

## PLOT
library(ggplot2)
library(cowplot)

# linear regression
dat.normumi.avgploidy.xpress[, summary(lm(V1~log2(as.numeric(copies))+clone ))]
dat.normumi.avgploidy.xpress[ploidy == 2, summary(lm(V1~log2(as.numeric(copies)) ))]
dat.normumi.avgploidy.xpress[ploidy == 4, summary(lm(V1~log2(as.numeric(copies)) ))]

# difference between observed and expected slopes
res.slope <- dat.normumi.avgploidy.xpress[, compare_slopes(lm(V1~log2(as.numeric(copies))), lm(exp~log2(as.numeric(copies))) ), by=c("ploidy")]

# equivalence test for 2N and 4N slopes within 0.5 units
dat.normumi.avgploidy.xpress[, tost_slopes(lm(V1[ploidy==2]~log2(as.numeric(copies[ploidy==2]))), lm(V1[ploidy==4]~log2(as.numeric(copies[ploidy==4]))), eqb = 0.5 )]

# buffering effect
dat.normumi.avgploidy.xpress[copies != ploidy, mean(2^abs(log2(V1/exp))), by=c("ploidy")]

# expression vs copies
p.normumi.copies.xpress <-
ggplot(dat.normumi.avgploidy.xpress, aes(x=factor(copies), y=V1 )) + 
  geom_boxplot(outlier.stroke = NA, outlier.colour = "#999999", show.legend = F, aes(fill=factor(copies-ploidy))) +
  geom_point(data=dat.normumi.avgploidy.xpress[copies != ploidy, median(exp), by = c("ploidy", "copies")], show.legend = F) +
  geom_text(data=res.slope, hjust=0, aes(y=14.5, x=0.5, label=gsub(" ", "\n", V1) )) +
  #geom_text(data=dat.normumi.avgploidy.xpress[copies != ploidy, t.test(V1, mu=exp[1])$p.value, by=c("ploidy", "copies")], nudge_x = -0.25, angle = 90, hjust = 1, aes(y=15, label=paste0("P = ",signif(V1,3) ) )) +
  labs(y="Expression (Spike-normalised UMIs)", x="DNA copies") +
  facet_grid(~ploidy, scales="free", space = "free_x") +
  expand_limits(y=0) +
  coord_cartesian(ylim=c(0,15)) +
  scale_fill_manual(values=c("#0573B2", "#93C5DE", "#F8F7F7", "#F5A582")) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"))

ggsave2("plots/xpress_exprs_copy.pdf", width = 2, height = 3, p.normumi.copies.xpress)

# expression vs copies per clone
p.normumi.copiesclone.xpress <-
ggplot(dat.normumi.avgploidy.xpress, aes(x=factor(copies), y=V1, fill=factor(copies-ploidy) )) + 
  geom_boxplot(outlier.stroke = NA, outlier.colour = "#999999", show.legend = F) +
  geom_point(data=dat.normumi.avgploidy.xpress[, list(mean(V1[copies == ploidy])*(copies/ploidy), copies), by = c("ploidy", "clone")][copies != ploidy, median(V1), by=c("ploidy", "clone", "copies")], show.legend = F) +
  labs(y="Expression (Spike-normalised UMIs)", x="DNA copies") +
  facet_grid(~ploidy+clone, scales="free", space = "free_x") +
  expand_limits(y=0) +
  coord_cartesian(ylim=c(0,25)) +
  scale_fill_manual(values=c("#0573B2", "#93C5DE", "#F8F7F7", "#F5A582", "#CA0020")) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"))

ggsave2("plots/xpress_exprs_copy_clone.pdf", width = 9, height = 3, p.normumi.copiesclone.xpress)

# expression per segment (rasterised)
p.normumi.segment.xpress <-
ggplot(dat.normumi.avgseg.xpress[variable != "total" & chrom != "Y"], aes(x=factor(segment), y=log2(value+1))) +
  ggrastr::geom_point_rast(stroke=NA, size=0.75, position = position_dodge2(0.85), aes(col=variable)) +
  geom_boxplot(outlier.shape = NA, width=0.5, position = position_dodge(0.85), aes(group=paste(variable, segment) )) + 
  facet_wrap(ploidy+clone~chrom, scales="free_x", ncol=20) +
  labs(x="Chromosome segment", y="Spike-normalised UMIs (log2+1)", lty=NULL) +
  coord_cartesian(ylim=c(0,2)) +
  scale_color_brewer(palette="Set1", name="Allele") +
  theme_cowplot() +
  theme(strip.background = element_blank(), axis.text.x = element_blank(), panel.border = element_rect(colour="black"))

ggsave2("plots/xpress_exprs_segment.pdf", width = 18, height = 20, p.normumi.segment.xpress)

# expression per segment (vectorised)
p.normumi.segment.xpress.v <-
  ggplot(dat.normumi.avgseg.xpress[variable != "total" & chrom != "Y"], aes(x=factor(segment), y=log2(value+1))) +
  geom_point(stroke=NA, size=0.75, position = position_dodge2(0.85), aes(col=variable)) +
  geom_boxplot(outlier.shape = NA, width=0.5, position = position_dodge(0.85), aes(group=paste(variable, segment) )) + 
  facet_wrap(ploidy+clone~chrom, scales="free_x", ncol=20) +
  labs(x="Chromosome segment", y="Spike-normalised UMIs (log2+1)", lty=NULL) +
  coord_cartesian(ylim=c(0,2)) +
  scale_color_brewer(palette="Set1", name="Allele") +
  theme_cowplot() +
  theme(strip.background = element_blank(), axis.text.x = element_blank(), panel.border = element_rect(colour="black"))

ggsave2("plots/xpress_exprs_segment_vectorised.pdf", width = 18, height = 20, p.normumi.segment.xpress.v)

##
