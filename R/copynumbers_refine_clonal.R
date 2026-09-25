library(data.table)

## Read data and corrected copy numbers
cnt.depth.seg <- fread("data/dat_clonal_cnv_fixed.tsv")
cnt.depth <- fread("data/dat_clonal_dnaseq.tsv")

# calculate aneuploidy score
cnt.depth.seg[, aneuploidy_score := sum(abs(round(ploidy[chr %in% 1:19])-round(copies[chr %in% 1:19])),na.rm=T), by="variable"] # +1 per segment deviating from ploidy

# add corrected copy number information and rescale segments
cnt.depth[cnt.depth.seg, c("ploidy", "copies") := list(ploidy, copies), on = c("variable == variable", "chr == chr", "segment == segment")]
cnt.depth[, copynumbers := adjustedcopynumbers*(copies/adjustedsegments)]
cnt.depth[variable == "MEF_Bm4" & chr == "Y", copynumbers := adjustedcopynumbers*(copies/median(adjustedcopynumbers))] # fix chrY
cnt.depth[cnt.depth.seg, aneuploidy_score := aneuploidy_score, on = c("variable == variable")]

fwrite(cnt.depth[, list(variable, chr, start,end,segment,ploidy,copies,copynumbers,ratio, aneuploidy_score)], "data/dat_clonal_dnaseq_fixed.tsv", sep = "\t", quote = F)
fwrite(cnt.depth.seg[, list(variable, chr, start,end,segment,ploidy,aneuploidy_score,ratio,copies,note)], "data/dat_clonal_cnv_fixed.tsv", sep = "\t", quote = F)

##