library(UMIcountR)

##

get.spikes <- function(x, spike_info = spike_info, max_pattern_dist = 2, nBCs = 5, N = 100){
  require(UMIcountR)
  require(data.table)
  
  spikedat <- extract_complex_spike_dat(x, bc_df = spike_info, max_pattern_dist = 2)
  spikedat[, UB_directional := return_corrected_umi(UX, editham = 1, collapse_mode = "adjacency_directional"), by = BC]
  
  spikedat.n <- spikedat[, .(.N, nBCs = length(unique(BC))), by= "spikeUMI_hd2"]
  
  # filter and transform to matrix
  spikedat.cast <- dcast(spikedat, spikeBC_correct~BC, fun.aggregate = function(x) length(unique(x)), value.var = "spikeUMI_hd2")
  spikedat.filt.cast <- dcast(spikedat[!spikeUMI_hd2 %in% spikedat.n[nBCs > nBCs | N > N, spikeUMI_hd2]], spikeBC_correct~BC, fun.aggregate = function(x) length(unique(x)), value.var = "spikeUMI_hd2")
  
  out <- list(
    "spike_data" = spikedat,
    "spike_occurrence" = spikedat.n,
    "spike_matrixfull" = spikedat.cast,
    "spike_matrixfilt" = spikedat.filt.cast
  )
  return(out)
}

##

## get spikeUMI data
load("bin/UMIcountR/data/molspike_barcodes_infos_fivep_final.RData")

# index bam file with samtools before running
ls.spike.batch1 <- get.spikes("out_mRNA/clonal_xpress/clonal_xpress.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam")
ls.spike.batch2 <- get.spikes("out_mRNA/clonal_xpress_batch2/clonal_xpress_batch2.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam")

# export tables
fwrite(ls.spike.batch1$spike_data, "out_mRNA/clonal_xpress/zUMIs_output/molspikes/clonal_xpress.spikes.data.txt.gz", quote = F, sep = "\t")
fwrite(ls.spike.batch1$spike_occurrence, "out_mRNA/clonal_xpress/zUMIs_output/molspikes/clonal_xpress.spikes.occurence.txt.gz", quote = F, sep = "\t")
fwrite(ls.spike.batch1$spike_matrixfull, "out_mRNA/clonal_xpress/zUMIs_output/molspikes/clonal_xpress.spikes.full.txt.gz", quote = F, sep = "\t")
fwrite(ls.spike.batch1$spike_matrixfilt, "out_mRNA/clonal_xpress/zUMIs_output/molspikes/clonal_xpress.spikes.filtered.txt.gz", quote = F, sep = "\t")

fwrite(ls.spike.batch2$spike_data, "out_mRNA/clonal_xpress_batch2/zUMIs_output/molspikes/clonal_xpress_batch2.spikes.data.txt.gz", quote = F, sep = "\t")
fwrite(ls.spike.batch2$spike_occurrence, "out_mRNA/clonal_xpress_batch2/zUMIs_output/molspikes/clonal_xpress.spikes.occurence.txt.gz", quote = F, sep = "\t")
fwrite(ls.spike.batch2$spike_matrixfull, "out_mRNA/clonal_xpress_batch2/zUMIs_output/molspikes/clonal_xpress.spikes.full.txt.gz", quote = F, sep = "\t")
fwrite(ls.spike.batch2$spike_matrixfilt, "out_mRNA/clonal_xpress_batch2/zUMIs_output/molspikes/clonal_xpress.spikes.filtered.txt.gz", quote = F, sep = "\t")


## 834147