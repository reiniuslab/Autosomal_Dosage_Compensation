#!/usr/bin/env Rscript
# packages ----------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
#suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(yaml))

Sys.time()

##
## Original code by Christoph Ziegenhain (https://github.com/sandberg-lab/Smart-seq3/blob/5d5938475039f5c98d0d94faf89db917f66fe8ac/allele_level_expression/get_variant_overlap_CAST.R)
## Modified by Antonio Lentini
##

print("allele level expression v1.2 (AL EDIT v3 2021-05-24)")

# number crunching function -----------------------------------------------

load_reads <- function(is_UMI, cellBCs, filename, return_map, read_layout){
  if( is_UMI ){
    cols_to_read <- c(2,4,5,6,7,8)
    colname_vec <- c("pos","cigar","seq","BC","UB","GeneID")
  }else{
    cols_to_read <- c(2,4,5,6,7)
    colname_vec <- c("pos","cigar","seq","BC","GeneID")
  }

  reads <- fread(file = filename,
                 sep = "\t",
                 header = F, fill = T,
                 select = cols_to_read, #read only necessary cls
                 col.names = colname_vec)[ BC %in% cellBCs ][ ! GeneID == "" ] #directly drop unnecessary rows

  if(return_map == FALSE){
    system(paste("pigz -p ",ncores,filename))
  }
  if(read_layout == "SE"){
    reads[, readID := paste0("r_",1:nrow(reads))]
  }else{
    reads[, readID := paste0("r_",1:nrow(reads))] #for now, keep PE reads as individual reads, because I did not make sure that we load proper pairs adjacent to each other!
  }
  return(reads)
}

variant_parsing <- function(reads, variant_positions, is_UMI){
  #parse all cigars to reference seq
  ops <- c("M", "=", "X")
  ranges_on_ref <- cigarRangesAlongReferenceSpace(reads$cigar, pos=reads$pos, ops=ops)
  ranges_on_query <- cigarRangesAlongQuerySpace(reads$cigar, ops=ops)
  gc(verbose = F)
  range_group <- togroup(PartitioningByWidth(ranges_on_ref))
  ranges_on_ref <- unlist(ranges_on_ref, use.names=FALSE)
  ranges_on_query <- unlist(ranges_on_query, use.names=FALSE)
  query2ref_shift <- start(ranges_on_ref) - start(ranges_on_query)

  var_pos <- variant_positions
  hits <- findOverlaps(var_pos, ranges_on_ref)
  hits_at_in_x <- var_pos[queryHits(hits)] - query2ref_shift[subjectHits(hits)]
  hits_group <- range_group[subjectHits(hits)]
  fetched_bases <- subseq(reads[hits_group,]$seq, start=hits_at_in_x, width=1L)

  #now add everything together in the output data.table
  out_vars <- data.table(
    obs_base = fetched_bases,
    pos = var_pos[queryHits(hits)]
  )
  out_vars[, c("BC","GeneID","readID") := reads[hits_group, c("BC", "GeneID", "readID"), with = F] ]
  if( is_UMI ){
    out_vars[, UB := reads[hits_group]$UB ]
  }

  out_vars <- out_vars[obs_base %in% c("A","C","G","T") ]
  setnames(out_vars,"pos","POS")

  return(out_vars)
}

calc_coverage_new_return_map <- function(vcf_chunk, out, cellBCs, type, read_layout){
  chr <- unique(vcf_chunk$CHROM)
  is_UMI <- any(grepl("UMIs", type))
  print(paste("Starting to read data for chr ", chr))
  Sys.time()
  
  if(!is.null(nreads)){
    reads <- load_reads(is_UMI = is_UMI, cellBCs = cellBCs, filename = paste0(out,chr,".var_overlap.readsout_downsampled_",nreads), return_map = TRUE, read_layout = read_layout)
  }else{
    reads <- load_reads(is_UMI = is_UMI, cellBCs = cellBCs, filename = paste0(out,chr,".var_overlap.readsout"), return_map = TRUE, read_layout = read_layout)
  }

  print("Reading complete, processing reads & cigar values...")
  Sys.time()

  out_vars <- variant_parsing(reads, variant_positions = as.integer(vcf_chunk$POS), is_UMI = is_UMI)

  #crunch the numbers :-)
  out_vars <- merge(out_vars,vcf_chunk,by = "POS" )

  out_vars[          , basecall := "other"][
      obs_base == REF, basecall := "c57"][
      obs_base == ALT, basecall := "cast"]

  out_reads <- out_vars[, .(readcall = read_decision(basecall)), by = c("BC","GeneID","readID")]
  if( is_UMI ){
    out_UMIs <- out_vars[! UB == "" , .(UMIcall = read_decision(basecall)), by = c("BC","GeneID","UB")]
    return(out_UMIs)
  }else{
    return(out_reads)
  }
}

calc_coverage_new <- function(vcf_chunk, out, cellBCs, type, read_layout){
  chr <- unique(vcf_chunk$CHROM)
  is_UMI <- any(grepl("UMIs", type))
  print(paste("Starting to read data for chr ", chr))
  Sys.time()
  
  # AL EDIT: add subsampling
  if(!is.null(nreads)){
    reads <- load_reads(is_UMI = is_UMI, cellBCs = cellBCs, filename = paste0(out,chr,".var_overlap.readsout_downsampled_",nreads), return_map = FALSE, read_layout = read_layout)
  }else{
    reads <- load_reads(is_UMI = is_UMI, cellBCs = cellBCs, filename = paste0(out,chr,".var_overlap.readsout"), return_map = FALSE, read_layout = read_layout)
  }
  

  print("Reading complete, processing reads & cigar values...")
  Sys.time()

  out_vars <- variant_parsing(reads, variant_positions = as.integer(vcf_chunk$POS), is_UMI = is_UMI)

  #crunch the numbers :-)
  out_vars <- merge(out_vars,vcf_chunk,by = "POS" )

  out_vars[              , basecall := "other"][
          obs_base == REF, basecall := "c57"][
          obs_base == ALT, basecall := "cast"]

  # AL EDIT: check if reads are mapped, otherwise return empty tables
  if( nrow(out_vars) == 0 ){
    out_dat <- data.table(chr = character(), BC = character(), GeneID = character(), c57 = numeric(), cast = numeric(), other = numeric(), total = numeric(), CAST_fraction = numeric() )
    
    if ( is_UMI ){
      out_list <- list(reads = out_dat, UMIs = out_dat)
      return(out_list)
    }else{
      return(out_dat)
    }

  }

  out_reads <- out_vars[, .(readcall = read_decision(basecall)), by = c("BC","GeneID","readID")]

  if( is_UMI ){
    # AL EDIT: check if UMIs are mapped, otherwise return empty table
    if( out_vars[, sum(! UB == "") == 0 ] ){
      skip_UMIs <- TRUE
      out_dat_UMIs <- data.table(chr = character(), BC = character(), GeneID = character(), c57 = numeric(), cast = numeric(), other = numeric(), total = numeric(), CAST_fraction = numeric() )
    }else{
      skip_UMIs <- FALSE
      out_UMIs <- out_vars[! UB == "" , .(UMIcall = read_decision(basecall)), by = c("BC","GeneID","UB")]
    }
  }

  rm(out_vars)

  out_dat <- out_reads[
    , .N, by=.(BC,GeneID,readcall)][
      , chr := chr]

  rm(out_reads)

  out_dat[, readcall := factor(readcall, levels=c("c57", "cast", "other"))] # AL EDIT: This fixes issues where only one allele is detected
  
  out_dat <- dcast(out_dat, formula = chr+BC+GeneID ~ readcall, value.var = "N", fill = 0, drop = F) # AL EDIT: drop = F to keep unused alleles  out_dat[, total := c57+cast+other]
  out_dat[, total := c57+cast+other]

  out_dat <- out_dat[other/total < 0.33]

  out_dat[, CAST_fraction := cast/(cast+c57), by = c("BC","GeneID")]

  print("Done!")

  if( is_UMI ){
    # AL EDIT: check if UMI counting should be skipped
    if( !skip_UMIs ){
      out_dat_UMIs <- out_UMIs[
      , .N, by=c("BC","GeneID","UMIcall")][
        , chr := chr]
    
      out_dat_UMIs[, UMIcall := factor(UMIcall, levels=c("c57", "cast", "other"))] # AL EDIT: This fixes issues where only one allele is detected
    
      rm(out_UMIs)
    
      out_dat_UMIs <- dcast(out_dat_UMIs, formula = chr+BC+GeneID ~ UMIcall, value.var = "N", fill = 0, drop = F) # AL EDIT: drop = F to keep unused alleles
      out_dat_UMIs[, total := c57+cast+other]
    
      out_dat_UMIs <- out_dat_UMIs[other/total < 0.33]
    
      out_dat_UMIs[, CAST_fraction := cast/(cast+c57), by = c("BC","GeneID")]
    }

    out_list <- list(reads = out_dat,
                     UMIs = out_dat_UMIs)
    return(out_list)
  }else{
    return(out_dat)
  }
}

makeWide <- function(allele_dat, metric = c("cast","c57","CAST_fraction")){
  dat <- allele_dat[, c("BC","GeneID",metric), with = F]
  fill_val <- ifelse(metric %in% c("cast","c57"), 0, NA)
  dat_w <- dcast(dat, formula = GeneID ~ BC, fill=fill_val, value.var = metric)
  return(dat_w)
}

makeUMIs <- function(dge_path, CASTfracts){
  dge <- readRDS(dge_path)
  
  # AL EDIT: add subsampling
  if(!is.null(nreads)){
    ex <- as.matrix(dge$umicount$exon$downsampling[[paste0("downsampled_",nreads)]])
  }else{
    ex <- as.matrix(dge$umicount$exon$all)
  }
  
  fract_mat <- as.matrix(CASTfracts)
  row.names(fract_mat) <- fract_mat[,1]
  fract_mat <- fract_mat[,-1]
  class(fract_mat) <- "numeric"

  shared_genes <- intersect(row.names(fract_mat),row.names(ex))
  shared_cells <- intersect(colnames(fract_mat),colnames(ex))

  fract_mat <- fract_mat[shared_genes,shared_cells]
  ex <- ex[shared_genes,shared_cells]

  no_expr <- (ex == 0)
  umis_CAST <- round(fract_mat*ex,0)
  umis_BL6 <- round((1-fract_mat)*ex,0)

  umis_CAST[no_expr] <- 0
  umis_BL6[no_expr] <- 0

  outlist <- list(
    umis_CAST = umis_CAST,
    umis_BL6 = umis_BL6
  )

  return(outlist)
}

read_decision <- function(basecalls){
  if(length(basecalls) == 1){
    return(basecalls)
  }else{
    ux <- unique(basecalls)
    basecall_summary <- tabulate(match(basecalls, ux))
    names(basecall_summary) <- ux
    majority_basecall <- ux[which.max(basecall_summary)]
    if(basecall_summary[majority_basecall]/sum(basecall_summary) >= 0.66){
      return(majority_basecall)
    }else{
      return("other")
    }
  }
}

check_nonUMIcollapse <- function(seqfiles){
  #decide wether to run in UMI or no-UMI mode
  UMI_check <- lapply(seqfiles,
                      function(x) {
                        if(!is.null(x$base_definition)) {
                          if(any(grepl("^UMI",x$base_definition))) return("UMI method detected.")
                        }
                      })

  umi_decision <- ifelse(length(unlist(UMI_check))>0,"UMI","nonUMI")
  return(umi_decision)
}


# startup variables -------------------------------------------------------
# AL EDIT: fix formatting
option_list <- list(
  make_option(c("-y", "--yaml"), type="character",
              help="Coordinate sorted bam file. Mandatory"),
  make_option(c("-v", "--vcf"), type="character",
              help="SNP position list (VCF file) with variant annotation. Mandatory"),
  make_option(c("-t","--tagBC"), type="character",
                help="Bam tag containing cell barcodes. Default: BC", default="BC"),
  make_option(c("-m","--minCount"), type="integer",
                help="Cutoff for minimum coverage in a Cell/Gene pair. Default: 0", default=0),
  make_option(c("-u", "--umi_map"), action="store_true", default=FALSE,
              help="Print UMI-allele mapping table"),
  make_option(c("-s", "--subsample"), type="integer",
              help="Number of reads to subsample."), # AL EDIT: Add subsample argument
  make_option(c("-e", "--exclude"), type="character",
              help="Chromosomes to exclude, comma-separated.") # AL EDIT: Add chromosome exclusion argument
)
opt <- parse_args(OptionParser(option_list=option_list))

if (any(is.null(opt$yaml),is.null(opt$vcf))) {
  stop("All mandatory parameters must be provided. See script usage (--help)")
}


#####
#testing
#####
#BCtag <- "BC"
#path_snps <- "/home/chrisz/resources/genomes/Mouse/old_validated_cast_c57_snps.mm10.vcf"
#path_snps <- "/home/chrisz/resources/genomes/Mouse/CAST.SNPs.superset.vcf.gz"
#minC <- 0
#opt   <- read_yaml("/home/perj/moved_data/mmu/per_fibroblasts_final/zUMIs_rerun/zUMIs_rerun.yaml")
#outpath <- paste0(opt$out_dir,"/zUMIs_output/allelic/")
#####
#/testing
#####


BCtag <- opt$tagBC
path_snps <- opt$vcf
minC <- opt$minCount
map_flag <- opt$umi_map

# AL EDIT: add downsampling and chromosome exclusion
nreads <- opt$subsample
if(!is.null(nreads)){print(paste("Downsampling to",nreads,"reads"))}
if(!is.null(opt$exclude)){
  skip <- strsplit(opt$exclude,",")[[1]]
  print(paste("Skipping chromosomes:", skip))
}

opt   <- read_yaml(opt$yaml)
outpath <- paste0(opt$out_dir,"/zUMIs_output/allelic/")

  if(!dir.exists(outpath)){
    try(system(paste("mkdir",outpath)))
  }

  outpath <- paste0(outpath,opt$project,".")
  ncores <- opt$num_threads
  cellBCs <- paste0(opt$out_dir,"/zUMIs_output/",opt$project,"kept_barcodes.txt")

  setwd(opt$out_dir)
  setDTthreads(ncores)


  UMIdata_flag <- check_nonUMIcollapse(opt$sequence_files)


# read stuff --------------------------------------------------------------

cellBCs <- fread(cellBCs)
cellBCs <- cellBCs$XC

print("Reading Variants...")
if(grepl(path_snps, pattern = ".gz$")){
  vcf <- fread(cmd = paste("zcat",path_snps," | grep -v '^#'","| cut -f1,2,4,5"), col.names = c("CHROM","POS","REF","ALT"))
}else{
  vcf <- fread(cmd = paste("grep -v '^#'",path_snps,"| cut -f1,2,4,5"), col.names = c("CHROM","POS","REF","ALT"))
}

print("Done!")
Sys.time()

chroms_todo <- unique(vcf$CHROM)
if(!is.null(opt$exclude)){
  chroms_todo <- chroms_todo[! chroms_todo %in% skip] # AL EDIT: exclude chromosomes
}

# detect if zUMIs >= 2.6.0 is used ----------------------------------------
if( file.exists(paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.sorted.bam")) || file.exists(paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.UBcorrected.sorted.bam")) ){
  genetag <- "GE"
  if( file.exists(paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.UBcorrected.sorted.bam")) ){
    hammingflag <- TRUE
    path_bam <- paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.UBcorrected.sorted.bam")
  }else{
    hammingflag <- FALSE
    path_bam <- paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.sorted.bam")
  }
}else{
  genetag <- "XT"
  if( file.exists( paste0(opt$out_dir,"/",opt$project,".filtered.tagged.Aligned.out.bam.ex.featureCounts.UBfix.bam")) ){
    hammingflag <- TRUE
    path_bam <- file.exists( paste0(opt$out_dir,"/",opt$project,".filtered.tagged.Aligned.out.bam.ex.featureCounts.UBfix.bam"))
  }else{
    hammingflag <- FALSE
    path_bam <- paste0(opt$out_dir,"/",opt$project,".filtered.tagged.Aligned.out.bam.ex.featureCounts.bam")
  }
}

# extract unique maps per chromosome  -------------------------------------

# AL EDIT: check for downsampled files
if(
  is.null(nreads) & (file.exists( paste0(outpath,chroms_todo[[1]],".var_overlap.readsout") ) | file.exists( paste0(outpath,chroms_todo[[1]],".var_overlap.readsout.gz") )) |
  !is.null(nreads) & (file.exists( paste0(outpath,chroms_todo[[1]],".var_overlap.readsout_downsampled_",nreads) ) | file.exists( paste0(outpath,chroms_todo[[1]],".var_overlap.readsout_downsampled",nreads,".gz") )) 
  ){
  
  if(!is.null(nreads)){
    zipped_files <- list.files(path=paste0(opt$out_dir,"/zUMIs_output/allelic/"), pattern=paste0(".var_overlap.readsout_downsampled_",nreads,".gz"), full.names=T)
  }else{
    zipped_files <- list.files(path=paste0(opt$out_dir,"/zUMIs_output/allelic/"), pattern=".var_overlap.readsout.gz", full.names=T)
  }
  
  print("Decompressing reads...")
  for(f in zipped_files){
    system(paste("pigz -d -p",ncores,f))
  }
}else{
  print("Extracting reads...")
  samtoolsexc <- opt$samtools_exec
  if(UMIdata_flag == "UMI"){
    if(hammingflag){
      samtools_cmd1 <- "view -@2 -x QB -x QU -x ES -x IS -x EN -x IN -x GI -x BX -x UX -x NH -x AS -x nM -x HI -x IH -x NM -x uT -x MD -x jM -x jI -x XN -x XS -x vA -x vG -x vW"
      
      # AL EDIT: add downsampling
      if(!is.null(nreads)){
        samtools_cmd2 <- paste0(" | cut -f3,4,5,6,10,12,13,14 | shuf -n ",nreads," | grep '",genetag,"' | sed 's/",genetag,":Z://' | sed 's/UB:Z://' | sed 's/",BCtag,":Z://' | awk 'BEGIN{IFS=\"\t\";OFS=\"\t\";}{print $1,$2,$3,$4,$5,$6,$8,$7;}' | awk '{if($3 == \"255\"){print > \"",outpath,"\"$1\".var_overlap.readsout_downsampled_",nreads,"\"}}'")
      }else{
        samtools_cmd2 <- paste0(" | cut -f3,4,5,6,10,12,13,14 | grep '",genetag,"' | sed 's/",genetag,":Z://' | sed 's/UB:Z://' | sed 's/",BCtag,":Z://' | awk 'BEGIN{IFS=\"\t\";OFS=\"\t\";}{print $1,$2,$3,$4,$5,$6,$8,$7;}' | awk '{if($3 == \"255\"){print > \"",outpath,"\"$1\".var_overlap.readsout\"}}'")
      }
      
    }else{
      samtools_cmd1 <- "view -@2 -x QB -x QU -x ES -x IS -x EN -x IN -x GI -x BX -x UX -x NH -x AS -x nM -x HI -x IH -x NM -x uT -x MD -x jM -x jI -x XN -x XS -x vA -x vG -x vW"
      
      # AL EDIT: add downsampling
      if(!is.null(nreads)){
        samtools_cmd2 <- paste0(" | cut -f3,4,5,6,10,12,13,14 | shuf -n ",nreads," | grep '",genetag,"' | sed 's/",genetag,":Z://' | sed 's/UB:Z://' | sed 's/",BCtag,":Z://' | awk '{if($3 == \"255\"){print > \"",outpath,"\"$1\".var_overlap.readsout_downsampled_",nreads,"\"}}'")
      }else{
        samtools_cmd2 <- paste0(" | cut -f3,4,5,6,10,12,13,14 | grep '",genetag,"' | sed 's/",genetag,":Z://' | sed 's/UB:Z://' | sed 's/",BCtag,":Z://' | awk '{if($3 == \"255\"){print > \"",outpath,"\"$1\".var_overlap.readsout\"}}'")
      }
    }
  }else{
    samtools_cmd1 <- "view -@2 -x QB -x QU -x ES -x IS -x EN -x IN -x GI -x BX -x UX -x NH -x AS -x nM -x HI -x IH -x NM -x uT -x MD -x jM -x jI -x XN -x XS -x vA -x vG -x vW -x UB"
    
    # AL EDIT: add downsampling
    if(!is.null(nreads)){
      samtools_cmd2 <- paste0(" | cut -f3,4,5,6,10,12,13 | shuf -n ",nreads," | grep '",genetag,"' | sed 's/",genetag,":Z://' | sed 's/",BCtag,":Z://' | awk '{if($3 == \"255\"){print > \"",outpath,"\"$1\".var_overlap.readsout_downsampled_",nreads,"\"}}'")
    }else{
      samtools_cmd2 <- paste0(" | cut -f3,4,5,6,10,12,13 | grep '",genetag,"' | sed 's/",genetag,":Z://' | sed 's/",BCtag,":Z://' | awk '{if($3 == \"255\"){print > \"",outpath,"\"$1\".var_overlap.readsout\"}}'")
    }
  }
  samtools_cmd <- paste(samtoolsexc,samtools_cmd1,path_bam,samtools_cmd2)
  system(samtools_cmd)
}

print("Done")
Sys.time()


# crunch data ---------------------------------------------------------------

vcf_list <- split(vcf[CHROM %in% chroms_todo], by = "CHROM")

if(UMIdata_flag == "UMI"){
  if(map_flag){
    print("Producing molecule assignment map...")
    map_out_list <- lapply(vcf_list, function(x) calc_coverage_new_return_map(vcf_chunk = x, out = outpath, cellBCs = cellBCs, type = "UMIs", read_layout = opt$read_layout ))
    map_out <- rbindlist(map_out_list)
    fwrite(map_out, file = paste0(outpath,"molecule_assignments.txt" ), sep= "\t", quote = F)
    print("Continuing with allelic expression tables...")
  }

  out_list <- lapply(vcf_list, function(x) calc_coverage_new(vcf_chunk = x, out = outpath, cellBCs = cellBCs, type = c("reads","UMIs"), read_layout = opt$read_layout ))
  read_list <- lapply(out_list, function(x) x$reads)
  UMI_list <- lapply(out_list, function(x) x$UMIs)

  out_reads <- rbindlist(read_list)
  out_UMIs <- rbindlist(UMI_list)
}else{
  out_list <- lapply(vcf_list, function(x) calc_coverage_new(vcf_chunk = x, out = outpath, cellBCs = cellBCs, type = "reads", read_layout = opt$read_layout ))
  out_reads <- rbindlist(out_list)
}


print("Finalizing converting & output ...")
Sys.time()
out_reads <- out_reads[ (cast+c57) >= minC ]

CAST_reads <- makeWide(allele_dat = out_reads, metric = "cast")
BL6_reads <- makeWide(allele_dat = out_reads, metric = "c57")
fract_CAST <- makeWide(allele_dat = out_reads, metric = "CAST_fraction")


print("Processing complete, writing output...")
Sys.time()
# AL EDIT: add subsample
if(!is.null(nreads)){
  fwrite(CAST_reads, file = paste0(outpath,"CAST_reads_downsampled_",nreads,".txt" ), sep= "\t", quote = F)
  fwrite(BL6_reads, file = paste0(outpath,"BL6_reads_downsampled_",nreads,".txt" ), sep= "\t", quote = F)
  fwrite(fract_CAST, file = paste0(outpath,"fract_CAST_reads_downsampled_",nreads,".txt" ), sep= "\t",na = "NA", quote = F)
}else{
  fwrite(CAST_reads, file = paste0(outpath,"CAST_reads.txt" ), sep= "\t", quote = F)
  fwrite(BL6_reads, file = paste0(outpath,"BL6_reads.txt" ), sep= "\t", quote = F)
  fwrite(fract_CAST, file = paste0(outpath,"fract_CAST_reads.txt" ), sep= "\t",na = "NA", quote = F)
}


if(UMIdata_flag == "UMI"){
  out_UMIs <- out_UMIs[ (cast+c57) >= minC ]
  #get directly counted UMIs and write them
  CAST_UMIs <- makeWide(allele_dat = out_UMIs, metric = "cast")
  BL6_UMIs <- makeWide(allele_dat = out_UMIs, metric = "c57")
  fract_CAST_UMIs <- makeWide(allele_dat = out_UMIs, metric = "CAST_fraction")
  
  # AL EDIT: add subsampling
  if(!is.null(nreads)){
    fwrite(CAST_UMIs, file = paste0(outpath,"CAST_direct_UMIs_downsampled_",nreads,".txt" ), sep= "\t", quote = F)
    fwrite(BL6_UMIs, file = paste0(outpath,"BL6_direct_UMIs_downsampled_",nreads,".txt" ), sep= "\t", quote = F)
    fwrite(fract_CAST_UMIs, file = paste0(outpath,"fract_CAST_direct_UMIs_downsampled_",nreads,".txt" ), sep= "\t", na = "NA", quote = F)
  }else{
    fwrite(CAST_UMIs, file = paste0(outpath,"CAST_direct_UMIs.txt" ), sep= "\t", quote = F)
    fwrite(BL6_UMIs, file = paste0(outpath,"BL6_direct_UMIs.txt" ), sep= "\t", quote = F)
    fwrite(fract_CAST_UMIs, file = paste0(outpath,"fract_CAST_direct_UMIs.txt" ), sep= "\t", na = "NA", quote = F)
  }
  
  #also convert total UMI counts into fractional allele counts with read count derived allele fractions
  dge <- paste(opt$out_dir,"/zUMIs_output/expression/",opt$project,".dgecounts.rds",sep="")
  UMIs <- makeUMIs(dge_path = dge, fract_CAST)
  
  # AL EDIT: add subsampling
  if(!is.null(nreads)){
    write.table(UMIs$umis_CAST, file = paste0(outpath,"CAST_fractional_UMIs_downsampled_",nreads,".txt" ), sep= "\t", quote = F)
    write.table(UMIs$umis_BL6, file = paste0(outpath,"BL6_fractional_UMIs_downsampled_",nreads,".txt" ), sep= "\t", quote = F)
  }else{
    write.table(UMIs$umis_CAST, file = paste0(outpath,"CAST_fractional_UMIs.txt" ), sep= "\t", quote = F)
    write.table(UMIs$umis_BL6, file = paste0(outpath,"BL6_fractional_UMIs.txt" ), sep= "\t", quote = F)
  }
}


print("DONE")
Sys.time()
