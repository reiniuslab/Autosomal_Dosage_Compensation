library(data.table)
set.seed("999")

## functions

mad_filter <- function(x, nmads = 3) {
  med <- median(x, na.rm = TRUE)
  m <- mad(x, constant = 1, na.rm = TRUE)
  x > (med - nmads * m) & x < (med + nmads * m)
}


download_data <- function(url, dest) {
  if (!file.exists(dest)) {
    message("Downloading: ", dest)
    download.file(url, destfile = dest, mode = "wb")
  } else {
    message("Already exists, skipping: ", dest)
  }
}


make_ratio_plot <- function(min_frac = 0, plot_title = "No filter") {
  require(ggplot2)
  
  keep_genes <- names(gene_detection_fraction)[gene_detection_fraction >= min_frac]
  
  combined_sub <- subset(combined, features = keep_genes)
  
  pb_mat <- AggregateExpression(combined_sub,
                                group.by = "condition",
                                assays = "RNA",
                                slot = "counts",
                                return.seurat = FALSE)$RNA
  
  cpm <- sweep(pb_mat, 2, colSums(pb_mat), "/") * 1e6
  
  ratio_dt <- data.table(gene_id = rownames(cpm),
                         diploid = cpm[, "diploid"],
                         trisomy21 = cpm[, "trisomy21"])
  
  ratio_dt[gtf, "chr" := seqnames, on = "gene_id == gene_id"]
  ratio_dt <- ratio_dt[!is.na(chr)]
  ratio_dt <- ratio_dt[!(chr %in% c("chrM", "chrX", "chrY"))]
  ratio_dt[, "chr" := factor(chr, levels = paste0("chr", 1:22))]
  ratio_dt <- ratio_dt[diploid > 0 & trisomy21 > 0]
  
  ratio_dt[, "ratio" := trisomy21 / diploid]
  ratio_dt[, "chr_group" := ifelse(chr == "chr21", "chr21", "autosome")]
  
  auto_med <- median(ratio_dt[chr_group == "autosome"]$ratio, na.rm = TRUE)
  ratio_dt[, "ratio_scaled" := ratio / auto_med]
  
  expected_chr21_scaled <- 1.5 / auto_med
  
  label_dt <- ratio_dt[, .(n_genes = .N), by = chr]
  label_dt[, "label" := paste0("n=", n_genes)]
  label_dt[, "y" := 0] # bottom of plot position
  
  chr21_ratios <- ratio_dt[chr == "chr21"]$ratio_scaled
  
  wt <- wilcox.test(chr21_ratios, mu = expected_chr21_scaled, alternative = "two.sided")
  
  stats_dt <- data.table(filter = plot_title,
                         min_fraction = min_frac,
                         n_genes_total = nrow(ratio_dt),
                         n_chr21_genes = nrow(ratio_dt[chr == "chr21"]),
                         auto_median_raw_ratio = auto_med,
                         expected_chr21_scaled = expected_chr21_scaled,
                         chr21_median_scaled = median(chr21_ratios, na.rm = TRUE),
                         wilcox_p = wt$p.value)
  
  p <- ggplot(ratio_dt, aes(x = chr, y = ratio_scaled, fill = chr_group)) +
    geom_boxplot(width = 0.6, outlier.alpha = 0.08, outlier.color = "grey60") +
    geom_text(data = label_dt,
              aes(x = chr, y = y, label = label),
              inherit.aes = FALSE,
              size = 2.2,
              angle = 90,
              vjust = 0.5,
              hjust = 0) +
    coord_cartesian(ylim = c(0, 2)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_segment(x = 20.5,
                 xend = 21.5,
                 y = expected_chr21_scaled,
                 yend = expected_chr21_scaled,
                 linewidth = 1,
                 col = "red") +
    scale_fill_manual(values = c("autosome" = "white", "chr21" = "#4DBBD5")) +
    labs(title = plot_title,
         x = "Chromosome",
         y = "Trisomy21 / diploid CPM ratio",
         fill = "Chromosome") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  list(plot = p, stats = stats_dt, data = ratio_dt)
}


## directories

data_dir <- "data"
pmid_dir <- file.path(data_dir, "pmid_31582743")

dir.create(pmid_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("rev_plots", recursive = TRUE, showWarnings = FALSE)

## download data

data_urls <- c(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4013nnn/GSM4013118/suppl/GSM4013118%5Ftwin1B%2E10x%2Eprocessed%2Ecsv%2Egz",
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4013nnn/GSM4013119/suppl/GSM4013119%5Ftwin2A%2E10x%2Eprocessed%2Ecsv%2Egz"
)

gencode_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh37_mapping/gencode.v37lift37.annotation.gtf.gz"

data_dest <- file.path(pmid_dir, basename(URLdecode(data_urls)))
gencode_dest <- file.path(data_dir, basename(gencode_url))

for (i in seq_along(data_urls)) {
  download_data(data_urls[i], data_dest[i])
}

download_data(gencode_url, gencode_dest)

rm(data_urls, gencode_url, data_dest, gencode_dest, i)

## load data
# load annotation
library(rtracklayer)
gtf <- as.data.table(rtracklayer::import(file.path(data_dir, "gencode.v37lift37.annotation.gtf.gz")))
gtf <- gtf[type == "gene"]
gtf[, gene_id := tstrsplit(gene_id, "\\.", fixed = FALSE, keep = 1)]

# load count matrices
sc_trisomy21 <- fread(file.path(pmid_dir, "GSM4013118_twin1B.10x.processed.csv.gz"))
sc_diploid   <- fread(file.path(pmid_dir, "GSM4013119_twin2A.10x.processed.csv.gz"))

## process data
# matrices & gene names as rownames
gene_ids_tri <- sc_trisomy21$V1
gene_ids_dip <- sc_diploid$V1

mat_trisomy <- as.matrix(sc_trisomy21[, -1, with = FALSE])
mat_diploid <- as.matrix(sc_diploid[, -1, with = FALSE])

rownames(mat_trisomy) <- gene_ids_tri
rownames(mat_diploid) <- gene_ids_dip

# create Seurat objects
library(Seurat)
library(sctransform)
diploid_obj <- CreateSeuratObject(mat_diploid)
diploid_obj$condition <- "diploid"

trisomy_obj <- CreateSeuratObject(mat_trisomy)
trisomy_obj$condition <- "trisomy21"

# remove outlier cells
keep_diploid <-
  mad_filter(diploid_obj$nFeature_RNA, 3) &
  mad_filter(diploid_obj$nCount_RNA, 3)

keep_trisomy <-
  mad_filter(trisomy_obj$nFeature_RNA, 3) &
  mad_filter(trisomy_obj$nCount_RNA, 3)

cat("Diploid removed:", sum(!keep_diploid), "of", length(keep_diploid), "\n") # Diploid removed: 404 of 3801 
cat("Trisomy removed:", sum(!keep_trisomy), "of", length(keep_trisomy), "\n") # Trisomy removed: 523 of 4939 

diploid_obj <- subset(diploid_obj, cells = colnames(diploid_obj)[keep_diploid])
trisomy_obj <- subset(trisomy_obj, cells = colnames(trisomy_obj)[keep_trisomy])

# combine datasets
combined <- merge(diploid_obj,
                  y = trisomy_obj,
                  add.cell.ids = c("diploid", "trisomy21"))

# SCTransform normalization
combined <- SCTransform(combined, 
                        assay = "RNA", 
                        new.assay.name = "SCT", 
                        vst.flavor = "v2", 
                        vars.to.regress = NULL, 
                        verbose = FALSE)


# SCT-based single-cell
DefaultAssay(combined) <- "SCT"
combined <- RunPCA(combined, assay = "SCT", verbose = FALSE)
combined <- RunUMAP(combined, dims = 1:30, reduction = "pca", verbose = FALSE)
combined <- FindNeighbors(combined, dims = 1:30, reduction = "pca", verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.5, verbose = FALSE)

# pseudobulk
DefaultAssay(combined) <- "RNA"
combined[["RNA"]] <- JoinLayers(combined[["RNA"]])
counts_mat <- LayerData(combined, assay = "RNA", layer = "counts")
total_cells <- ncol(counts_mat)
gene_detection_fraction <- Matrix::rowSums(counts_mat > 0) / total_cells

# run the different filters
res_no_filter <- make_ratio_plot(0, "No filter")
res_10 <- make_ratio_plot(0.10, "Genes detected in >= 10% of cells")
res_15 <- make_ratio_plot(0.15, "Genes detected in >= 15% of cells")
res_20 <- make_ratio_plot(0.20, "Genes detected in >= 20% of cells")

## plotting
library(ggplot2)
library(cowplot)
library(patchwork)

ratio_plot_combined <-
  res_no_filter$plot +
  res_10$plot +
  res_15$plot +
  res_20$plot +
  plot_layout(ncol = 4, guides = "collect")

ggsave2(
  filename = "plots/pmid_31582743_aratio_detection_filters.pdf",
  plot = ratio_plot_combined,
  width = 22,
  height = 5
)

## statistics
stats_table <- rbindlist(list(res_no_filter$stats, res_10$stats, res_15$stats, res_20$stats))
stats_table
fwrite(stats_table, file = "plots/pmid_31582743_aratio_detection_filter_stats.tsv", sep = "\t")
