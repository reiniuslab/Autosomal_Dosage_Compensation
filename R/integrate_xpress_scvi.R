library(scater)
library(scran)

sce.xpress.merge <- readRDS("data/sce_xpress_merge.rds")

hvgs <- getTopHVGs(modelGeneVar(sce.xpress.merge, block=sce.xpress.merge$dataset), n = 500, fdr.threshold = 0.05)

## Integrate data
library(reticulate)

sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)

sce.xpress.merge.red <- sce.xpress.merge[hvgs,]

# import data
adata <- sc$AnnData(
  X   = t(as.matrix(assay(sce.xpress.merge.red,slot='counts'))),
  obs = as.data.frame(colData(sce.xpress.merge.red))
)

# setup and train scVI model
scvi$model$SCVI$setup_anndata(adata, batch_key = 'dataset')
scvi_model <- scvi$model$SCVI(adata, gene_likelihood = 'nb')

scvi_model$train()
#scvi_model$save("scvi_model/", save_anndata = 'True', overwrite = 'True')

## add scVI data to Seurat object
# get the latent represenation
scvi_latent <- as.matrix(scvi_model$get_latent_representation(adata))
rownames(scvi_latent) <- colnames(sce.xpress.merge)
colnames(scvi_latent) <- paste0("scvi_",1:ncol(scvi_latent))
reducedDim(sce.xpress.merge, "scvi") <- scvi_latent

# Rerun clustering on scVI latent
library(magrittr)
sce.xpress.merge %<>% runTSNE(dimred="scvi")
sce.xpress.merge %<>% runUMAP(dimred="scvi")
  
# export
saveRDS(sce.xpress.merge, "data/sce_xpress_scvi.rds")