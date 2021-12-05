## the command below is a one-line shortcut for:
## library(BiocManager)
## install("SingleCellExperiment")
BiocManager::install("SingleCellExperiment")
BiocManager::install(c('scater', 'scran', 'uwot'))
BiocManager::install("AnnotationHub")
library(SingleCellExperiment)

counts_matrix <- data.frame(cell_1 = rpois(10, 10), 
                            cell_2 = rpois(10, 10), 
                            cell_3 = rpois(10, 30))
rownames(counts_matrix) <- paste0("gene_", 1:10)
counts_matrix <- as.matrix(counts_matrix) # must be a matrix object!

sce <- SingleCellExperiment(assays = list(counts = counts_matrix))
sce <- scater::logNormCounts(sce)
logcounts(sce)

cell_metadata <- data.frame(batch = c(1, 1, 2))
rownames(cell_metadata) <- paste0("cell_", 1:3)
sce <- SingleCellExperiment(assays = list(counts = counts_matrix),
                            colData = cell_metadata)
sce$more_stuff
sce <- scater::addPerCellQC(sce)
colData(sce)

sce$more_stuff <- runif(ncol(sce))
colnames(colData(sce))

sce <- scater::addPerFeatureQC(sce)
rowData(sce)

library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]] # Human, Ensembl v97.
genes(edb)[,2]














































































































