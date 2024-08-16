if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnsDb.Hsapiens.v86")
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(BiocGenerics)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(dplyr)
library(IRanges)
# library(chromVAR)
# library(JASPAR2020)
# library(TFBSTools)
# library(motifmatchr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(irlba)
install.packages("Matrix", type = "source")
install.packages("irlba", type = "source")
set.seed(1234)

setwd("~/Desktop/有事/Harvard/Research/COPD/Data/ATAC/655_cellranger-atac200_count_5891STDY8038655_GRCh38-2020-A-2_0_0")
################## ATAC files ################## 
# part1 <- Read10X(data.dir = './part1/')
# part2 <- Read10X(data.dir = './part2/', gene.column = 1)
# test <- Read10X(data.dir = ".")
# ATAC_data <- Read10X(data.dir = './part1_new.mtx')
# metadata_atac <- read.csv('./ATAC/metadata.csv')
# rownames(metadata_atac) <- metadata_atac$X

counts <- Read10X_h5(filename = "./filtered_peak_bc_matrix.h5")
metadata <- read.csv(file = "./singlecell.csv", header = TRUE, row.names = 1)
chrom_assay <- CreateChromatinAssay(counts = counts, sep = c(":", "-"),
                                    fragments = './fragments.tsv.gz', 
                                    min.cells = 10, min.features = 200)
atac_655 <- CreateSeuratObject(counts = chrom_assay, assay = "peaks",
                               meta.data = metadata)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
Annotation(atac_655) <- annotations

atac_655
################## END IMPORT  ################## 
# compute nucleosome signal score per cell
atac_655 <- NucleosomeSignal(atac_655)
# compute TSS enrichment score per cell
atac_655 <- TSSEnrichment(atac_655)

saveRDS(atac_655, "../RData/atac_655_1.rds")
atac_655 <- readRDS("../RData/atac_655_1.rds")
# add blacklist ratio and fraction of reads in peaks
atac_655$pct_reads_in_peaks <- atac_655$peak_region_fragments / atac_655$passed_filters * 100
atac_655$blacklist_ratio <- atac_655$blacklist_region_fragments / atac_655$peak_region_fragments

atac_655$high.tss <- ifelse(atac_655$TSS.enrichment > 3, 'High', 'Low')
table(atac_655$high.tss)

VlnPlot(atac_655, features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 
                               'nucleosome_signal', 'pct_reads_in_peaks'),
        pt.size = 0.1, ncol = 5)

atac_655 <- subset(x = atac_655,
                   subset = nCount_peaks > 3000 &
                     nCount_peaks < 30000 &
                     pct_reads_in_peaks > 15 &
                     blacklist_ratio < 0.05 &
                     nucleosome_signal < 4 &
                     TSS.enrichment > 3)
atac_655

# Normalization and linear dimensional reduction
atac_655 <- RunTFIDF(atac_655)
atac_655 <- FindTopFeatures(atac_655, min.cutoff = 'q0')
atac_655 <- RunSVD(atac_655)

atac_655 <- RunUMAP(object = atac_655, reduction = 'lsi', dims = 2:30)
atac_655 <- FindNeighbors(object = atac_655, reduction = 'lsi', dims = 2:30)
atac_655 <- FindClusters(object = atac_655, verbose = FALSE, algorithm = 3)
DimPlot(object = atac_655, label = TRUE) + NoLegend()

saveRDS(atac_655, "../RData/atac_655_2.rds")
atac_655 <- readRDS("../RData/atac_655_2.rds")
##
Rawlin_220511 <- readRDS("~/Desktop/有事/Harvard/Research/COPD/Data/RData/Rawlin_220511_normalize2.rds")


multiomics <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

######  Normalizing the data ###### 
Rawlin_220511 <- NormalizeData(Rawlin_220511, normalization.method = "LogNormalize", 
                               scale.factor = 10000)
Rawlin_220511 <- FindVariableFeatures(Rawlin_220511, selection.method = "vst", 
                                      nfeatures = 2000)
Rawlin_220511 <- ScaleData(Rawlin_220511)
Rawlin_220511 <- RunPCA(Rawlin_220511, features = VariableFeatures(object=Rawlin_220511))
saveRDS(Rawlin_220511, file = "./RData/Rawlin_220511_normalize2.rds")
Rawlin_220511 <- readRDS("./RData/Rawlin_220511_normalize2.rds")

###### Cluster ######  
Rawlin_220511 <- FindNeighbors(Rawlin_220511, dims = 1:10)
Rawlin_220511 <- FindClusters(Rawlin_220511, resolution = 0.3)
# saveRDS(Rawlin_220511, file = "./RData/Rawlin_220511_cluster5.rds")
saveRDS(Rawlin_220511, file = "./RData/Rawlin_220511_cluster5+.rds")
Rawlin_220511 <- readRDS("./RData/Rawlin_220511_cluster5.rds")

# Look at cluster IDs of the first 5 cells
head(Idents(Rawlin_220511), 5)

# Run non-linear dimensional reduction (UMAP/tSNE)
Rawlin_220511 <- RunUMAP(Rawlin_220511, dims = 1:10)

pdf("./Plot/cluster1.pdf")
DimPlot(Rawlin_220511, reduction = "umap")
dev.off()

######  Finding differentially expressed features (cluster biomarkers) ###### 
# find all markers of cluster 2
cluster2.markers <- FindMarkers(Rawlin_220511, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(Rawlin_220511, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
Rawlin_220511.markers <- FindAllMarkers(Rawlin_220511, only.pos = TRUE)
Rawlin_220511.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster0.markers <- FindMarkers(Rawlin_220511, ident.1 = 0, logfc.threshold = 0.25, 
                                test.use = "roc", only.pos = TRUE)
VlnPlot(Rawlin_220511, features = c("UBE2C", "CENPE"))

FeaturePlot(Rawlin_220511, features = c("UBE2C", "CENPE"))

# ###### Cell Cycle Scoring and Regression ###### 
# exp.mat <- read.table(file = "./cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt",
#                       header = TRUE, as.is = TRUE, row.names = 1)
# 
# # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# # segregate this list into markers of G2/M phase and markers of S phase
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# 
# Rawlin_220511 <- CellCycleScoring(Rawlin_220511, s.features = s.genes, 
#                                   g2m.features = g2m.genes, set.ident = TRUE)
# saveRDS(Rawlin_220511, file = "./RData/Rawlin_220511_cellcycle3.rds")
# Rawlin_220511 <- readRDS("./RData/Rawlin_220511_cellcycle3.rds")
# # view cell cycle scores and phase assignments
# head(Rawlin_220511[[]])
# # RidgePlot(Rawlin_220511, features = c("CD74", "COL1A1", "IGLV5-37", "COL1A2"), ncol = 2)
# 
# # Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by phase
# Rawlin_220511 <- RunPCA(Rawlin_220511, features = c(s.genes, g2m.genes))
# DimPlot(Rawlin_220511)
# 
# # Regress out cell cycle scores during data scaling
# Rawlin_220511 <- ScaleData(Rawlin_220511, vars.to.regress = c("S.score", "G2M.score"), 
#                            features = rownames(Rawlin_220511))
# # saveRDS(Rawlin_220511, file = "./RData/Rawlin_220511_cellcycle4.rds")
# # Rawlin_220511 <- readRDS("./RData/Rawlin_220511_cellcycle4.rds")
# 
# # Warning: Error: vector memory exhausted (limit reached?)?
# # Need cloud service
# 
# # Now, a PCA on the variable genes no longer returns components associated with cell cycle
# Rawlin_220511 <- RunPCA(Rawlin_220511, features = VariableFeatures(Rawlin_220511), 
#                         nfeatures.print = 10)
# # When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
# Rawlin_220511 <- RunPCA(Rawlin_220511, features = c(s.genes, g2m.genes))
# DimPlot(Rawlin_220511)
# 
# # Examine and visualize PCA results a few different ways
# print(Rawlin_220511[["pca"]], dims = 1:5, nfeatures = 5)
# VizDimLoadings(Rawlin_220511, dims = 1:2, nfeatures = 10, reduction = "pca")
# DimPlot(Rawlin_220511, reduction = "pca") + NoLegend()
# 
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(Rawlin_220511), 10)
# 
# # plot variable features with and without labels
# # plot1 <- VariableFeaturePlot(Rawlin_220511)
# # plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# # plot2
# VlnPlot(Rawlin_220511, c("UBE2C", "CENPE", "AGBL4"))







