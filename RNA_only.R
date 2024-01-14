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
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
#library(monocle3)
library(clusterProfiler)
library(org.Mm.eg.db)
library(Seurat)
library(SeuratDisk)
library(Matrix)
set.seed(1234)

setwd("~/Desktop/有事/Harvard/Research/COPD/Data")

################## RNA files ################## 
# Load in data, already processed h5ad file
RNA_data <- Read10X(data.dir = './220511_filtered_RNA/')
metadata_RNA <- read.csv('./220511_filtered/metadata.csv')
rownames(metadata_RNA) <- metadata_RNA$X

# Create Seurat object
Rawlin_220511 <- CreateSeuratObject(counts = RNA_data,
                                    meta.data = metadata_RNA)
saveRDS(Rawlin_220511, file = "./RData/Rawlin_220511.rds")
Rawlin_220511 <- readRDS("./RData/Rawlin_220511.rds")

## Data Overview
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Rawlin_220511[["percent.mt"]] <- PercentageFeatureSet(Rawlin_220511, pattern = "^MT-")
# Visualize QC metrics as a violin plot
which(is.na(Rawlin_220511$percent.mt)) # 0
which(is.na(Rawlin_220511$nFeature_RNA)) # 0
which(is.na(Rawlin_220511$nCount_RNA)) # 0
# Remove NA "X"
# Rawlin_220511 <- Rawlin_220511[,!colnames(Rawlin_220511) %in% "X"]
# length(colnames(Rawlin_220511)) ## 33539 -> 33538

# zero_nCount_RNA <- sum(Rawlin_220511$nCount_RNA <= 0)
# all_nCount_RNA <- sum(Rawlin_220511$nCount_RNA)
# To print the result
print(zero_nCount_RNA)
print(zero_nCount_RNA/all_nCount_RNA)

small <- sum(Rawlin_220511$nFeature_RNA <= 1000)
cnt <- sum(Rawlin_220511$nFeature_RNA)
# To print the result
print(small)
print(small/cnt)

# saveRDS(Rawlin_220511, file = "./Data/Rawlin_220511_preFilter.rds")
# readRDS("./Data/Rawlin_220511_preFilter.rds")

VlnPlot(Rawlin_220511, features = c("nFeature_RNA","nCount_RNA"),pt.size = 0)
FeatureScatter(Rawlin_220511, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Rawlin_220511 <- subset(Rawlin_220511, subset = nFeature_RNA > 100 &
                          nFeature_RNA < 9000 & nCount_RNA > 1000)
saveRDS(Rawlin_220511, file = "./RData/Rawlin_220511_filter1.rds")
Rawlin_220511 <- readRDS("./RData/Rawlin_220511_filter1.rds")

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
