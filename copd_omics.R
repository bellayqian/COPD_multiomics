library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(BiocGenerics)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(dplyr)
library(biovizBase)
library(clusterProfiler)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
# library(IRanges)
#library(monocle3)
set.seed(1234)

setwd("/udd/reyqi/COPD")

# load the RNA and ATAC data, all in./atac_fragments.tsv.gz
# counts <- Read10X_h5("/proj/rerefs/reref00/rna/whole/singlecell_lung/Rawlins/filtered_tf_bc_matrix.h5")
counts <- Read10X_h5("/proj/rerefs/reref00/rna/whole/singlecell_lung/Rawlins/data/raw/cellranger-atac200_count_5891STDY8038655_GRCh38-2020-A-2_0_0/filtered_peak_bc_matrix")
fragpath <- "./RData/fragments.tsv.gz"

# create a Seurat object containing the RNA data
multi.Rawlins <- CreateSeuratObject(counts = counts, assay = "RNA")

# get gene annotations for hg38
# annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotation) <- "UCSC"
# genome(annotation) <- "hg38"
# saveRDS(annotation, "./RData/annotation.rds")
annotation <- readRDS("./RData/annotation.rds")

# create ATAC assay and add it to the object
# counts_atac <- Read10X_h5("/proj/rerefs/reref00/rna/whole/singlecell_lung/Rawlins/filtered_peak_bc_matrix.h5")
counts_atac <- Read10X_h5("/proj/rerefs/reref00/rna/whole/singlecell_lung/Rawlins/data/raw/cellranger-atac200_count_5891STDY8038655_GRCh38-2020-A-2_0_0/filtered_tf_bc_matrix")

multi.Rawlins[["ATAC"]] <- CreateChromatinAssay(counts = counts_atac, sep = c(":", "-"),
                                                fragments = fragpath, annotation = annotation)
multi.Rawlins
saveRDS(multi.Rawlins, "./RData/multi.Rawlins_orig.rds")
multi.Rawlins <- readRDS("./RData/multi.Rawlins_orig.rds")

## Quality Control
# change defaultAssay from RNA to ATAC
DefaultAssay(multi.Rawlins) <- "ATAC"

# Quality control
multi.Rawlins <- NucleosomeSignal(multi.Rawlins)
multi.Rawlins <- TSSEnrichment(multi.Rawlins)
VlnPlot(
  object = multi.Rawlins,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
saveRDS(multi.Rawlins, "./RData/multi.Rawlins_2.rds")
multi.Rawlins <- readRDS("./RData/multi.Rawlins_2.rds")
# # filter out low quality cells
# multi.Rawlins <- subset(
#   x = multi.Rawlins,
#   subset = nCount_RNA < 750000 &
#     nCount_RNA > 1000 &
#     nCount_ATAC > 1000 &
#     nCount_ATAC < 150000 &
#     nucleosome_signal < 5 &
#     TSS.enrichment > 1 & 
#     TSS.enrichment < 12
# )

# pdf(file = "./Plots/Vlnplot_QC.pdf", width = 15, height = 15)
# VlnPlot(
#   object = multi.Rawlins,
#   features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
#   ncol = 4,
#   pt.size = 0
# )
# dev.off()

# QC
saveRDS(multi.Rawlins, "./RData/multi.Rawlins_3.rds")
multi.Rawlins <- readRDS("./RData/multi.Rawlins_3.rds")

# ## Gene expression data processing
# ## Cell-Cycle Scoring and Regression
# DefaultAssay(multi.Rawlins) <- "RNA"
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes

## Joint UMAP visualization
# build a joint neighbor graph using both assays
# RNA analysis
DefaultAssay(multi.Rawlins) <- "RNA"
multi.Rawlins <- NormalizeData(multi.Rawlins)
multi.Rawlins <- FindVariableFeatures(multi.Rawlins, selection.method = "vst")
multi.Rawlins <- ScaleData(multi.Rawlins, features = rownames(multi.Rawlins))
multi.Rawlins <- RunPCA(multi.Rawlins, features = VariableFeatures(multi.Rawlins))
multi.Rawlins <- RunUMAP(multi.Rawlins, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
# DimHeatmap(multi.Rawlins, dims = c(8, 10))
saveRDS(multi.Rawlins, "./RData/multi.Rawlins_4.rds")
multi.Rawlins <- readRDS("./RData/multi.Rawlins_4.rds")
# Download multi.Rawlins_4.rds from cloud service

# multi.Rawlins <- CellCycleScoring(multi.Rawlins, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# # Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by phase
# multi.Rawlins <- RunPCA(multi.Rawlins, features = c(s.genes, g2m.genes))
# p1 <- DimPlot(multi.Rawlins)
# multi.Rawlins <- ScaleData(multi.Rawlins, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(multi.Rawlins))
# saveRDS(multi.Rawlins, "./IntermediateFile/multi.Rawlins_CellCycle_RregressS_G2M.rds")
# multi.Rawlins <- readRDS("./IntermediateFile/multi.Rawlins_CellCycle_RregressS_G2M.rds")
# multi.Rawlins <- RunPCA(multi.Rawlins)
# p2 <- DimPlot(multi.Rawlins)
# p1+p2

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
# install.packages("Matrix", type = "source")
# install.packages("irlba", type = "source")
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(multi.Rawlins) <- "ATAC"
multi.Rawlins <- RunTFIDF(multi.Rawlins)
multi.Rawlins <- FindTopFeatures(multi.Rawlins, min.cutoff = 'q0')
multi.Rawlins <- RunSVD(multi.Rawlins)
multi.Rawlins <- RunUMAP(multi.Rawlins, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
multi.Rawlins <- FindMultiModalNeighbors(object = multi.Rawlins,reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
multi.Rawlins <- RunUMAP(multi.Rawlins, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
multi.Rawlins <- FindClusters(multi.Rawlins, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
saveRDS(multi.Rawlins, "./RData/multi.Rawlins_5.rds")
# obtained in local R
multi.Rawlins <- readRDS("./RData/multi.Rawlins_5.rds")

pdf(file = "./Plots/Dimplot_cluster.pdf", width = 12, height = 6)
p1 <- DimPlot(multi.Rawlins, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(multi.Rawlins, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(multi.Rawlins, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(file = "./Plots/WNN_cluster.pdf", width = 6, height = 6)
p3 <- DimPlot(multi.Rawlins, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

## Feature Plots

## Find Marker Genes for multi.Rawlins dataset
DefaultAssay(multi.Rawlins) <- "RNA"
all.markers <- FindAllMarkers(object = multi.Rawlins)
saveRDS(all.markers, "./RData/all.markers.rds")
all.markers <- readRDS("./RData/all.markers.rds")

for (i in 0:19){
  name <- paste("clusterMarker.", i, sep = "")
  markers <- assign(name, dplyr::filter(all.markers, all.markers$cluster ==i))
  write.table(markers, file  = paste0("./Results/Cluster",i,"_Marker.txt"), sep = "\t")
}
## Pathway and Gene Ontology (GO) Enrichment Analysis
# Input to GO_Enrichment_Analysis.R to perform GO analysis

# Feb 25, 2024 END
################################################################################
# Apr 9, 2024 Begin
# add annotations
# Use Slingshot.R to build up trajectory DONE
# leave cluster annotation for future reference
# multi.Rawlins <- RenameIdents(multi.Rawlins, '19' = 'pDC','20' = 'HSPC','15' = 'cDC')
# multi.Rawlins <- RenameIdents(multi.Rawlins, '0' = 'CD14 Mono', '9' ='CD14 Mono', '5' = 'CD16 Mono')
# multi.Rawlins$celltype <- Idents(multi.Rawlins)

# Use stackViolin.R to plot marker gene in each cluster, and feature plot DONE

## Linking peaks to genes
DefaultAssay(multi.Rawlins) <- "RNA"
multi.Rawlins <- SCTransform(multi.Rawlins)

DefaultAssay(multi.Rawlins) <- "ATAC"
multi.Rawlins <- RegionStats(multi.Rawlins, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
multi.Rawlins <- LinkPeaks(
  object = multi.Rawlins,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  genes.use = c("SOX9", "ASCL1"))

# from local cuz no BSgenome package
multi.Rawlins <- readRDS("./RData/multi.Rawlins_6.rds")

# idents.plot <- c("2", "18", "15", "11", "14", "10", "0", "13")
# idents.plot <- c(2, 18, 15, 11, 14, 10, 0, 13)

p1 <- CoveragePlot(object = multi.Rawlins,
  region = "SOX9",
  features = "SOX9",
  expression.assay = "SCT",
  extend.upstream = 500,
  extend.downstream = 10000)
p2 <- CoveragePlot(object = multi.Rawlins,
  region = "ASCL1",
  features = "ASCL1",
  expression.assay = "SCT",
  extend.upstream = 8000,
  extend.downstream = 5000)

pdf(file = "./Plots/peakgeneLink.pdf", width = 15, height = 10)
patchwork::wrap_plots(p1, p2, ncol = 1)
dev.off()

saveRDS(multi.Rawlins, "../ATAC/multi.Rawlins_7.rds")
multi.Rawlins <- readRDS("../ATAC/multi.Rawlins_7.rds")

# Apr 15, 2024 END







######################### Backup Codes ######################### 
## Create Peaks assay
metadata <- read.csv(file = "./Inputs/per_barcode_metrics.csv", header = TRUE, row.names = 1)
annotation <- readRDS("./IntermediateFile/annotation.rds")

# call peaks using MACS2
## Feb 24. Unable to install MACS2 or MACS3
# Error
peaks <- CallPeaks(multi.Rawlins, macs2.path = "./env/bin/macs3") 
# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)
saveRDS(peaks, "./RData/peaks.rds")
peaks <- readRDS("./RData/peaks.rds")

# quantify counts in each peak
macs3_counts <- FeatureMatrix(fragments = Fragments(multi.Rawlins), features = peaks, cells = colnames(multi.Rawlins))
saveRDS(macs3_counts, "./RData/macs3_counts.rds")
macs3_counts <- readRDS("./RData/macs3_counts.rds")

# create a new assay using the MACS2 peak set and add it to the Seurat object
multi.Rawlins[["peaks"]] <- CreateChromatinAssay(counts = macs3_counts,
  fragments = './RData/fragments.tsv.gz',annotation = annotation)

multi.Rawlins

DefaultAssay(multi.Rawlins) <- "peaks"

# Add metadata info to multi.Rawlins peaks assay
multi.Rawlins <- AddMetaData(multi.Rawlins, metadata)
saveRDS(multi.Rawlins, "./RData/multi.Rawlins_6.rds")
multi.Rawlins <- readRDS("./RData/multi.Rawlins_6.rds")

# Make Heatmap for each cluster
all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10

pdf(file = "./Plots/Heatmap_multi.RawlinsClusters.pdf", width = 12, height = 12)
DoHeatmap(multi.Rawlins, features = top10$gene)
dev.off()

# multi.Rawlins.subset <- readRDS("./IntermediateFile/multi.Rawlins_subset.rds")
# 
# levels(multi.Rawlins.subset@active.ident) <- c("Fib_1","iFib_1","iCM_1","Fib_2","iFib_2","iFib_3",
#                                             "iFib_4","iFib_5","iFib_6","Fib_3","iCM_2")

##############################
# # If we want to subset
# multi.Rawlins.subset <- subset(x = multi.Rawlins, idents = c("0","1","2","3","4","5","6","7","8","11","12"))
# saveRDS(multi.Rawlins.subset, "./IntermediateFile/multi.Rawlins_subset.rds")
# multi.Rawlins.subset <- readRDS("./IntermediateFile/multi.Rawlins_subset.rds")
# 
# levels(multi.Rawlins.subset@active.ident) <- c("Fib_1","iFib_1","iCM_1","Fib_2","iFib_2","iFib_3",
#                                             "iFib_4","iFib_5","iFib_6","Fib_3","iCM_2")
# 
# pdf(file = "./Plots/Dimplot_subsetCluster.pdf", width = 12, height = 12)
# p1 <- DimPlot(multi.Rawlins.subset, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
# p2 <- DimPlot(multi.Rawlins.subset, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
# p3 <- DimPlot(multi.Rawlins.subset, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
# p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
# dev.off()
# 
# pdf(file = "./Plots/Dimplot_subsetWNN.pdf", width = 5, height = 5)
# DimPlot(multi.Rawlins.subset, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
# dev.off()
# 
# pdf(file = "./Plots/Dimplot_subsetrna.pdf", width = 5, height = 5)
# DimPlot(multi.Rawlins.subset, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) 
# dev.off()
# 
# pdf(file = "./Plots/Dimplot_subsetatac.pdf", width = 5, height = 5)
# DimPlot(multi.Rawlins.subset, reduction = "umap.atac", label.size = 2.5, repel = TRUE)
# dev.off()
##############################

## Compared between cluster 7 and 12 to determine cluster identity
# DefaultAssay(multi.Rawlins.subset) <- "RNA"
# cluster7.12 <- FindMarkers(multi.Rawlins.subset, ident.1 = "7", ident.2 = "12")
# cluster7.12 <- filter(cluster7.12, p_val_adj < 0.05)
# cluster7 <- filter(cluster7.12, avg_log2FC > 0)
# cluster12 <- filter(cluster7.12, avg_log2FC < 0)
# Input to GO_Enrichment_Analysis.R

## Find Marker Genes for multi.Rawlins dataset
DefaultAssay(multi.Rawlins) <- "RNA"
all.markers <- FindAllMarkers(object = multi.Rawlins)
saveRDS(all.markers, ".//all.markers.rds")
all.markers <- readRDS("./IntermediateFile/all.markers.rds")

for (i in 0:12){
  name <- paste("clusterMarker.", i, sep = "")
  markers <- assign(name, dplyr::filter(all.markers, all.markers$cluster ==i))
  write.table(markers, file  = paste0("./Results/Cluster",i,"_Marker.txt"), sep = "\t")
}
# Input to GO_Enrichment_Analysis.R to perform GO analysis

## Find Marker Genes for multi.Rawlins.subset dataset
DefaultAssay(multi.Rawlins.subset) <- "RNA"
clusters <- FindAllMarkers(multi.Rawlins.subset)
# clusters <- filter(clusters, p_val_adj < 0.05)
saveRDS(clusters, "./IntermediateFile/clustersMarkers.rds")
clusters <- readRDS("./IntermediateFile/clustersMarkers.rds")

for (i in 0:12){
  name <- paste("clusterMarker.", i, sep = "")
  markers <- assign(name, dplyr::filter(clusters, clusters$cluster ==i))
  write.table(markers, file  = paste0("./Results/Final GO/Cluster",i,"_Marker.txt"), sep = "\t")
}
# Input to GO_Enrichment_Analysis.R to perform GO analysis


# ###### Does not work ###########
# ## Find differentially accessible peaks between clusters
# DefaultAssay(multi.Rawlins.subset) <- 'peaks'
# da_peaks <- FindMarkers(object = multi.Rawlins.subset, ident.1 = c("12","2"), ident.2 = c("3"),
#                         min.pct = 0.05, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
# 
# head(da_peaks)
# 
# plot1 <- VlnPlot(
#   object = multi.Rawlins.subset,
#   features = rownames(da_peaks)[1],
#   pt.size = 0.1,
#   idents = c("0","3","11","2","12"))
# 
# plot2 <- FeaturePlot(
#   object = multi.Rawlins.subset,
#   features = rownames(da_peaks)[2],
#   pt.size = 0.1)
# 
# plot1 | plot2
# 
# open_iCM <- rownames(da_peaks[da_peaks$avg_log2FC > 0.5, ])
# open_Fibs <- rownames(da_peaks[da_peaks$avg_log2FC < -0.5, ])
# 
# # Empty list for both
# closest_genes_iCM <- ClosestFeature(multi.Rawlins.subset, regions = open_iCM)
# closest_genes_Fibs <- ClosestFeature(multi.Rawlins.subset, regions = open_Fibs)
# 
# head(closest_genes_iCM)
# head(closest_genes_Fibs)
# ###### Probably does not work ###########

## Find Differentially expressed gene list
DefaultAssay(multi.Rawlins.subset) <- 'RNA'
da_motif <- FindMarkers(multi.Rawlins.subset,ident.1 = c("2","12"),ident.2 = c("3","11"))

# saveRDS(da_motif, "./IntermediateFile/enriched_Motifs.rds")
# da_motif <- readRDS("./IntermediateFile/enriched_Motifs.rds")

saveRDS(da_motif, "./IntermediateFile/enriched_Motifs.rds")
da_motif <- readRDS("./IntermediateFile/enriched_Motifs.rds")

da_motif <- filter(da_motif, p_val_adj < 0.05)
clusteriCM <- filter(da_motif, avg_log2FC > 0)
clusterFibs <- filter(da_motif, avg_log2FC < 0)

## Find significant differentially accessible gene-peak linkage
## Linking peaks to genes
DefaultAssay(multi.Rawlins.subset) <- "peaks"

# first compute the GC content for each peak
multi.Rawlins.subset <- RegionStats(multi.Rawlins.subset, genome = BSgenome.Mmusculus.UCSC.mm10)

m.Fibs <- LinkPeaks(object = multi.Rawlins.subset, peak.assay = "peaks",
                    expression.assay = "RNA",
                    genes.use = rownames(clusterFibs))

m.iCM <- LinkPeaks(object = multi.Rawlins.subset, peak.assay = "peaks",
                   expression.assay = "RNA",
                   genes.use = rownames(clusteriCM))

df_Fibs <- as_data_frame(m.Fibs[["peaks"]]@links)
df_iCM <- as_data_frame(m.iCM[["peaks"]]@links)

# write.table(df_Fibs, "./Results/Fibs_PeakLinkage.txt", sep="\t", col.names = NA, row.names = TRUE)
# write.table(df_iCM, "./Results/iCM_PeakLinkage.txt", sep="\t", col.names = NA, row.names = TRUE)

write.table(df_Fibs, "./Results/Fibs_PeakLinkage_.txt", sep="\t", col.names = NA, row.names = TRUE)
write.table(df_iCM, "./Results/iCM_PeakLinkage_.txt", sep="\t", col.names = NA, row.names = TRUE)

# Input to DORC.R to compute DORCs

# Using Joint Analysis significant peak-gene linkage to perform Motif Enrichment Analysis
# iCM.old.Linkage <- read.table("./Results/iCM_PeakLinkage.txt", sep = '\t',header = T)
iCM_Linkage <- read.table("./Results/iCM_PeakLinkage_.txt", sep = '\t',header = T)
# 2022-08-11: Difference between iCM_PeakLinkage_.txt and iCM_PeakLinkage.txt could be FindMarkers for
# _ is between ("2","12") and ("3","11"), while later was between ("2","12") and ("0",3","11"). 

gene <- iCM_Linkage$gene
peakID <- iCM_Linkage$peak
chromosome <- iCM_Linkage$seqnames
start <- iCM_Linkage$start
end <- iCM_Linkage$end
strand <- iCM_Linkage$strand

dup.linkage <- data.frame(gene, peakID, chromosome, start, end, strand)
dedup.linkage <- list()

# trim_dedup is generated from DORC.R line 19
for (i in 1:length(trim_dedup$genes)){
  # print(i)
  for (j in 1: length(dup.linkage$gene)){
    if (trim_dedup[i,1] == dup.linkage[j,1]){
      dedup.linkage <- rbind(dedup.linkage,dup.linkage[j,])
    }
  }
}

write.table(dedup.linkage, "./Coembedding/Results/dedupLinkage.txt", sep="\t", 
            quote=FALSE, col.names = FALSE, row.names = FALSE)

## Output for Homer in Linux
# ssh yqian9@longleaf.unc.edu
# scp dedupLinkage.txt yqian9@longleaf.unc.edu:/nas/longleaf/home/yqian9/R
# 
# cut -f 3,4,5 dedupLinkage.txt > Linkage.bed
# 
# perl /nas/longleaf/home/yqian9/homer/bin/findMotifsGenome.pl /nas/longleaf/home/yqian9/R/Linkage.bed mm10 /nas/longleaf/home/yqian9/R/MotifOutputwparse -size given -preparse
# 
# scp -r yqian9@longleaf.unc.edu:/nas/longleaf/home/yqian9/R/MotifOutputwparse ./
#   
# # Install Homer
#   
# scp configureHomer.pl yqian9@longleaf.unc.edu:/nas/longleaf/home/yqian9/homer
# perl /nas/longleaf/home/yqian9/homer/configureHomer.pl -install
# perl /nas/longleaf/home/yqian9/homer/configureHomer.pl -list
# perl /nas/longleaf/home/yqian9/homer/configureHomer.pl -install mm10
#
# cd /nas/longleaf/home/yqian9/R
# 
# cd /nas/longleaf/home/yqian9/homer/bin

# link peaks to genes
## Ryr2
multi.Rawlins.subset <- LinkPeaks(object = multi.Rawlins.subset, peak.assay = "peaks",
                               expression.assay = "RNA",
                               genes.use = c("Ryr2"))

pdf(file = "./Plots/CoveragePlot_Ryr2.pdf", width = 12, height = 12)
CoveragePlot(object = multi.Rawlins.subset, region = "Ryr2", features = "Ryr2", expression.assay = "RNA",
             extend.upstream = 500, extend.downstream = 10000)
dev.off()

# FeaturePlot for DORCs
DefaultAssay(multi.Rawlins.subset) <- "peaks"
multi.Rawlins.subset <- RegionStats(multi.Rawlins.subset, genome = BSgenome.Mmusculus.UCSC.mm10)
DORCs <- c("Hand2","Atf7","Fos","Foxo1","FOXA1","Sox21","Sox6","Nkx2-5","Smad3","Nkx6-1",
           "Smad4","Mef2c","Tbx5","Gata4")
peakLinkage <- LinkPeaks(object = multi.Rawlins.subset, peak.assay = "peaks",
                         expression.assay = "RNA",
                         genes.use = DORCs)

CoveragePlot(object = peakLinkage, region = "Smad3", features = "Smad3", expression.assay = "RNA",
             extend.upstream = 500, extend.downstream = 10000)

pdf(file = "./Plots/FeaturePlot_DORCs.pdf", width = 12, height = 12)
DefaultAssay(multi.Rawlins.subset) <- "RNA"
FeaturePlot(multi.Rawlins.subset, reduction = "wnn.umap", 
            features = DORCs)
dev.off()


