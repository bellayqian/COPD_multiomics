library(Seurat)
library(patchwork)
library(ggplot2)

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

## Run Functions
# multi.Rawlins <- readRDS("./RData/multi.Rawlins_5.rds")
DefaultAssay(multi.Rawlins) <- "RNA"

features_Epithelial <- c('SOX9','ETV5','ASCL1','NEUROD1','SPDEF','TP63')
pdf(file = "./Plots/StackedVlnPlots_Epithelial_Markers.pdf" , width = 7, height = 7)
StackedVlnPlot(obj = multi.Rawlins, features = features_Epithelial)
dev.off()

# features_Mesenchymal <- c('NOTUM')
# pdf(file = "./Plots/StackedVlnPlots_Mesenchymal_Markers.pdf" , width = 7, height = 7)
# StackedVlnPlot(obj = multi.Rawlins, features = features_Mesenchymal)
# dev.off()

features_Endothelial <- c('PROX1')
pdf(file = "./Plots/StackedVlnPlots_Endothelial_Markers.pdf" , width = 7, height = 7)
StackedVlnPlot(obj = multi.Rawlins, features = features_Endothelial)
dev.off()

## Feature plot
features_all <- c('SOX9','ETV5','ASCL1','NEUROD1','SPDEF','TP63','PROX1')
pdf(file = "./Plots/featurePlot_Markers.pdf" , width = 15, height = 15)
FeaturePlot(multi.Rawlins, features = features_all)
dev.off()


####################################################################################
#### Check gene name in multi.Rawlins object
# Check the dimensions of the multi.Rawlins object
dim(multi.Rawlins)

# Retrieve the gene names
gene_names <- rownames(multi.Rawlins)

# Print the first few gene names
head(gene_names)

features_Epithelial <- c('SOX9','ETV5','SCGB1A1','TESC','TPPP3','STC1','PDPN',
                         'HOPXLO','CYTL1','PCP4','SCGB3A','GRP','CHGA','SYP',
                         'ASCL1','TTR','GHRL','CHGA','SYP','NEUROD1','SCGB3A2',
                         'SCGB1A1','SCGB3A1','SPDEF','MUC16','TP63','F3','SPOCK2',
                         'SFTPC','SFTPC','SFTPA','NAPSA')
features_Mesenchymal <- c('SFRP2','PI16','AGTR2','S100A4','WNT2','FGFR4','CXCL14',
                          'KCNK17','PDGFRA','THBD','NOTUM')
features_Endothelial <- c('THY1','CD24','GJA5','DKK2','PTGIS','PVLAP','ACKR3',
                          'HDAC9','PROX1','STAB1','UCP2LO','CA4LO','S100A3')

# Check if the feature genes are present in the gene names
features_Epithelial %in% gene_names
features_Mesenchymal %in% gene_names
features_Endothelial %in% gene_names

# Get the indices of the feature genes in the multi.Rawlins object
feature_indices <- match(features_Epithelial, gene_names)

# Print the indices (NA indicates the gene is not found)
print(feature_indices)

# Subset the multi.Rawlins object to keep only the feature genes
multi.Rawlins_subset <- multi.Rawlins[feature_indices, ]

# Remove any rows with NA values (i.e., genes not found)
multi.Rawlins_subset <- multi.Rawlins_subset[complete.cases(multi.Rawlins_subset), ]

# Create the stacked violin plot using the subsetted object
StackedVlnPlot(obj = multi.Rawlins_subset, features = rownames(multi.Rawlins_subset))

