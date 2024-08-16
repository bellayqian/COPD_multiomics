## GO Enrichment Analysis for Cluster 7,12
library(clusterProfiler)
library(org.Mm.eg.db)

################################################################################
# enrichGO for each cluster
pdf(file = "../Plots/AllClusters_GOplot_BP.pdf", width = 12, height = 10)
for (i in 0:19) {
  fileName <- paste0("../Results/Cluster", i, "_Marker.txt")
  df <- read.table(fileName, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  geneList <- df$gene
  # Perform enrichGO analysis
  enrichResult <- enrichGO(gene          = geneList,
                           OrgDb         = org.Hs.eg.db,
                           keyType       = 'SYMBOL',
                           ont           = "BP", # For Biological Processes
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05)
  # Prepare data for plotting
  top10BP <- enrichResult@result %>%
    mutate(LOG10 = -log10(p.adjust)) %>%
    arrange(desc(LOG10)) %>%
    slice(1:20)
  
  p <- ggplot(top10BP, aes(x=LOG10, y=reorder(Description, LOG10), fill=LOG10)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(name = "-log(P-value)", low="#0066CC", high="#F8766D") +
    labs(title=paste("Cluster", i, "Top 10 Biological Processes"),
         x = "", y = "") +
    theme_classic() +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  print(p)
}
# barplot(enrichResult, showCategory=20
dev.off()