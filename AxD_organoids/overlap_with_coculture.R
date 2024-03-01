library("Seurat")
library("dplyr")
library("clusterProfiler")
library("reshape2")
library("ggplot2")
library("viridis")
library("gridExtra")


# overlap with co-culture data
# precursors = broad markers, log2FC > 1, adj p-val < 0.05
# iA subclusters = markers log2FC > 0.65, adj p-val < 0.05 => filtered additionally to log2FC > 1
# iN subclusters = only few markers above log2FC 1, so they are above 0.65

# organoid markers - broad signature defining clusters as they are as well as possible: filter all for log2FC > 0.25, adj p-val < 0.05

# create list of marker genes
create_gene_list <- function(gene.lists.path){
  dir <- dir(gene.lists.path, pattern = ".txt")
  ref.list <- list()
  for (i in 1:length(dir)){
    ref.list[[gsub(".txt","",dir[i])]] <- as.vector(read.table(paste0(gene.lists.path,dir[i]), header = F)[,1])  
  }
  return(ref.list)
}

# addmodulescore for UMAP projection
custom_addmodule <- function(so, ref.genes, plot.name = NULL){
  print("AddModuleScore running...")
  ref.genes <- head(ref.genes, n=10)
  print(ref.genes)
  seurat.scored <- AddModuleScore(so, features = list(ref.genes), assay = "RNA",
                                  name = "reference_score_")
  p <- FeaturePlot(seurat.scored, features = "reference_score_1", 
                   min.cutoff = "q10", max.cutoff = "q90", reduction = "umap", cols = c("#e0e0e0","#8c0052")) + 
    ggtitle(plot.name) + NoAxes()
  return(p)
  print("DONE")
}

# co-culture cluster markers in a list
dataset1_clusters <- create_gene_list("co-culture_overlap/")
names(dataset1_clusters) <- lapply(names(dataset1_clusters), function(x) gsub("_", " ", x))

#RUN ONCE write marker lists into txt, pick them out from csv table
# write_txt_markers <- function(markers, directory){
#   for (cluster.id in levels(as.factor(markers$cluster))){
#   write(x = markers %>% subset(cluster == cluster.id & avg_log2FC > 0.25 & p_val_adj < 0.05) %>% arrange(desc(avg_log2FC)) %>% rownames,
#         file = paste0(directory,"/markers_txt/",gsub("/","_",cluster.id),".txt"))
#   }
# }
# 
# write_txt_markers(read.csv("cerebral/markers_cerebral_namedclusters.csv", header = T, row.names = 1),
#                   "cerebral")
# write_txt_markers(read.csv("cortical_integrated/markers_integrated.csv", header = T, row.names = 1),
#                   "cortical_integrated")

# load co-culture Seurat to set background genes for ORA
seurat.object <- readRDS("seurat_all_merged_filtered_mt_rib.rds") # co-culture dataset
seurat.object@active.assay <- "RNA"
background <- rownames(seurat.object)

organoids <- c("cerebral","cortical_integrated")

for (org in organoids){
  dataset2_clusters <- create_gene_list(paste0(org,"/markers_txt/"))
  
  # prepare organoid markers in Term2Gene table
  Term2Gene <- data.frame(ID = character(0),
                          Gene = character(0))
  for (ref in 1:length(dataset2_clusters)){
    ref.df <- data.frame(ID = rep(names(dataset2_clusters)[ref], length(dataset2_clusters[[ref]])),
                         Gene = dataset2_clusters[[ref]])
    Term2Gene <- rbind(Term2Gene,ref.df)
  }
  
  # perform enrichment analysis
  enr.result <- lapply(dataset1_clusters, function(x) enricher(gene = x,
                                                               universe = background,
                                                               pvalueCutoff = 1,
                                                               minGSSize = 1,
                                                               maxGSSize = 1000,
                                                               qvalueCutoff = 1,
                                                               TERM2GENE = Term2Gene))
  
  enr.res.merge <- merge_result(enr.result)
  saveRDS(enr.res.merge, paste0("co-culture_overlap/", org, "_ORAresult.rds"))
  
  pval.matrix <- matrix(data=NA, nrow = length(dataset1_clusters), ncol = length(dataset2_clusters), 
                        dimnames = list(names(dataset1_clusters),names(dataset2_clusters)))
  
  for (i in 1:nrow(enr.res.merge@compareClusterResult)){
    selection <- enr.res.merge@compareClusterResult[i,]
    pval.matrix[selection$Cluster,selection$ID] <- enr.res.merge@compareClusterResult[i,"p.adjust"]
  }
  
  # correct matrix, so that it can be plotted in a "heatmap"
  #p-values that were not calculated because of lack of overlap are 1 and number of genes enriched is 0
  pval.matrix[is.na(pval.matrix)] <- 1
  
  pval.df <- melt(pval.matrix)
  count <- c()
  for (j in 1:nrow(pval.df)){
    if (length(which(enr.res.merge@compareClusterResult$Cluster == pval.df$Var1[j] & enr.res.merge@compareClusterResult$ID == pval.df$Var2[j])) != 0) {
      count.num <- enr.res.merge@compareClusterResult[enr.res.merge@compareClusterResult$Cluster == pval.df$Var1[j] & enr.res.merge@compareClusterResult$ID == pval.df$Var2[j], "Count"]
    } else {
      count.num <- 0
    }
    count <- append(count, count.num)
  }
  
  pval.df$count <- count 
  
  if (org == "cerebral"){
    levels(pval.df$Var2) <- c("cardiac/skeletal muscle cells","choroid plexus cells"  ,"astrocytes","epithelial cells","immature neurons NEUROD2+ FABP7+","immature neurons","intermediate progenitors","mesenchymal-like cells/fibroblasts IGFBP3+ COL15A1+","mesenchymal-like cells/fibroblasts","mesothelial cells","neuronal+mesenchymal-like cells", "pre-OPCs","pancreatic-like cells","peripheral neurons","proliferating radial glia","proliferating satellite cells","radial glia","satellite cells","upper layer neurons", "stressed VEGFA+ cells")
  } else if (org == "cortical_integrated"){
    levels(pval.df$Var2) <- c("CRABP1+ cells","pre-OPCs","intermediate progenitors","mesenchymal-like cells","muscle cells","neurons 1","neurons 2","oRG/astrocytes","proliferating mesenchymal-like cells","proliferating radial glia","radial glia","stressed VEGFA+ cells")                                                            
  }
  
  # ggplot "heatmap" showing number of enriched genes in each overlap and significance of this overlap. adj p-value is relative to each row - each enrichment analysis
  p <- ggplot(pval.df, aes(x = Var2, y = factor(Var1, levels = c("precursors", "ASTRO 1","ASTRO 2","ASTRO 3","ASTRO 4","NEURO 1","NEURO 2","NEURO 3","NEURO 4")), fill = value)) +
      geom_tile() +
      geom_text(aes(label = count), color = "black", size = 4) +
      coord_fixed() + 
      scale_fill_gradient(low = "#09B852", high = "#e7e7e7") +
      #scale_fill_viridis(direction = -1, begin = 0.5, end = 1) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1, size = 12),
            axis.text.y = element_text(size = 12)) +
      xlab("Organoid clusters") + ylab("Co-culture clusters") +
      guides(fill = guide_colourbar(title= bquote("p"["adj"]), reverse = T))
  p
  assign(paste0("ht_",org),p)
  ggsave(paste0("co-culture_overlap/", org, "_heatmap2.png"), plot = p, device = "png", bg = "white", width = 9, height = 6)
  
  # co-culture markers projected on organoid UMAPs
  # load organoid seurat objects
  
  if (org == "cerebral") {
    seurat.org <- readRDS("cerebral/seurat_subset_cerebral.rds")  
  } else if (org == "cortical_integrated"){
    seurat.org <- readRDS("cortical_integrated/seurat_cortical_integrated.rds")  
  }
  
  dataset1_clusters
  scored.plots <- lapply(seq_along(dataset1_clusters), function(x) custom_addmodule(seurat.org, dataset1_clusters[[x]],
                                                                           plot.name = names(dataset1_clusters)[x]))
  
  names(scored.plots) <- names(dataset1_clusters)
  g <- grid.arrange(grobs = scored.plots[c("precursors","ASTRO 1","ASTRO 2","NEURO 1","NEURO 2","NEURO 4")])
  g
  assign(paste0("umap_",org),g)
  ggsave(paste0("co-culture_overlap/", org, "_UMAP.png"), plot = g, device = "png", height = 8, width = 6)
}


# dimplot with cluster labels to make orientation in the supp figure easier
ggsave("co-culture_overlap/UMAP_cerebral.png", DimPlot(seurat.org, label = T, repel = T) + NoLegend() + NoAxes(), device = "png", width = 6, height = 6)


