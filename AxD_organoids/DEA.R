#DEA
library(Seurat)
library(ggplot2)
library(xlsx)
library(EnhancedVolcano)
library(ggrepel)
library(org.Hs.eg.db)
library(clusterProfiler)
library(xlsx)
library(enrichplot)
library(dplyr)



DEtests <- function(
    seuratObject,
    conditions,
    test = "t",
    #treshold = 0.05,
    verbose = TRUE,
    ...
) {
  DefaultAssay(seuratObject) <- "RNA"
  print(seuratObject)
  dfs <- list()
  for (i in conditions){
    DEGs <- FindMarkers(seuratObject, ident.1 = i[2], ident.2 = i[1], test.use = test, logfc.threshold = 0, min.pct = 0.1)
    DEGs <- DEGs[order(DEGs$avg_log2FC, decreasing = FALSE),]
    dfs[[paste0(i[1],"-",i[2])]] <- DEGs
  }
  return(dfs)
}

file.name <- "DEA"

seurat.object <- readRDS("cortical_integrated/seurat_oRG_forDEA.rds")
seurat.object <- seurat.DEA

levels(as.factor(seurat.object@meta.data[["Sample"]]))
Idents(seurat.object) <- seurat.object@meta.data[["Sample"]]

# adjustable
conditions <- list(levels(seurat.object)) 
DE_tests_name <- paste0(conditions[[1]][1],"-",conditions[[1]][2])

DE_tests <- DEtests(seurat.object, conditions)
saveRDS(DE_tests, paste0(file.name,DE_tests_name, ".rds"))
  
file.xlsx <- paste0(file.name,DE_tests_name, "_DEGs.xlsx")
write.xlsx("empty", file.xlsx, sheetName = "1", append = FALSE)
for (test in names(DE_tests)){
  top.genes <- DE_tests[[test]][DE_tests[[test]]$p_val_adj < 0.05,]
  top.genes <- top.genes[order(top.genes$p_val_adj),]
  if (nrow(top.genes) == 0){
    write.xlsx("no genes <0.05 found", file.xlsx, sheetName = test, append = TRUE)
  }
  else {
    write.xlsx(top.genes, file.xlsx, sheetName = test, append = TRUE)
    res <- DE_tests[[test]]
    
    keyvals.colour <- ifelse(
      res$avg_log2FC < -0.65 & res$p_val_adj < 0.05, '#5A7BC6',
      ifelse(res$avg_log2FC > 0.65 & res$p_val_adj < 0.05, '#EB7321',
             '#E7E7E7')) 
             names(keyvals.colour)[keyvals.colour == '#5A7BC6'] <- 'downregulated'
             names(keyvals.colour)[keyvals.colour == '#EB7321'] <- 'upregulated'
             names(keyvals.colour)[keyvals.colour == '#E7E7E7'] <- 'insignificant'
             
   labs <- c(res[res$avg_log2FC < -0.65 & res$p_val_adj < 0.05,] %>% arrange(p_val_adj) %>% head(10) %>% rownames(),
             res[res$avg_log2FC > 0.65 & res$p_val_adj < 0.05,] %>% arrange(p_val_adj) %>% head(10) %>% rownames())

   p <- EnhancedVolcano(toptable = res,
                        lab = rownames(res),
                        selectLab = labs,
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        axisLabSize = 8,
                        titleLabSize = 10,
                        title = "CER_cr-CER_AxD",
                        pCutoff = 0.05, 
                        colAlpha = 0.7,
                        colCustom = keyvals.colour,
                        subtitleLabSize = 9,
                        FCcutoff = 0.65, 
                        labSize = 3,
                        pointSize = 2.5,
                        drawConnectors = T,
                        maxoverlapsConnectors = Inf,
                        min.segment.length = 0.1,
                        arrowheads = FALSE,
                        boxedLabels = F,
                        cutoffLineType = "blank",
                        caption = "|log2FC| > 0.65 and p-adj < 0.05") +
     theme_light() +
     theme(plot.subtitle = element_blank(), 
           plot.caption = element_text(size = 10), 
           plot.title = element_text(hjust = 0.5, size = 12),
           panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
           axis.line = element_line(colour = "black", linewidth = 0.7, linetype = "solid"), 
           axis.text = element_text(size = 12),
           axis.title = element_text(size = 12),
           axis.ticks = element_line(linewidth = 0.7, colour = "black"),
           axis.ticks.length = unit(0.2, "cm"),
           legend.position = "none") +
     xlab(expression("Log"[2]*"FC")) +
     guides(color = guide_legend(override.aes = list(size = 4)))
    ggsave(plot = p, paste0(file.name,DE_tests_name,"_volcano.png"), device = "png", 
           height = 3.5, width = 4.5)
  }
}
  
#### ORA
seurat.object <- readRDS("cortical_integrated/seurat_oRG_forDEA.rds")
seurat.object@active.assay <- "RNA"; seurat.object@active.ident <- as.factor(seurat.object$Condition)
table(seurat.object@active.ident)
DE.object <- readRDS("cortical_integrated/DEA_oRG/Ctx_cr-Ctx_AxD.rds")
DE.object <- DE_tests
file.name <- "cortical_integrated/DEA_oRG/"

background.genes <- rownames(seurat.object)
DEA.list <- list(DE.object[[1]][DE.object[[1]]$avg_log2FC > 0.65 & DE.object[[1]]$p_val_adj < 0.05,],
                 DE.object[[1]][DE.object[[1]]$avg_log2FC < -0.65 & DE.object[[1]]$p_val_adj < 0.05,])
sapply(DEA.list, function(x) nrow(x))
names(DEA.list) <- c("UP","DOWN")
ego.results <- list()

for (gene.list in 1:length(DEA.list)){    
  ego <- enrichGO(gene       = rownames(DEA.list[[gene.list]]),
                  universe     = background.genes,
                  keyType       = "ALIAS",
                  OrgDb         = org.Hs.eg.db,
                  ont           = "ALL",
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 0.1,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE,
                  minGSSize     = 3)
  ego <- simplify(ego)
  ego.results[[names(DEA.list)[gene.list]]] <- ego
}  
saveRDS(ego.results, paste0(file.name,"ORA_results.rds"))

file.xlsx <- paste0(file.name,"ORA_results_sheet.xlsx")
write.xlsx("empty", file.xlsx, sheetName = "0", append = FALSE)
for (i in 1:length(ego.results)){
  if (dim(ego.results[[i]])[1] != 0) {
    write.xlsx(ego.results[[i]], file.xlsx, sheetName = names(ego.results)[i], append = TRUE)}
}
  
out <- merge_result(ego.results)
out@compareClusterResult <- out@compareClusterResult[order(out@compareClusterResult$p.adjust),]
out@compareClusterResult <- out@compareClusterResult[out@compareClusterResult$Count > 3,]
p <- dotplot(out, showCategory = 5, by = "GeneRatio", includeAll = F, label_format = 35) +
  theme_minimal() +
  theme(axis.text = element_text(size = 16, colour = "black"),
        axis.text.x = element_text(face = "bold"),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 16, face = "bold", hjust = 0)) + 
  ggtitle("GO Terms Enrichment") +
  scale_color_viridis(begin = 0.6, direction = -1)
p
ggsave(plot = p, filename = paste0(file.name,"ORA_plot",".png"), device = "png", bg = "white",
       width = 7, height = 4)
  
