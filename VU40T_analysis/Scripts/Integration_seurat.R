

###libs
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)
library(SingleCellExperiment)
library(clustree)
library(pagoda2) ## GSEA

git_dir <- "~/OralMucosa/VU40T_analysis"

### funcs

format_gene_symbol <- function(gene) {
  # If Ensembl ID (starts with ENS) → leave unchanged
  if (grepl("^ENS", gene)) {
    return(gene)
  }
  
  # If all uppercase → leave unchanged (standard convention for human gene symbols)
  if (toupper(gene) == gene) {
    return(gene)
  }
  
  # Else → capitalize first letter, lowercase rest (standard mouse gene symbol)
  paste0(toupper(substring(gene, 1, 1)), tolower(substring(gene, 2)))
}




preIntegrationSeuratList <- readRDS("annotated_VU40T_sepsamples_singlets_only_mouse.rds")

levels(preIntegrationSeuratList[[4]]$sample)[
  levels(preIntegrationSeuratList[[4]]$sample) == "P22_LPS-N"
] <- "P22_LPS-P"

## Integrate

Anchors <- FindIntegrationAnchors(preIntegrationSeuratList, 
                                  reduction = "cca",  
                                  normalization.method = "LogNormalize",
                                  scale = T
)

VU40T.combined <- IntegrateData(anchorset = Anchors, dims = 1:30)

VU40T.combined$condition <- ifelse(grepl("LPS-N", VU40T.combined$sample), "LPS-N", "LPS-P")
VU40T.combined$passage <- gsub("_.*", "", VU40T.combined$sample)

DefaultAssay(VU40T.combined) <- "integrated"

rm(Anchors)
rm(preIntegrationSeuratList)
VU40T.combined <- ScaleData(VU40T.combined, verbose = F)

VU40T.combined <- RunPCA(VU40T.combined, npcs = 50, verbose = F)

     

#checkpoint
saveRDS(VU40T.combined, file = "VU40T_combined_preRECLUSTER_mouse.rds")

ElbowPlot(VU40T.combined, ndims = 30)

# Run JackStraw (slow - use num.replicate=100 as a perm subsample)
VU40T.combined <- JackStraw(VU40T.combined, num.replicate = 100)

# Score PCs:
VU40T.combined <- ScoreJackStraw(VU40T.combined, dims = 1:20)

# Plot:
png(file = file.path(git_dir, "/Integrated/Mouse/Plots/VU40T_Jackstraw_integrated.png"), width = 16, height = 12, units = "in", res = 300)

JackStrawPlot(VU40T.combined, dims = 1:20) ## jackstraw shows optimum pc of 16
dev.off()
VU40T.combined <- readRDS("VU40T_combined_preRECLUSTER_mouse.rds")
## recluster

optimal_pc <- 15  ## 15 human, 16 mouse - as per jackstraw plot      

VU40T.combined <- FindNeighbors(VU40T.combined, reduction = "pca", dims = 1:optimal_pc)

## find best resolution for integration
resolutions <- seq(0.1, 1.5, by = 0.1)

for (res in resolutions) {
  VU40T.combined <- FindClusters(
    VU40T.combined,
    resolution = res,
    verbose = FALSE, 
    random.seed = 666
  )
  
  # Save clustering results manually
  VU40T.combined[[paste0("RNA_snn_res.", res)]] <- Idents(VU40T.combined)
  
}

VU40T.combined <- FindClusters(
  VU40T.combined,
  resolution = 0.5,
  verbose = FALSE, 
  random.seed = 666
)

png(file = file.path(git_dir, "/Integrated/Mouse/Plots/16_PCs_ClustTree_VU40T_stabilityscores_integrated.png"), width = 16, height = 12, units = "in", res = 300)
p1 <- clustree(VU40T.combined, prefix = "RNA_snn_res.", node_colour = "sc3_stability") ## best resolution lies between 0.6-0.8
print(p1)
dev.off()

png(file = file.path(git_dir, "/Integrated/Mouse/Plots/16_PCs_ClustTree_VU40T_integrated.png"), width = 16, height = 12, units = "in", res = 300)
p1 <- clustree(VU40T.combined, prefix = "RNA_snn_res.") ## best resolution lies between 0.6-0.8
print(p1)
dev.off()


### mouse optimum resolution: now indicated as 1 from trees (Seurat default), human optimum resolution = 0.5

# VU40T.combined <- FindClusters(VU40T.combined, resolution = 0.4)

saveRDS(VU40T.combined, file = "VU40T_combined_RECLUSTER_preUMAP.rds")



VU40T.combined <- RunUMAP(VU40T.combined, reduction = "pca", dims = 1:optimal_pc, verbose = F,  seed.use = 666)
VU40T.combined <- RunTSNE(VU40T.combined, reduction = "pca", dims = 1:optimal_pc, verbose = F, seed.use = 666)
saveRDS(VU40T.combined, file = "VU40T_combined_RECLUSTER_postUMAP.rds")

Idents(VU40T.combined) <- "RNA_snn_res.1" ## 0.5 for human, 1 for mouse
# DefaultAssay(VU40T.combined) <- "RNA"
VU40T.combined[["RNA"]] <- JoinLayers(VU40T.combined[['RNA']])
VU40T.combined$seurat_clusters <- VU40T.combined$RNA_snn_res.1 ##update with idents



Idents(VU40T.combined) <- "RNA_snn_res.0.5" ## 0.5 for human, 1.0 for mouse
DefaultAssay(VU40T.combined) <- "RNA"
VU40T.combined <- JoinLayers(VU40T.combined)
VU40T.combined$seurat_clusters <- VU40T.combined$RNA_snn_res.0.5 ##update with idents

saveRDS(VU40T.combined, file = "VU40T_combined_joined_1.0_res_mouse.rds")

##mouse
VU40T.combined<- readRDS("VU40T_combined_joined_1.0_res_mouse.rds")
##human
VU40T.combined<- readRDS("~/R/VU40T_joined_res0.5.rds")
# VU40T.combined <- NormalizeData(VU40T.combined)
# VU40T.combined <- ScaleData(VU40T.combined)
##find markers
markers <- FindAllMarkers(VU40T.combined,
                          only.pos = T,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)
sig_markers_0.5 <- markers[markers$p_val_adj <= 0.05, ]
write.csv(markers, file = file.path(
  git_dir, 
  "/Integrated/Mouse/15_PCs_ClusterMarkers_res_0.5_again_VU40T_combined.csv"), 
  row.names = FALSE)

## convert back to camelcase if mouse.

markers$gene <- sapply(markers$gene, stringr::str_to_title)



p1 <- DimPlot(VU40T.combined, reduction = "tsne", split.by = "condition", label = T)
p2 <- DimPlot(VU40T.combined, reduction = "umap", split.by = "condition", label = T) 
png(file = file.path(git_dir, "Integrated/Mouse/Plots/16_PCs_VU40T_UMAP_tSNE_Resolution_1.png"), width = 16, height = 10, units = "in", res = 300)
print(p1 + p2)
dev.off()

# print(p2)
library(SingleR) ## doesn't work with human

ref <- readRDS("~/R/HPCA_reference.rds") 
cluster_labels <- VU40T.combined$seurat_clusters
bp <- BiocParallel::MulticoreParam(workers = 8) 
# Run SingleR annotation
sce <- Seurat::as.SingleCellExperiment(VU40T.combined)

### doesn't work with human-aligned integrated. kinda works with mouse
pred <- SingleR::SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main,
  clusters =  cluster_labels,
  assay.type.test = 1,
  BPPARAM=bp
)



# Apply cluster annotations to each cell
VU40T.combined[["SingleR_HTA_cluster_label"]] <- pred$labels[
  match(VU40T.combined$seurat_clusters, rownames(pred))
]

png(file = file.path(git_dir, "Integrated/Mouse/Plots/16_PCs_SingleR_annotated_VU40T_combined_res_1_labelMain.png"), width = 10, height = 8, units = "in", res = 300)

# Visualize cell type annotation on UMAP
p <- Seurat::DimPlot(VU40T.combined, group.by='SingleR_HTA_cluster_label', split.by = "condition")
p2 <- DimPlot(VU40T.combined, reduction = "umap", split.by = "condition", label = T) 
print(p + p2)
# print (p)
dev.off()


###GSEA

### this requires a lot of memory - do not run.
# library(irGSEA)
# ##species Homo sapiens or Mus musculus
# VU40T.combined_GSEA <- irGSEA.score(object = VU40T.combined, assay = "RNA", 
#                              slot = "data", seeds = 123, ncores = 8,
#                              min.cells = 3, min.feature = 0,
#                              group.by = "seurat_clusters",
#                              custom = F, geneset = NULL, msigdb = T, 
#                              species = "Mus musculus", 
#                              category = "H",  
#                              subcategory = NULL, geneid = "symbol",
#                              method = c("AUCell", "UCell", "singscore", 
#                                         "ssgsea", "JASMINE", "viper"),
#                              aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
#                              kcdf = 'Gaussian')

library(msigdbr)
library(clusterProfiler)
library(dplyr)

# Load REACTOME
go_mouse_bp <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol)

# Needed format for clusterProfiler GSEA:
go_mouse_bp_list <- split(go_mouse_bp$gene_symbol, go_mouse_bp$gs_name)



# List of unique clusters
clusters <- unique(markers$cluster)

# Initialize list to store results
gsea_results_list <- list()

for (clust in clusters) {
  cat("Running GSEA for cluster:", clust, "\n")
  
  # Subset markers for this cluster
  markers_clust <- markers %>% filter(cluster == clust)
  
  # Prepare ranked gene list
  gene_list <- markers_clust$avg_log2FC
  names(gene_list) <- markers_clust$gene
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Run GSEA
  gsea_result <- GSEA(
    geneList = gene_list,
    TERM2GENE = go_mouse_bp,
    verbose = FALSE
  )
  
  # Save result
  gsea_results_list[[as.character(clust)]] <- gsea_result
}

# Loop over all GSEA results and plot only those with significant pathways
for (clust in names(gsea_results_list)) {
  cat("Plotting GSEA for Cluster", clust, "...\n")
  
  gsea_result <- gsea_results_list[[clust]]
  
  if (nrow(gsea_result@result) > 0) {
    png(
      filename = file.path(
        git_dir, paste0("/Integrated/Mouse/Plots/GSEA_reactome_cluster_", clust, ".png")
      ),
      width = 5, height = 10, units = "in", res = 300
    )
    print(
      dotplot(gsea_result, showCategory = 20) + 
        ggtitle(paste("Cluster", clust, "GSEA REACTOME"))
    )
  } else {
    message("No significant pathways for Cluster ", clust, " — skipping plot.")
  }
}


## human
go_human_RM <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol)

# Needed format for clusterProfiler GSEA:
go_human_RM_list <- split(go_human_RM$gene_symbol, go_human_RM$gs_name)


  
# List of unique clusters
clusters <- unique(markers$cluster)

# Initialize list to store results
gsea_results_list <- list()

for (clust in clusters) {
  cat("Running GSEA for cluster:", clust, "\n")
  
  # Subset markers for this cluster
  markers_clust <- markers %>% filter(cluster == clust)
  
  # Prepare ranked gene list
  gene_list <- markers_clust$avg_log2FC
  names(gene_list) <- markers_clust$gene
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Run GSEA
  gsea_result <- GSEA(
    geneList = gene_list,
    TERM2GENE = go_human_RM,
    verbose = FALSE
  )
  
  # Save result
  gsea_results_list[[as.character(clust)]] <- gsea_result
}

# Loop over all GSEA results and plot only those with significant pathways
for (clust in names(gsea_results_list)) {
  cat("Plotting GSEA for Cluster", clust, "...\n")
  
  gsea_result <- gsea_results_list[[clust]]
  
  if (nrow(gsea_result@result) > 0) {
    png(
      filename = file.path(
        git_dir, paste0("/Integrated/Plots/GSEA_reactome_cluster_", clust, ".png")
      ),
      width = 5, height = 10, units = "in", res = 300
    )
    print(
      dotplot(gsea_result, showCategory = 20) + 
        ggtitle(paste("Cluster", clust, "GSEA REACTOME"))
    )
  } else {
    message("No significant pathways for Cluster ", clust, " — skipping plot.")
  }
}

