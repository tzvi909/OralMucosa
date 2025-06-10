

###libs
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)
library(SingleCellExperiment)
library(clustree)

git_dir <- "~/OralMucosa/VU40T_analysis"

### funcs
## don't use for integration

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

VU40T.combined<- readRDS("VU40T_combined_joined_1.0_res_mouse.rds")
VU40T.combined <- NormalizeData(VU40T.combined)
VU40T.combined <- ScaleData(VU40T.combined)
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


p1 <- DimPlot(VU40T.combined, reduction = "tsne", split.by = "condition", label = T)
p2 <- DimPlot(VU40T.combined, reduction = "umap", split.by = "condition", label = T) 
png(file = file.path(git_dir, "Integrated/Mouse/Plots/16_PCs_VU40T_UMAP_tSNE_Resolution_1.png"), width = 16, height = 10, units = "in", res = 300)
print(p1 + p2)
dev.off()

print(p2)
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
library(irGSEA)
##species Homo sapiens or Mus musculus
VU40T.combined_GSEA <- irGSEA.score(object = VU40T.combined, assay = "RNA", 
                             slot = "data", seeds = 123, ncores = 8,
                             min.cells = 3, min.feature = 0,
                             custom = F, geneset = NULL, msigdb = T, 
                             species = "Mus musculus", 
                             category = "H",  
                             subcategory = NULL, geneid = "symbol",
                             method = c("AUCell", "UCell", "singscore", 
                                        "ssgsea", "JASMINE", "viper"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')

