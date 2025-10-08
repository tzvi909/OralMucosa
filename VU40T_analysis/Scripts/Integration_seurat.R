Sys.setenv(BIOCFILECACHE_DIR = tempfile())

###libs
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)
library(SingleCellExperiment)
library(clustree)
library(biomaRt)
library(msigdbr)
library(enrichplot)
library(clusterProfiler)
library(ComplexHeatmap)
library(circlize)  # for colors)
library(RColorBrewer)
library(optparse)

###---Opts----

option_list <- list(
  make_option(c("-s", "--species"),
              type = "character",
              default = "human",
              help = "Species to use [default = %default]. Options: human or mouse",
              metavar = "character"),
  make_option(c("-i", "--input"),
              type = "character",
              default = ,
              help = "input RDS dir to load",
              metavar = "character"),
  make_option(c("-n", "--n_cores"),
              type = "integer",
              default = 1,
              help = "Number of cores to use to run  [default = %default]",
              metavar = "character")
)


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
print(opt)

n_cores <- opt$n_cores

species <- opt$species

## to avoid hitting that 20gb quota on home dir,
proj_dir <- "/rds/projects/g/gendood-3dmucosa/"
analysis_dir <- file.path(proj_dir, "scRNAseqAnalysis/")
git_dir <- file.path(analysis_dir, "OralMucosa/VU40T_analysis")
out_prefix_dir <- file.path(git_dir, "Integrated")

cache_dir <- file.path(proj_dir, "rds_cache")

## check for dirs recursively

chk_dir_list <- list(analysis_dir, git_dir, cache_dir, out_prefix_dir)

for (path in chk_dir_list){
  if(!(dir.exists(path))){
    dir.create(path, recursive = T)
  }
}

### funcs



make_marker_dotplots <- function(seurat_obj, genelist, outPrefix, export = T, plot_dir = "."){
  out_dir <- file.path(plot_dir, "MarkerDotplots")
  if(!(dir.exists(out_dir))){
    dir.create(out_dir, recursive = T)
  } 
  
  colnames(genelist) <- sapply(colnames(genelist), function(x) {
    parts <- stringr::str_split(x, "_")[[1]]
    if (length(parts) >= 2 && parts[2] != toupper(parts[2])) {
      parts[2] <- stringr::str_to_title(parts[2])
    }
    paste(parts, collapse = "_")
  })
  colnames(genelist) <- gsub("[._]", " ", colnames(genelist))
  ## revert to whitespaces if whitespaces in original colname or change underscores to whitespace
  if (all(seurat_obj$species == "Mouse")){
    
    # Connect to Ensembl human mart (homolog info is already stored here)
    human <- useEnsembl(
      biomart = "genes", 
      dataset = "hsapiens_gene_ensembl", 
      host = "https://www.ensembl.org"  # pick a mirror to avoid overload
    )
    
    # Your input gene list
    all_genes <- unique(unlist(genelist))
    
    # ✅ These attributes are from the *same attribute page* (homologs)
    conversion <- getBM(
      attributes = c("external_gene_name", "mmusculus_homolog_associated_gene_name"),
      filters    = "external_gene_name",
      values     = all_genes,
      mart       = human
    )
    
    # Build mapping dict: human → mouse
    human_to_mouse <- setNames(
      toupper(conversion$mmusculus_homolog_associated_gene_name),
      conversion$external_gene_name
    )
    
    # Clean NAs or empty values
    human_to_mouse <- human_to_mouse[human_to_mouse != "" & !is.na(human_to_mouse)]
    # sanity check
    # head(human_to_mouse)
  }
  
  
  for (i in seq_along(genelist)) {
    # Extract non-NA genes from this column
    gene_set <- unique(genelist[[i]][!is.na(genelist[[i]]) & genelist[[i]] != ""])
    # Skip if empty
    if (length(gene_set) == 0) next
    
    # Replace genes if object is mouse and a dictionary is provided
    # Convert human genes to mouse if species is Mouse and a dictionary is provided
    if (all(seurat_obj$species == "Mouse") && !is.null(human_to_mouse)) {
      gene_set <- human_to_mouse[gene_set]                  # map human → mouse
      gene_set <- gene_set[!is.na(gene_set)]           # drop unmatched
      gene_set <- unique(gene_set)                     # remove duplicates
      gene_set <- intersect(gene_set, rownames(seurat_obj))  # keep genes in dataset
    }
    
    # For human data, keep only genes present in dataset
    if (all(seurat_obj$species == "Human")) {
      gene_set <- intersect(gene_set, rownames(seurat_obj))
    }
    
    # Skip if final gene set is empty
    if (length(gene_set) == 0) next
    
    # Generate dot plot
    p <- DotPlot(seurat_obj, features = gene_set, group.by = "seurat_clusters") +
      RotatedAxis() +
      ggplot2::ggtitle(colnames(genelist)[i]) +
      scale_color_gradient(low = "lightgrey", high = "red") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) + 
      labs(x = "Genes", y = "Clusters")  # Relabel axes
    
    # Modify x-axis labels to Title Case if using mouse data
    if (all(seurat_obj$species == "Mouse")){
      p <- p + scale_x_discrete(labels = function(x) tools::toTitleCase(tolower(x)))
    }
    
    
    
    # Export or print
    if (export) {
      
      num_genes <- length(na.omit(gene_set[[i]]))
      if (num_genes <= 6) {
        plot_width <- 6  # Small gene sets: fixed small width
      } else {
        plot_width <- min(max(2 * num_genes + 5, 12), 50)
      }
      colnames(genelist) <- gsub(" ", "_", colnames(genelist))
      
      
      if (all(seurat_obj$species == "Human")) {
        png(file = file.path(
          out_dir ,
          paste0("HumanONLY_Dotplot_clusterResolution0.3_",
                 outPrefix, "_", colnames(genelist)[i], "_VU40T_combined.png")),
          width = plot_width, height = 4, units = "in", res = 300)
      } else {
        png(file = file.path(
          out_dir,
          paste0("Dotplot_clusterResolution0.6_",
                 outPrefix, "_", colnames(genelist)[i], "_VU40T_combined.png")),
          width = plot_width, height = 4, units = "in", res = 300)
      }
      colnames(genelist) <- gsub("[._]", " ", colnames(genelist))
      print(p)
      dev.off()
    } else {
      colnames(genelist) <- gsub("[._]", " ", colnames(genelist))
      print(p)
    }
  }
}

if (is.null(opt$input)){
  preIntegrationSeuratList <- readRDS(file.path(cache_dir, paste0("VU40T_singlets_only_sep_samples_",species,".RDS")))
} else {
  preIntegrationSeuratList <- readRDS(opt$input)
  ### sanity check if using opt$input arg to make sure species matches seurat list RDS
  if (any(sapply(preIntegrationSeuratList, function(obj) any(grepl("^ENSMUS", rownames(obj[["RNA"]])))))) {
    species <- "mouse"
  } else {
    species <- "human"
  }
}

plots_dir <- file.path(out_prefix_dir, paste0(stringr::str_to_title(species), "/Plots"))
res_dir <- file.path(out_prefix_dir, paste0(stringr::str_to_title(species)))

out_dir_list <- list(plots_dir, res_dir)

for (path in out_dir_list){
  if(!(dir.exists(path))){
    dir.create(path, recursive = T)
  }
}

## species_df

species_df <- read.csv(file.path(git_dir, "/VU40T_species_assignment_by_reads.csv"))[,-c(1)]

preIntegrationSeuratList <- readRDS("annotated_VU40T_sepsamples_singlets_only_mouse.rds")

preIntegrationSeuratList <- readRDS("R/VU40T_singlets_only_sep_samples.RDS")

species <- "Human"

levels(preIntegrationSeuratList[[4]]$sample)[
  levels(preIntegrationSeuratList[[4]]$sample) == "P22_LPS-N"
] <- "P22_LPS-P"

## Integrate

## filter out non-aligned species


## add in species data

for (sample in names(preIntegrationSeuratList)) {
  
  # Get the Seurat object for the current sample
  seurat_obj <- preIntegrationSeuratList[[sample]]
  
  # Filter species_df for this sample
  df_sample <- subset(species_df, sample == sample)
  
  # Make sure barcodes match
  df_sample <- df_sample[!duplicated(df_sample$barcode), ]
  rownames(df_sample) <- df_sample$barcode
  
  # Match the species to the Seurat object's barcodes
  matched_species <- df_sample[colnames(seurat_obj), "species_label"]
  
  # Add to metadata
  seurat_obj$species <- matched_species
  
  # Store back in list
  preIntegrationSeuratList[[sample]] <- seurat_obj
}


### filter in automated fashion


##> table(preIntegrationSeuratList[[1]]$species)

##Human Mixed Mouse 
#2004    91  7062 
#> table(preIntegrationSeuratList[[2]]$species)

#Human Mixed Mouse 
#3284   114  4315 
#> table(preIntegrationSeuratList[[3]]$species)

#Human Mixed Mouse 
#2802    86  4916 
#> table(preIntegrationSeuratList[[4]]$species)

#Human Mixed Mouse 
#1441    90  6022 

## for human this will be reversed

## so we remove the one with the least number of cells - this will only work with mouse so we have to be selective

## Cells assigned as "Mixed" seem to be evenly spread around all clusters 
## so should be removed anyway in case of being ambient RNA or doublets not removed in earlier steps

if(species == "Mouse"){
  
  for (i in seq_along(preIntegrationSeuratList)) {
    
    seurat_obj <- preIntegrationSeuratList[[i]]
    
    # Keep only cells from dominant species
    cells_to_keep <- colnames(seurat_obj)[seurat_obj$species == "Mouse"]
    
    # Subset Seurat object
    seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
    
    # Store it back
    preIntegrationSeuratList[[i]] <- seurat_obj
  }
}else{
  for (i in seq_along(preIntegrationSeuratList)) {
    
    seurat_obj <- preIntegrationSeuratList[[i]]
    
    # Keep only cells from dominant species
    cells_to_keep <- colnames(seurat_obj)[seurat_obj$species == "Human"]
    
    # Subset Seurat object
    seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
    
    # Store it back
    preIntegrationSeuratList[[i]] <- seurat_obj
  }
}


## now we can integrate and not worry about mess from cross-species interaction

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
saveRDS(VU40T.combined, file = file.path(cache_dir, "VU40T_combined_preRECLUSTER_human_ONLY.rds"))

# VU40T.combined <- readRDS("VU40T_combined_preRECLUSTER_mouse_ONLY.rds")

ElbowPlot(VU40T.combined, ndims = 50) ## around 15 for mouse around 12 for human


# Run JackStraw (slow - use num.replicate=100 as a perm subsample)
VU40T.combined <- JackStraw(VU40T.combined, num.replicate = 100)

# Score PCs:

VU40T.combined <- ScoreJackStraw(VU40T.combined, dims = 1:20)

# Plot:
# png(file = file.path(git_dir, "/Integrated/Mouse/Plots/Mouse_only_VU40T_Jackstraw_integrated.png"), width = 16, height = 12, units = "in", res = 300)
png(file = file.path(git_dir, "/Integrated/Plots/Human_only_VU40T_Jackstraw_integrated.png"), width = 16, height = 12, units = "in", res = 300)

JackStrawPlot(VU40T.combined, dims = 1:20) ## jackstraw shows optimum pc of 16, 20 when mouse only
dev.off()
VU40T.combined <- readRDS("VU40T_combined_preRECLUSTER_mouse.rds")
## recluster

optimal_pc <- 20  ## 15 human, 16 mouse - as per jackstraw plot when using all data together. 20 for mouse seperately, 10 for human seperately    

optimal_pc <- 10  ## 15 human, 16 mouse - as per jackstraw plot when using all data together. 20 for mouse seperately, 10 for human seperately    

optimal_pc <- 15  ## 15 human, 16 mouse - as per jackstraw plot when using all data together. 20 for mouse seperately, 10 for human seperately    


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
  resolution = 0.6,
  verbose = FALSE, 
  random.seed = 666
)

# VU40T.combined <- FindClusters(
#   VU40T.combined,
#   resolution = 0.5,
#   verbose = FALSE, 
#   random.seed = 666
# )

png(file = file.path(git_dir, "/Integrated/Mouse/Plots/MouseOnly_20_PCs_ClustTree_VU40T_stabilityscores_integrated.png"), width = 16, height = 12, units = "in", res = 300)
p1 <- clustree(VU40T.combined, prefix = "RNA_snn_res.", node_colour = "sc3_stability") ## best resolution lies between 0.6-0.8
print(p1)
dev.off()

png(file = file.path(git_dir, "/Integrated/Mouse/Plots/MouseOnly_20_PCs_ClustTree_VU40T_integrated.png"), width = 16, height = 12, units = "in", res = 300)
p1 <- clustree(VU40T.combined, prefix = "RNA_snn_res.") ## best resolution lies between 0.6-0.8
print(p1)
dev.off()


png(file = file.path(git_dir, "/Integrated/Plots/HumanOnly_20_PCs_ClustTree_VU40T_stabilityscores_integrated.png"), width = 16, height = 12, units = "in", res = 300)
p1 <- clustree(VU40T.combined, prefix = "RNA_snn_res.", node_colour = "sc3_stability") ## best resolution lies between 0.6-0.8
print(p1)
dev.off()

png(file = file.path(git_dir, "/Integrated/Plots/HumanOnly_20_PCs_ClustTree_VU40T_integrated.png"), width = 16, height = 12, units = "in", res = 300)
p1 <- clustree(VU40T.combined, prefix = "RNA_snn_res.") ## best resolution lies between 0.6-0.8
print(p1)
dev.off()

### mouse optimum resolution: now indicated as 1 from trees (Seurat default), human optimum resolution = 0.5 (without species filtering)
### human optimum resolution now: 0.3 when human only
### when mouse only using 20 pcs, optimum resolution is 0.6

# VU40T.combined <- FindClusters(VU40T.combined, resolution = 0.4)

saveRDS(VU40T.combined, file = file.path(cache_dir, "VU40T_combined_RECLUSTER_preUMAP_human.rds"))



VU40T.combined <- RunUMAP(VU40T.combined, reduction = "pca", dims = 1:optimal_pc, verbose = F,  seed.use = 666)
VU40T.combined <- RunTSNE(VU40T.combined, reduction = "pca", dims = 1:optimal_pc, verbose = F, seed.use = 666)
saveRDS(VU40T.combined, file = file.path(cache_dir, "VU40T_combined_RECLUSTER_postUMAP_humanOnly.rds"))

VU40T.combined <- readRDS(file.path(cache_dir, "VU40T_combined_RECLUSTER_postUMAP_humanOnly.rds"))
Idents(VU40T.combined) <- "RNA_snn_res.0.6" ## 0.5 for human, 1 for mouse (combined). 0.6 for mouse Only

Idents(VU40T.combined) <- "RNA_snn_res.0.3" ## 0.5 for human, 1 for mouse (combined). 0.6 for mouse Only, 0.3 for human only

DefaultAssay(VU40T.combined) <- "RNA"
VU40T.combined[["RNA"]] <- JoinLayers(VU40T.combined[['RNA']])
VU40T.combined$seurat_clusters <- Idents(VU40T.combined) ##update with idents



# Idents(VU40T.combined) <- "RNA_snn_res.0.5" ## 0.5 for human, 1.0 for mouse
# DefaultAssay(VU40T.combined) <- "RNA"
# VU40T.combined <- JoinLayers(VU40T.combined)
# VU40T.combined$seurat_clusters <- VU40T.combined$RNA_snn_res.0.5 ##update with idents

saveRDS(VU40T.combined, file = file.path(cache_dir, "VU40T_combined_joined_1_res_humanOnly.rds"))
saveRDS(VU40T.combined, file = file.path(cache_dir, "VU40T_combined_joined_0.3_res_humanOnly.rds"))

##mouse
# VU40T.combined<- readRDS("VU40T_combined_joined_1.0_res_mouse.rds")
##human
VU40T.combined<- readRDS("~/R/VU40T_joined_res0.5.rds")
##VU40T.combined<- readRDS(file.path(cache_dir, "VU40T_combined_joined_1_res_humanOnly.rds")) ## overclustered

VU40T.combined <- readRDS(file.path(cache_dir, "VU40T_combined_joined_0.3_res_humanOnly.rds"))
# VU40T.combined <- NormalizeData(VU40T.combined)
# VU40T.combined <- ScaleData(VU40T.combined)
##find markers
##mouse only
VU40T.combined <- readRDS(file.path(cache_dir, "VU40T_combined_joined_0.6_res_mouseOnly.rds"))

if (all(VU40T.combined$species == "Mouse")){
  Species <- "Mouse"
} else{
  Species <- "Human"
}

markers <- FindAllMarkers(VU40T.combined,
                          only.pos = T,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)
sig_markers_0.3 <- markers[markers$p_val_adj <= 0.05, ]


# 
# markers <- FindAllMarkers(VU40T.combined,
#                           only.pos = F,
#                           min.pct = 0,
#                           logfc.threshold = 0)
# 
# c6_markers <- markers[markers$cluster == 6, ]

# c6_markers

write.csv(markers, file = file.path(
  git_dir, 
  "/Integrated/HumanOnly10_PCs_ClusterMarkers_res_0.3_VU40T_combined.csv"), 
  row.names = FALSE)

## convert back to camelcase if mouse.
if (all(VU40T.combined$species == "Mouse")){
  markers$gene <- sapply(markers$gene, stringr::str_to_title)
}




p1 <- DimPlot(VU40T.combined, reduction = "tsne", split.by = "condition", label = T)
p2 <- DimPlot(VU40T.combined, reduction = "umap", split.by = "condition", label = T) 

if (all(VU40T.combined$species == "Mouse")){
  png(file = file.path(git_dir, "Integrated/Mouse/Plots/MouseOnly_20_PCs_VU40T_UMAP_tSNE_Resolution_0.6.png"), width = 16, height = 10, units = "in", res = 300)
}else{
  png(file = file.path(git_dir, "Integrated/Plots/HumanOnly_10_PCs_VU40T_UMAP_tSNE_Resolution_0.3.png"), width = 16, height = 10, units = "in", res = 300)
}

print(p1 + p2)
dev.off()

# print(p2)
library(HCATonsilData)


library("BiocFileCache")

Sys.setenv(BFC_CACHE=file.path(proj_dir, ".cache"))
# Example: Set the cache to a specific directory
bfc = BiocFileCache(cache = file.path(proj_dir, ".cache"), ask = FALSE)
bfccache(bfc)

bfccache(bfc)

tools::R_user_dir("BiocFileCache", which="cache")

if (all(VU40T.combined$species == "Human")){
  sce_HCATonsil_epithelial <- HCATonsilData(assayType = "RNA", cellType = "epithelial")
  # Get counts and metadata manually
  counts_epi <- counts(sce_HCATonsil_epithelial)
  counts_epi <- as(counts_epi, "dgCMatrix")
  logcounts_epi <- logcounts(sce_HCATonsil_epithelial)
  logcounts_epi <- as(logcounts_epi, "dgCMatrix")
  metadata_epi <- as.data.frame(colData(sce_HCATonsil_epithelial))
  seu_epi <- CreateSeuratObject(counts = counts_epi, meta.data = metadata_epi)
  seu_epi <- SetAssayData(seu_epi, slot = "data", new.data = logcounts_epi)
  seu_ref <- seu_epi
  rm(list = c(seu_epi, counts_epi, logcounts_epi, metadata_epi))
}

if (all(VU40T.combined$species == "Mouse")){
  sce_HCATonsil_mesenchymal <- HCATonsilData(assayType = "RNA", cellType = "FDC")
  counts_fdc <- counts(sce_HCATonsil_mesenchymal)
  logcounts_fdc <- logcounts(sce_HCATonsil_mesenchymal)
  counts_fdc <- as(counts_fdc, "dgCMatrix")
  logcounts_fdc <- as(logcounts_fdc, "dgCMatrix")
  metadata_fdc <- as.data.frame(colData(sce_HCATonsil_mesenchymal))
  seu_fdc <- CreateSeuratObject(counts = counts_fdc, meta.data = metadata_fdc)
  seu_fdc <- SetAssayData(seu_fdc, slot = "data", new.data = logcounts_fdc)
  seu_ref <- seu_fdc
  rm(list = c(seu_fdc, counts_fdc, logcounts_fdc, metadata_fdc))
}

seu_ref <- NormalizeData(seu_ref)
seu_ref <- FindVariableFeatures(seu_ref)
seu_ref <- ScaleData(seu_ref)
seu_ref <- RunPCA(seu_ref)

anchors <- FindTransferAnchors(
  reference = seu_ref,
  query = VU40T.combined,
  dims = 1:30,
  reference.assay = "RNA",
  query.assay = "RNA"
)

predictions <- TransferData(
  anchorset = anchors,
  refdata = seu_fdc$annotation_20230508,  # or ref_seurat$celltype if that's what it's called
  dims = 1:30
)
VU40T.combined <- AddMetaData(object = VU40T.combined, metadata = predictions)

p1 <- DimPlot(VU40T.combined, group.by = "predicted.id", split.by = "condition" , label = TRUE, label.size = 3)
p2 <- DimPlot(VU40T.combined, group.by = "seurat_clusters", split.by = "condition")
if(all(VU40T.combined$species == "Mouse")){
  png(file = paste0(git_dir,"/Integrated/Mouse/Plots/Azimuth_annotated_UMAP_VU40T_combined.png"), width = 12, height = 10, units = "in", res = 300)
}else{
  png(file = paste0(git_dir,"/Integrated/Plots/Azimuth_annotated_UMAP_VU40T_combined.png"), width = 12, height = 10, units = "in", res = 300)
}
print(p1 + p2)
dev.off()

## for mouse
if (all(VU40T.combined$species == "Mouse")) {
  colnames(VU40T.combined) = sapply(colnames(VU40T.combined), stringr::str_to_title)
  ref <- MouseRNAseqData(cell.ont = "nonna")
}else{
  ref <- readRDS("~/R/HPCA_reference.rds") 
}


cluster_labels <- VU40T.combined$seurat_clusters
bp <- BiocParallel::MulticoreParam(workers = 8) 
# Run SingleR annotation
sce <- Seurat::as.SingleCellExperiment(VU40T.combined)

### doesn't work with human-aligned integrated. kinda works with mouse
## doesn't work with mouse only
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

png(file = file.path(git_dir, "Integrated/Plots/HumanOnly_10_PCs_SingleR_annotated_VU40T_combined_res_0.3_labelMain.png"), width = 10, height = 8, units = "in", res = 300)

# Visualize cell type annotation on UMAP
p <- Seurat::DimPlot(VU40T.combined, group.by='SingleR_HTA_cluster_label', split.by = "condition")
p2 <- DimPlot(VU40T.combined, reduction = "umap", split.by = "condition", group.by = "seurat_clusters", label = T) 
print(p + p2)
#print (p2)
dev.off()


### marker dotplots
genelist <- read.csv("~/Markers_for_dotplots_2_set1.csv", header = T) ## both human and mouse markers

make_marker_dotplots(VU40T.combined, genelist, export = T, outPrefix = paste0(Species, "_onlySET1"))

## genelist 2 -> human only
if (all(VU40T.combined$species == "Human")){
  genelist <- read.csv("~/markers_secondset_human_ep.csv", header = T)
  make_marker_dotplots(VU40T.combined, genelist, export = T, outPrefix = paste0(Species, "_onlySET2_Epithelial"))
  ### for heatmap
  heatmap_markers <- genelist[, c(1:5)]
  colnames(heatmap_markers) <- stringr::str_replace(colnames(heatmap_markers), "\\.", " ")
  colnames(heatmap_markers) <- stringr::str_replace(colnames(heatmap_markers), "junctions", "Junctions")
  # Assuming your marker df is named marker_df
  marker_long <- heatmap_markers %>%
    tidyr::pivot_longer(cols = everything(), names_to = "Group", values_to = "Gene") %>%
    dplyr::filter(!is.na(Gene) & Gene != "")
  # Filter for genes that are present
  marker_long <- marker_long %>% dplyr::filter(Gene %in% rownames(VU40T.combined))
  
  ### line wrapping func.
  marker_long$Group <- sapply(marker_long$Group, function(x) {
    if (nchar(x) > 15 && grepl(" ", x)) {
      gsub(" ", "\n", x)  # replace . or _ with a line break
    } else {
      x  # keep as is
    }
  })
  # Compute average expression
  avg_expr <- AverageExpression(VU40T.combined, features = unique(marker_long$Gene), return.seurat = FALSE, assays = "RNA")$RNA
  avg_expr@Dimnames[[2]] <- stringr::str_replace(avg_expr@Dimnames[[2]], "g", "")
  # Order genes as in the original df
  marker_long <- marker_long %>%
    distinct(Gene, .keep_all = TRUE)  # drop duplicate entries
  
  # Reorder avg_expr to match marker order
  avg_expr <- avg_expr[marker_long$Gene, ]
  
  scaled_expr <- t(scale(t(avg_expr)))
  # Clip values to range [-2, 2] for visual clarity
  scaled_expr[scaled_expr > 2] <- 2
  
  # Annotation (for row group labels)
  row_annot <- rowAnnotation(
    Group = marker_long$Group,
    col = list(Group = structure(
      circlize::rand_color(length(unique(marker_long$Group))),
      names = unique(marker_long$Group)
    )),
    show_annotation_name = TRUE,
    show_legend = FALSE 
  )
  
  # Select 3-color palette from RColorBrewer
  color_palette <- RColorBrewer::brewer.pal(n = 3, name = "YlOrRd")
  
  # Build color function
  col_fun <- colorRamp2(
    breaks = c(-2, 0, 2),
    colors = color_palette
  )
  png(file = file.path(
    git_dir,
    "Integrated/Plots/HeatMap_epiMarkers_clusterResolution0.6_VU40T_combined.png"),
    width = 8, height = 10, units = "in", res = 300)
  p <- Heatmap(
    matrix = as.matrix(scaled_expr),
    name = "Z-scaled Expr",
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10),
    column_names_rot = 0,
    row_split = marker_long$Group,
    left_annotation = row_annot,
    col = col_fun,
    heatmap_legend_param = list(title = "Z-scaled\nExpression", fontsize = 8),
    use_raster = TRUE,
  )
  draw(p,   padding = unit(c(10, 20, 10, 10), "mm"))  # prevent clipping
  dev.off()
  
  ### UMAP overlays in viridis
  
  umap_markers <- c("KRT6A",
                    "DSC3",
                    "CDH1",
                    "ITGB1",
                    "LAMB3",
                    "ITGAV",
                    "HAS2",
                    "TRIO",
                    "RHOA",
                    "KRT81",
                    "CLDN1",
                    "EPCAM",
                    "GCLM"
  )
  
  png(file = file.path(
    git_dir,
    "Integrated/Plots/Marker_Overlay_UMAPs_EpiMarkers_Viridis_clusterResolution0.3_VU40T_combined.png"),
    width = 14, height = 14, units = "in", res = 300)
  p <- FeaturePlot(VU40T.combined, 
                   features = umap_markers, 
                   label = F,# label.size = 3, repel = T,
  ) & 
    scale_color_viridis_c()
  print(p)
  dev.off()
  ### new human dotplots
  genelist <- readxl::read_excel(file.path(proj_dir,"scRNAseqAnalysis/Markers_for_dotplots_2.xlsx"), sheet = 8)
  make_marker_dotplots(VU40T.combined, genelist, export = T, outPrefix = paste0(Species, "_onlySET3_Epithelial"))
  
  ### and violin plots
  genelist <- readxl::read_excel(file.path(proj_dir,"scRNAseqAnalysis/Markers_for_dotplots_2.xlsx"), sheet = 9,col_names = F)
  make_marker_dotplots(VU40T.combined, genelist, export = T, outPrefix = paste0(Species, "_onlySET3_Epithelial"))
  
  
  for (col in seq_len(ncol(genelist))) {
    png(file = file.path(
      git_dir,
      paste0("Integrated/Plots/ViolinPlot_EpiMarkers_set3_clusterResolution0.3_VU40T_combined",col,".png")),
      width = 6, height = 10, units = "in", res = 300)
    p <- VlnPlot(
      VU40T.combined,
      features = na.omit(genelist[[col]]),
      alpha = 1, pt.size = 0,
      combine = FALSE
    ) 
    
    # remove legends + tighten margins
    p <- lapply(p, function(pp) {
      pp + theme(
        legend.position = "none",
        plot.margin = margin(2, 0.5, 2, 1),
        axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 10, angle = 0, vjust = 0.5),
      )
    })
    
    for (i in seq_along(p)){
      p[[i]][["labels"]][["x"]] <- ""
      p[[i]][["labels"]][["y"]] <- ""
    } 
    
    p_stack <- wrap_plots(p, ncol = 2, nrow = 5,  guides = "collect")
    
    print(p_stack)
    dev.off
    
    ### new human dotplots -> HNSCC and EMT
    genelist <- readxl::read_excel(file.path(proj_dir,"scRNAseqAnalysis/Markers_for_dotplots_3.xlsx"), sheet = 11, skip = 1)
    make_marker_dotplots(VU40T.combined, genelist, export = T, outPrefix = paste0(Species, "_Fig4_Epithelial"), plot_dir = plots_dir)
    
    FeaturePlot(VU40T.combined, features=c("DDR1", "RHOA")) & scale_color_viridis()
    
    pemt_genes <- c(
      "VIM","FN1","ITGA5","ITGB1","SERPINE1","LAMC2","LAMB3","PDPN",
      "MMP14","MMP1","MMP10","TGFBI","PLAUR","FSCN1","COL7A1",
      "CXCL8","CXCL1","SNAI2","ZEB1","KRT14","ITGA6"
    )
    rho_acto_genes <- c(
      "RHOA","RHOC","ROCK1","ROCK2","MYL9","MYH9","MYH10","ACTN1",
      "TAGLN","DIAPH1","DIAPH3","LIMK1","LIMK2","CFL1","PFN1",
      "PPP1R12A","MYLK","DAAM1","ARHGEF11","ARHGEF12"
    )
    modules = list(pemt_genes, rho_acto_genes)
    names(modules) <- c("pemt_genes", "rho_acto_genes")
    seu <- VU40T.combined
    for (nm in names(modules)) {
      
      seu <- AddModuleScore(seu, features = list(modules[[nm]]), name = nm, nbin = 24, ctrl = 100)
      
    }
    png(file.path(plots_dir, "rho_acto_pemt_module_umap.png"), 
        width = 12, height = 8, units = "in", res = 300)
    
    
    p <- FeaturePlot(seu, features = c("pemt_genes1","rho_acto_genes1"),
                     order = TRUE, cols = viridisLite::viridis(n = 5))
    print(p)
    dev.off()
    
    png(file.path(plots_dir, "rho_acto_pemt_module_violin.png"), 
        width = 12, height = 8, units = "in", res = 300)
    
    p <- VlnPlot(seu, features = c("pemt_genes1","rho_acto_genes1"))
    print(p)
    dev.off()
  }
  table(VU40T.combined$seurat_clusters)
  
  clust6 <- subset(seu, idents = 6)
  clust6$pemt_genes1
  clust6$rho_acto_genes1
  
  # Create dataframe
  clust6_df <- data.frame(
    cell_barcodes = names(clust6$pemt_genes1),
    pemt_genes = clust6$pemt_genes1,
    rho_acto_genes = clust6$rho_acto_genes1
  )
  
  # Define color and alpha rules
  clust6_df <- clust6_df %>%
    mutate(
      label_color = case_when(
        pemt_genes > 0.8 ~ "red",
        rho_acto_genes > 0.6 ~ "blue",
        TRUE ~ "grey70"
      ),
      alpha_val = case_when(
        pemt_genes > 0.8 ~ 1,
        rho_acto_genes > 0.6 ~ 1,
        TRUE ~ 0.2
      )
    )
  
  # Pivot to long format for plotting
  clust6_long <- clust6_df %>%
    tidyr::pivot_longer(
      cols = c(pemt_genes, rho_acto_genes),
      names_to = "module",
      values_to = "expression"
    )
  
  
  # Plot
  png(file.path(plots_dir, "cluster6_rho_acto_pemt_module_scatter_plot.png"), 
      width = 12, height = 8, units = "in", res = 300)
  
  p <- ggplot(clust6_long, aes(x = cell_barcodes, y = expression, color = module, alpha = alpha_val)) +
    geom_point() +
    scale_alpha_identity() +
    guides(alpha = "none") +
    labs(x = "Cell barcodes", y = "Expression", color = "Module") +
    scale_x_discrete(
      labels = setNames(
        paste0("<span style='color:", clust6_df$label_color, "'>", clust6_df$cell_barcodes, "</span>"),
        clust6_df$cell_barcodes
      )
    ) + theme_minimal() +
    theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1))
  
  print(p)
  dev.off()
}
## genelist 3 -> mouse and human epithelial, fibroblast and EMT markers
genelist <- as.data.frame(read.csv("~/Markers_for_dotplots_2_ep_2nd_set.csv", header = T, skip = 1))
genelist <- as.data.frame(lapply(genelist, toupper))
genelist[] <- lapply(genelist, function(x) gsub("CDH2", "", x))

if (all(VU40T.combined$species == "Mouse")){
  genelist <- genelist[, c(4,5)] ## both human and mouse markers w/out emt
  colnames(genelist) <- gsub("_[HM]", "", colnames(genelist))
  
}else{
  genelist <- genelist[, c(1:2)] ## both human and mouse markers, no need for EMT here
  colnames(genelist) <- gsub("_[HM]", "", colnames(genelist))
}
colnames(genelist) <- gsub("_", " ", colnames(genelist))
colnames(genelist) <- gsub("markers", "Markers", colnames(genelist))
## convert to named vectors
genevectors <- lapply(genelist, function(x) unique(x[x != "" & !is.na(x)]))

if (all(VU40T.combined$species == "Human")){
  png(file = file.path(
    git_dir,
    "Integrated/Plots/CombinedDotplot_epi_fibro_Markers_clusterResolution0.3_VU40T_combined.png"),
    width = 8, height = 4, units = "in", res = 300)
  
} else {
  png(file = file.path(
    git_dir,
    "Integrated/Mouse/Plots/CombinedDotplot_epi_Fibro_Markers_clusterResolution0.6_VU40T_combined.png"),
    width = 8, height = 4, units = "in", res = 300)
}

p <- DotPlot(
  VU40T.combined,
  features = genevectors,
  group.by = "seurat_clusters",
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic")  # Optional: italics for gene names
  ) + labs(y = "Clusters")  # 👈 Relabel y-axis
print(p)  
dev.off()

if (all(VU40T.combined$species == "Mouse")){
  VU40T.combined <- subset(VU40T.combined, idents = "13", invert = TRUE)
  Species = "Mouse"
  
  ## make stacked barplot of num fibroblasts per clust
  clust_cell_df <- as.data.frame(table(VU40T.combined$seurat_clusters))
  colnames(clust_cell_df) <- c("Cluster", "No. Cells")
  write.csv(clust_cell_df, file = file.path(git_dir, "Integrated/Mouse/MouseOnly_fibroblast_cells_per_clust.csv"))
  clust_cell_df$`Cell Type` <- "Fibroblast"
  ## plot barplot
  
  png(file = file.path(
    git_dir,
    "Integrated/Mouse/Plots/MouseONLY_Stacked_Bar_numCells_per_cluster_Resolution0.6_VU40T_combined.png"),
    width = 4, height = 6, units = "in", res = 300)
  
  p <- clust_cell_df %>% 
    ggplot(aes(x = `Cell Type`, y = `No. Cells`, fill = Cluster)) + 
    geom_bar(position="stack", stat="identity") + 
    theme_minimal() + xlab("")
  print(p)
  dev.off()
  
  png(file = file.path(
    git_dir,
    "Integrated/Mouse/Plots/MouseONLY_Stacked_Bar_percent_Cells_per_cluster_Resolution0.6_VU40T_combined.png"),
    width = 4, height = 6, units = "in", res = 300)
  p <- clust_cell_df %>% 
    ggplot(aes(x = `Cell Type`, y = `No. Cells`, fill = Cluster)) + 
    geom_bar(position="fill", stat="identity") + 
    theme_minimal() +
    ylab("Proportion of Cells") + xlab("") + 
    scale_y_continuous(labels = scales::percent)
  print(p)
  dev.off()
  genelist <- read.csv("~/Fibroblasts_dotplots_markers.csv", header = F)
  names(genelist) <- "Fibroblast Markers"
  
  png(file = file.path(
    git_dir,
    "Integrated/Mouse/Plots/MouseONLY_MarkerDotplots/Dotplot_Fibro_Markers_clusterResolution0.6_VU40T_combined.png"),
    width = 12, height = 4, units = "in", res = 300)
  
  p <- DotPlot(
    VU40T.combined,
    features = genelist,
    group.by = "seurat_clusters",
  ) +
    scale_x_discrete(labels = function(x) stringr::str_to_title(x)) +
    RotatedAxis() +
    scale_color_gradient(low = "lightgrey", high = "red") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(face = "italic")  # Optional: italics for gene names
    ) + labs(y = "Clusters")  # 👈 Relabel y-axis
  print(p)  
  dev.off()
  
  ## fibroblast heatmap
  # heatmap_markers <- read.csv("~/fibroblast_heatmap.csv", header = T)
  fibro2_marker_ls <- list()
  for (i in seq(1,6)){
    fibro2_marker_ls[[i]] <- readxl::read_xlsx(file.path(proj_dir ,"/scRNAseqAnalysis/Fibroblasts_dotplots.xlsx"), sheet = i)
  }
  heatmap_markers <- fibro2_marker_ls[[1]]
  ### for heatmap
  # heatmap_markers <- genelist[, c(1:5)]
  colnames(heatmap_markers) <- stringr::str_replace(colnames(heatmap_markers), "_", " ")
  ## for fibro list 2
  colnames(heatmap_markers) <- stringr::str_replace(colnames(heatmap_markers), "Oxid phosph", "Oxidative phosphorylation")
  colnames(heatmap_markers) <- stringr::str_to_title(colnames(heatmap_markers))
  colnames(heatmap_markers) <- stringr::str_replace(colnames(heatmap_markers), "Ecm", "ECM")
  colnames(heatmap_markers) <- stringr::str_replace(colnames(heatmap_markers), "To", "to")
  marker_long <- heatmap_markers %>%
    tidyr::pivot_longer(cols = everything(), names_to = "Group", values_to = "Gene") %>%
    dplyr::filter(!is.na(Gene) & Gene != "")
  # Filter for genes that are present
  marker_long <- marker_long %>% dplyr::filter(Gene %in% rownames(VU40T.combined))
  
  ### line wrapping func.
  marker_long$Group <- sapply(marker_long$Group, function(x) {
    if (nchar(x) > 15 && grepl(" ", x)) {
      gsub(" ", "\n", x)  # replace . or _ with a line break
    } else {
      x  # keep as is
    }
  })
  # Compute average expression
  avg_expr <- AverageExpression(VU40T.combined, features = unique(marker_long$Gene), return.seurat = FALSE, assays = "RNA")$RNA
  avg_expr@Dimnames[[2]] <- stringr::str_replace(avg_expr@Dimnames[[2]], "g", "")
  # Order genes as in the original df
  marker_long <- marker_long %>%
    distinct(Gene, .keep_all = TRUE)  # drop duplicate entries
  
  # Reorder avg_expr to match marker order
  avg_expr <- avg_expr[marker_long$Gene, ]
  
  scaled_expr <- t(scale(t(avg_expr)))
  # Clip values to range [-2, 2] for visual clarity
  scaled_expr[scaled_expr > 2] <- 2
  
  grp <- as.character(marker_long$Group)
  lev <- sort(unique(grp))
  set.seed(123)
  pal <- setNames(
    circlize::rand_color(length(lev), distinct = TRUE),
    lev
  )
  
  
  row_annot <- rowAnnotation(
    Group = grp,                        
    annotation_legend_param = list(
      Group = list(title = "Group") # legend title
    ))
  
  # # Annotation (for row group labels)
  
  # Select 3-color palette from RColorBrewer
  color_palette <- RColorBrewer::brewer.pal(n = 3, name = "YlOrRd")
  
  # Build color function
  col_fun <- colorRamp2(
    breaks = c(-2, 0, 2),
    colors = color_palette
  )
  png(file = file.path(
    git_dir,
    "Integrated/Mouse/Plots/Heatmap_FibroMarkers_set2_clusterResolution0.6_VU40T_combined.png"),
    width = 8, height = 14, units = "in", res = 300)
  p <- Heatmap(
    matrix = as.matrix(scaled_expr),
    name = "Z-scaled Expr",
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_names = T,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10),
    column_names_rot = 0,
    row_split = marker_long$Group,      # still split by group
    row_title = NULL,                   # <- hide group slice labels
    row_title_gp = gpar(fontsize = 0),  # <- belt-and-braces: ensure invisible
    left_annotation = row_annot,
    col = col_fun,
    heatmap_legend_param = list(title = "Z-scaled\nExpression", fontsize = 8),
    use_raster = TRUE,
  )
  ## because it's mouse let's modify gene names so they fit convention: sentence case
  p@row_names_param[["labels"]] <- stringr::str_to_sentence(p@row_names_param[["labels"]])
  p@row_names_param[["anno"]]@var_env[["value"]] <- stringr::str_to_sentence(p@row_names_param[["anno"]]@var_env[["value"]])
  
  draw(p,   padding = unit(c(10, 20, 10, 10), "mm"))  # prevent clipping
  dev.off()
  
  ### violinplot
  
  ViolinGenes <- fibro2_marker_ls[[4]]
  # make sure header gets included in vector
  colname <- colnames(ViolinGenes)[1]
  ViolinGenes <- c(colname, as.character(ViolinGenes[[1]]))
  ## add extra markers
  ViolinGenes <- append(ViolinGenes, 
                        values = c( # "ACTA2", <- wasn't that informative, all show minimal expression
                          "P4HA1", "CTSD")
  )
  png(file = file.path(
    git_dir,
    "Integrated/Mouse/Plots/ViolinPlot_FibroMarkers_set2_clusterResolution0.6_VU40T_combined.png"),
    width = 6, height = 10, units = "in", res = 300)
  p <- VlnPlot(
    VU40T.combined,
    features = ViolinGenes,
    alpha = 1, pt.size = 0,
    combine = FALSE
  ) 
  
  # remove legends + tighten margins
  p <- lapply(p, function(pp) {
    pp + theme(
      legend.position = "none",
      plot.margin = margin(2, 0.5, 2, 1),
      axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 10, angle = 0, vjust = 0.5),
    )
  })
  
  for (i in seq_along(p)){
    p[[i]][["labels"]][["title"]] <- stringr::str_to_sentence(p[[i]][["labels"]][["title"]])
    p[[i]][["labels"]][["x"]] <- ""
    p[[i]][["labels"]][["y"]] <- ""
  } 
  
  p_stack <- wrap_plots(p, ncol = 2, nrow = 5,  guides = "collect")
  
  print(p_stack)
  dev.off()
  
  ### v2 violinplots
  
  
  ViolinGenes <- as.data.frame(fibro2_marker_ls[[5]])
  
  # desired height per subplot in inches
  per_gene_height <- 1   
  # 1) get a clean character vector of genes for this panel/column
  for (panel in colnames(ViolinGenes)) {
    # 1) get a clean character vector of genes for this panel/column
    genes <- ViolinGenes[[panel]]
    genes <- unique(na.omit(as.character(genes)))
    # keep only genes present in the object
    genes <- intersect(genes, rownames(VU40T.combined))
    if (length(genes) == 0) next
    
    plots <- lapply(seq_along(genes), function(i) {
      g <- genes[i]
      gp <- VlnPlot(
        VU40T.combined,
        features = g,
        group.by = "seurat_clusters",
        pt.size  = 0,
        combine  = TRUE
      ) + 
        NoLegend() +
        labs(title = "", x = "",
             y = stringr::str_to_sentence(tolower(g))) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.6, size = 14),
          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14),
          axis.text.y = element_text(size = 8, angle = 0),
          plot.margin = margin(2, 5, 2, 5)
        )
      
      # Remove x-axis elements for all but the last plot
      if (i < length(genes)) {
        gp <- gp +
          theme(
            axis.title.x = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank()
          )
      } else {
        gp <- gp +
          labs(x = "") +
          theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
      }
      gp
    })
    
    # 3) stack vertically and add a global title
    panel_chr <- as.character(panel)
    p_stack <- wrap_plots(plots, ncol = 1, guides = "collect") +
      plot_annotation(
        title = panel_chr,
        theme = theme(plot.title = element_text(hjust = 0.6, size = 16, face = "bold"))
      )
    
    # output height = per-gene height * number of genes + extra for title
    h_in <- per_gene_height * length(plots) + 1
    # 4) save with png()/dev.off()
    png(
      file = file.path(
        git_dir,
        paste0(
          "Integrated/Mouse/Plots/V2_ViolinPlot_FibroMarkers_",
          gsub("[^A-Za-z0-9._-]+", "_", panel_chr),
          "_clusterResolution0.6_VU40T_combined.png"
        )
      ),
      width = 6, height = h_in, units = "in", res = 300
    )
    print(p_stack)
    dev.off()
  }
  
  
  
  ### overlay umaps
  
  umap_markers <- fibro2_marker_ls[[6]]
  colname <- colnames(umap_markers)[1]
  umap_markers <- c(colname, as.character(umap_markers[[1]]))
  umap_markers<- append(umap_markers, "ZBTB7B")
  
  png(file = file.path(
    git_dir,
    "Integrated/Mouse/Plots/Marker_Overlay_UMAPs_FibroMarkers_set2_clusterResolution0.6_VU40T_combined.png"),
    width = 21, height = 21, units = "in", res = 300)
  p <- FeaturePlot(VU40T.combined, 
                   features = umap_markers, 
                   label = F,# label.size = 3, repel = T,
  ) & 
    scale_color_viridis_c()
  for (i in seq_along(p)){
    p[[i]][["labels"]][["title"]] <- stringr::str_to_sentence(p[[i]][["labels"]][["title"]])
  } 
  print(p)
  dev.off()
}
#test


### CSC EMT_geneList

if (all(VU40T.combined$species == "Human")){
  
  genelist <- read.csv(file.path(proj_dir, "scRNAseqAnalysis/CSC_EMT_markers.csv"), header = T)
  # genelist[] <- lapply(genelist, function(x) gsub("CDH2", "", x))
  colnames(genelist) <- gsub("_", " ", colnames(genelist))
  ## convert to named vectors
  genevectors <- lapply(genelist, function(x) unique(x[x != "" & !is.na(x)]))
  
  png(file = file.path(
    git_dir,
    "Integrated/Plots/CombinedDotplot_CSC_EMT_Markers_clusterResolution0.3_VU40T_combined.png"),
    width = 8, height = 4, units = "in", res = 300)
  
  p <- DotPlot(
    VU40T.combined,
    features = genevectors,
    group.by = "seurat_clusters",
  ) +
    RotatedAxis() +
    scale_color_gradient(low = "lightgrey", high = "red") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(face = "italic")  # Optional: italics for gene names
    ) + labs(y = "Clusters")  # 👈 Relabel y-axis
  print(p)  
  dev.off()
  
  
  ViolinGENEs <- c("CD44", "MET", "ALDH1A3", "COL17A1", "CDKN1A", "FOS", "JUN", "ID1", "ATF3", "BTG2", "DUSP1")
  
  png(file = file.path(
    git_dir,
    "Integrated/Plots/ViolinPlot_CSC_EPI_Markers_clusterResolution0.3_VU40T_combined.png"),
    width = 12, height = 8, units = "in", res = 300)
  
  p <- VlnPlot(
    VU40T.combined,
    features = ViolinGENEs,
    alpha = 1, pt.size = 0
  )
  
  print(p)
  dev.off()
  
  basalGENES <-  c("COL17A1", "CDKN1A", "FOS", "JUN", "ID1", "ATF3", "BTG2", "DUSP1")
  
  png(file = file.path(
    git_dir,
    "Integrated/Plots/ViolinPlot_Basal_Markers_clusterResolution0.3_VU40T_combined.png"),
    width = 6, height = 8, units = "in", res = 300)
  
  p <- VlnPlot(
    VU40T.combined,
    features = basalGENES,
    alpha = 1, pt.size = 0,
    idents = c(0,1,2,5)
  )
  
  print(p)
  dev.off()
}
for (i in seq_len(ncol(genelist))) {
  gene_vector <- na.omit(genelist[[i]])  # extract column as vector and remove NAs
  
  if (length(gene_vector) > 0) {  # skip empty gene sets
    p <- FeaturePlot(VU40T.combined, features = gene_vector, label = T)
    
    print(p)  # or save if needed
  }
}

for (i in genelist$Fibroblast_markers){
  # Choose your gene of interest
  gene <- i
  
  print(gene)
  # Get expression matrix (log-normalized by default)
  expr <- GetAssayData(VU40T.combined, slot = "data", assay = "RNA")[gene, ]
  
  expressing_cells <- expr > 0
  
  # Get cluster labels for those cells
  cluster_labels <- Idents(VU40T.combined)
  
  # Total cells per cluster
  total_cells <- table(cluster_labels)
  
  # Find cells where expression > 0
  expressing_cells <- expr > 0
  
  # Tabulate by cluster
  print(table(Idents(VU40T.combined)[expressing_cells]))
  
  # Expressing cells per cluster
  expressing_cells_per_cluster <- table(cluster_labels[expressing_cells])
  
  # Ensure all clusters are represented (even if 0 cells express the gene)
  all_clusters <- names(total_cells)
  percent_expressing <- sapply(all_clusters, function(clust) {
    total <- total_cells[clust]
    expressed <- expressing_cells_per_cluster[clust]
    if (is.na(expressed)) expressed <- 0
    round((expressed / total) * 100, 2)  # percentage rounded to 2 decimal places
  })
  
  # Convert to a named vector or data frame
  percent_expressing_df <- data.frame(
    Cluster = names(percent_expressing),
    PercentExpressing = percent_expressing
  )
  
  # View the result
  print(percent_expressing_df)
} 


for (i in genelist$Epithelial_markers){
  # Choose your gene of interest
  gene <- i
  
  print(gene)
  # Get expression matrix (log-normalized by default)
  expr <- GetAssayData(VU40T.combined, slot = "data", assay = "RNA")[gene, ]
  
  expressing_cells <- expr > 0
  
  # Get cluster labels for those cells
  cluster_labels <- Idents(VU40T.combined)
  
  # Total cells per cluster
  total_cells <- table(cluster_labels)
  
  # Find cells where expression > 0
  expressing_cells <- expr > 0
  
  # Tabulate by cluster
  print(table(Idents(VU40T.combined)[expressing_cells]))
  
  # Expressing cells per cluster
  expressing_cells_per_cluster <- table(cluster_labels[expressing_cells])
  
  # Ensure all clusters are represented (even if 0 cells express the gene)
  all_clusters <- names(total_cells)
  percent_expressing <- sapply(all_clusters, function(clust) {
    total <- total_cells[clust]
    expressed <- expressing_cells_per_cluster[clust]
    if (is.na(expressed)) expressed <- 0
    round((expressed / total) * 100, 2)  # percentage rounded to 2 decimal places
  })
  
  # Convert to a named vector or data frame
  percent_expressing_df <- data.frame(
    Cluster = names(percent_expressing),
    PercentExpressing = percent_expressing
  )
  
  # View the result
  print(percent_expressing_df)
} 

gene <- "RHOA"

print(gene)
# Get expression matrix (log-normalized by default)
expr <- GetAssayData(VU40T.combined, slot = "data", assay = "RNA")[gene, ]

expressing_cells <- expr > 0

# Get cluster labels for those cells
cluster_labels <- Idents(VU40T.combined)

# Total cells per cluster
total_cells <- table(cluster_labels)

# Find cells where expression > 0
expressing_cells <- expr > 0

# Tabulate by cluster
print(table(Idents(VU40T.combined)[expressing_cells]))

# Expressing cells per cluster
expressing_cells_per_cluster <- table(cluster_labels[expressing_cells])

# Ensure all clusters are represented (even if 0 cells express the gene)
all_clusters <- names(total_cells)
percent_expressing <- sapply(all_clusters, function(clust) {
  total <- total_cells[clust]
  expressed <- expressing_cells_per_cluster[clust]
  if (is.na(expressed)) expressed <- 0
  round((expressed / total) * 100, 2)  # percentage rounded to 2 decimal places
})

# Convert to a named vector or data frame
percent_expressing_df <- data.frame(
  Cluster = lapply(names(percent_expressing), integer),
  PercentExpressing = percent_expressing
)

# View the result
print(percent_expressing_df)

FeaturePlot(VU40T.combined, features = "RHOA", label = T)

make_marker_dotplots(VU40T.combined, genelist, export = T, outPrefix = paste0(Species, "_epi_fibro"))





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
        git_dir, paste0("/Integrated/Mouse/Plots/MouseOnly_GSEA_reactome_cluster_", clust, ".png")
      ),
      width = 8, height = 15, units = "in", res = 300
    )
    print(
      dotplot(gsea_result, showCategory = 20) + 
        ggtitle(paste("Cluster", clust, "GSEA REACTOME"))
    )
    dev.off()
  } else {
    message("No significant pathways for Cluster ", clust, " — skipping plot.")
  }
}


## human (can do bp too)
go_human_RM <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol)

# Needed format for clusterProfiler GSEA:
go_human_RM_list <- split(go_human_RM$gene_symbol, go_human_RM$gs_name)

# For each cluster, get top 200 genes by avg_log2FC
top_genes_per_cluster <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 200)


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
        git_dir, paste0("/Integrated/Plots/Top200_HumanOnlyGSEA_reactome_cluster_", clust, ".png")
      ),
      width = 5, height = 10, units = "in", res = 300
    )
    print(
      dotplot(gsea_result, showCategory = 20) + 
        ggtitle(paste("Cluster", clust, "GSEA REACTOME"))
    )
    dev.off()
  } else {
    message("No significant pathways for Cluster ", clust, " — skipping plot.")
  }
}


# Function to rename rownames across an assay
make_gene_names_uppercase <- function(seurat_obj, assay = "RNA") {
  # Get the current assay data
  assay_data <- seurat_obj[[assay]]
  
  # Create mapping from old to new names
  old_names <- rownames(assay_data@counts)
  new_names <- toupper(old_names)
  
  # Apply to counts, data, and scale.data
  rownames(assay_data@counts) <- new_names
  rownames(assay_data@data) <- new_names
  rownames(assay_data@scale.data) <- new_names
  
  # Apply to meta.features
  rownames(assay_data@meta.features) <- new_names
  
  # Replace back into the Seurat object
  seurat_obj[[assay]] <- assay_data
  
  return(seurat_obj)
}

# Apply to your object

if (all(VU40T.combined$species == "Mouse")){
  VU40T.combined <- make_gene_names_uppercase(VU40T.combined)
}

## signed averages

# Step 1: Define gene sets (convert to uppercase)
epithelial_genes <- toupper(c("Cdh1", "Cldn7", "Krt17", "Krt19", "Dsp"))
fibroblast_genes <- toupper(c("Vim", "Col1a1", "S100a4", "Pdgfra", "Fn1"))
emt_genes <- toupper(c("Acta2", "Cdh2", "Snai1", "Twist2", "Zeb1"))

# Step 2: Set correct assay and identity class
DefaultAssay(VU40T.combined) <- "RNA"
Idents(VU40T.combined) <- "seurat_clusters"

# Step 3: Compute average expression per cluster for selected genes
cluster_avgs <- AverageExpression(
  VU40T.combined, 
  features = unique(c(epithelial_genes, fibroblast_genes, emt_genes)), 
  slot = "data",
  assay = "RNA"
)$RNA

# Step 4: Compute mean scores for each gene set per cluster
epithelial_mean <- colMeans(cluster_avgs[intersect(epithelial_genes, rownames(cluster_avgs)), , drop = FALSE])
fibroblast_mean <- colMeans(cluster_avgs[intersect(fibroblast_genes, rownames(cluster_avgs)), , drop = FALSE])
emt_mean <- colMeans(cluster_avgs[intersect(emt_genes, rownames(cluster_avgs)), , drop = FALSE])

# Step 5: Combine into a data frame
cluster_scores <- data.frame(
  Cluster = names(epithelial_mean),
  Epithelial = epithelial_mean,
  Fibroblast = fibroblast_mean,
  EMT = emt_mean
)

# Step 6: Assign annotation based on dominant expression
cluster_scores$Phenotype <- apply(cluster_scores[, c("Epithelial", "Fibroblast", "EMT")], 1, function(x) {
  types <- c("Epithelial", "Fibroblast", "EMT")
  types[which.max(x)]
})

# Step 7: Add annotations back to Seurat metadata
cluster_annotations <- cluster_scores$Phenotype
names(cluster_annotations) <- cluster_scores$Cluster
# Make sure seurat_clusters is a factor or character
cluster_ids <- as.character(VU40T.combined$seurat_clusters)

# Fix the cluster names in the score table
cluster_scores$Cluster <- gsub("^g", "", cluster_scores$Cluster)

# Now map phenotype by cluster ID per cell
cluster_ids <- as.character(VU40T.combined$seurat_clusters)
VU40T.combined$Phenotype <- cluster_scores$Phenotype[match(cluster_ids, cluster_scores$Cluster)]


# Step 8: Plot UMAP with phenotype labels
if (all(VU40T.combined$species == "Human")){
  png(file.path(git_dir, "Integrated/Plots/HumanOnly_UMAP_cluster_annotations.png")
      , width = 8, height = 6, units = "in", res = 300)
}else{
  png(file.path(git_dir, "Integrated/Mouse/Plots/MouseOnly_UMAP_cluster_annotations.png"), 
      width = 8, height = 6, units = "in", res = 300)
}

p <- DimPlot(VU40T.combined, group.by = "Phenotype", label = TRUE, repel = TRUE) +
  ggtitle("UMAP Colored by Dominant Gene Set Phenotype")
print(p)
dev.off()

# Optional: View cluster scores table
print(cluster_scores)

write.csv(cluster_scores, file = file.path(git_dir, "Integrated/Epi_fib_EMT_cluster_signed_avs.csv"))

# Define gene sign vector: +1 for epithelial, -1 for fibroblast and EMT
gene_weights <- c(
  setNames(rep(1, length(epithelial_genes)), epithelial_genes),
  setNames(rep(-1, length(fibroblast_genes)), fibroblast_genes),
  setNames(rep(-1, length(emt_genes)), emt_genes)
)

# Filter genes present in the Seurat object
valid_genes <- intersect(names(gene_weights), rownames(VU40T.combined[["RNA"]]))

# Get scaled or normalized data (slot = "data")
expr_mat <- GetAssayData(VU40T.combined, assay = "RNA", slot = "data")[valid_genes, ]

# Apply gene weights
signed_expr <- expr_mat * gene_weights[valid_genes]

# Compute per-cell signed score
signed_score <- Matrix::colMeans(signed_expr)

# Add to metadata
VU40T.combined$signed_score <- signed_score

cluster_means <- tapply(VU40T.combined$signed_score, VU40T.combined$seurat_clusters, mean)

# Step 6: Annotate dominant phenotype per cluster
# Thresholds: >0 = Epithelial; <0 = Mesenchymal-like (Fibroblast/EMT);  ~0 = Mixed
cluster_labels <- ifelse(cluster_means > 0.25, "Epithelial",
                         ifelse(cluster_means < -0.25, "Mesenchymal", "Mixed"))

# Store cluster annotations in metadata
VU40T.combined$Phenotype <- cluster_labels[as.character(VU40T.combined$seurat_clusters)]

# Step 7: Prepare UMAP + overlay phenotype labels
# Extract UMAP coordinates
umap_df <- Embeddings(VU40T.combined, "umap") %>%
  as.data.frame() %>%
  `colnames<-`(c("UMAP_1", "UMAP_2")) %>%
  mutate(
    cluster = VU40T.combined$seurat_clusters,
    phenotype = VU40T.combined$Phenotype,
    score = VU40T.combined$signed_score  # <- use correct column
  )

# Get median UMAP coordinates for each cluster for labeling
label_df <- umap_df %>%
  group_by(cluster, phenotype) %>%
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")


p <- FeaturePlot(
  object = VU40T.combined,
  features = "signed_score",  # Use your actual score name
  label = FALSE,
  pt.size = 0.5
) +
  scale_color_gradient2(
    low = "blue",
    mid = "gray90",
    high = "red",
    midpoint = 0
  ) +
  ggtitle("Signed Score with Cluster Phenotype Labels") +
  theme(plot.title = element_text(hjust = 0.5))

# Step 2: Extract UMAP coordinates for labeling
umap_coords <- Embeddings(VU40T.combined, "umap")
umap_df <- as.data.frame(umap_coords)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cluster <- VU40T.combined$seurat_clusters
umap_df$phenotype <- VU40T.combined$Phenotype

# Step 3: Calculate cluster centers for labels
label_df <- umap_df %>%
  group_by(cluster, phenotype) %>%
  summarize(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )

# Step 4: Add phenotype labels to FeaturePlot
p + geom_text(
  data = label_df,
  aes(x = UMAP_1, y = UMAP_2, label = phenotype),
  size = 2.5,       # smaller than default 4
  fontface = "bold",
  color = "black"
)
ggsave(file.path(git_dir, "Integrated/Plots/HumanOnly_epi_signed_score_with_labels.png"), width = 8, height = 6, dpi = 300)


###pseudobulk

## isolate fibro clusters - all but 13

if(Species == "Mouse"){
  VU40T.combined <- subset(VU40T.combined, subset = seurat_clusters != 13)
}

library(edgeR)

#Get cluster IDs
cluster_ids <- unique(VU40T.combined$seurat_clusters)

pseudobulk_DE_results <- list()

# Loop over each cluster
for (clust in cluster_ids) {
  # Subset to cells in this cluster
  cells_in_cluster <- WhichCells(VU40T.combined, idents = clust)
  cluster_data <- subset(VU40T.combined, cells = cells_in_cluster)
  
  # Get raw counts and metadata
  counts <- GetAssayData(cluster_data, slot = "counts")
  meta <- cluster_data@meta.data
  
  # Aggregate counts by sample
  sample_ids <- unique(meta$sample)
  pb_counts <- sapply(sample_ids, function(sid) {
    cells <- rownames(meta)[meta$sample == sid]
    if (length(cells) == 1) {
      counts[, cells]
    } else {
      Matrix::rowSums(counts[, cells])
    }
  })
  
  # Set condition per sample
  condition <- meta %>%
    distinct(sample, condition) %>%
    filter(sample %in% colnames(pb_counts)) %>%
    arrange(match(sample, colnames(pb_counts))) %>%
    pull(condition)
  
  # edgeR pseudobulk
  group <- factor(condition)
  dge <- DGEList(counts = pb_counts, group = group)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~ group)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit)
  top <- topTags(qlf, n = Inf)
  
  # Save results
  pseudobulk_DE_results[[as.character(clust)]] <- top$table
}
for (i in seq_along(pseudobulk_DE_results)){
  if (Species == "Mouse"){
    write.csv(pseudobulk_DE_results[[i]], 
              file = file.path(git_dir, 
                               paste0("/Integrated/Mouse/Pseudobulk_DEs_MouseOnly_cluster", 
                                      names(pseudobulk_DE_results)[i], 
                                      ".csv" )
              )
    )
    
  }else{
    write.csv(pseudobulk_DE_results[[i]], 
              file = file.path(git_dir, 
                               paste0("/Integrated/Pseudobulk_DEs_HumanOnly_cluster", 
                                      names(pseudobulk_DE_results)[i], 
                                      ".csv" )
              )
    )
  }
  
}

library(EnhancedVolcano)

# Directory to save volcano plots
# out_dir <- "volcano_plots"
# dir.create(out_dir, showWarnings = FALSE)

# Loop through all clusters in your DE results
for (clust in names(pseudobulk_DE_results)) {
  de_table <- pseudobulk_DE_results[[clust]]
  
  # Skip if the table is empty
  if (is.null(de_table) || nrow(de_table) == 0) next
  
  # Build filename
  if (Species == "Mouse"){
    out_file <- file.path(git_dir, paste0("Integrated/Mouse/Plots/MouseOnly_PseudoBulk_volcano_Cluster_", clust, ".png"))
  }else{
    out_file <- file.path(git_dir, paste0("Integrated/Plots/HumanOnly_PseudoBulk_volcano_Cluster_", clust, ".png"))
  }
  
  message(paste("Generating volcano for cluster", clust, "→", out_file))
  # Create the pic
  png(out_file, width = 8, height = 6, units = "in", res = 300)
  
  select_labels <- head(rownames(de_table[order(de_table$FDR), ]), 10)
  select_labels <- select_labels[!is.na(select_labels) & select_labels != ""]
  
  print(EnhancedVolcano(
    de_table,
    lab = rownames(de_table),
    x = 'logFC',
    y = 'FDR',
    pCutoff = 0.1,
    FCcutoff = 1,
    title = paste("Epithelial Cluster", clust, ": LPS-P vs LPS-N"),
    pointSize = 2.0,
    labSize = 3,
    selectLab = select_labels,
    legendLabels = c('NS', 'Log2 FC', 'FDR', 'FDR & Log2 FC'),
    legendPosition = 'right',
    drawConnectors = TRUE,
    widthConnectors = 0.5
  ))
  
  dev.off()
}


## GSEA

if (Species == "Mouse"){
  # Mouse Reactome gene sets from MSigDB
  msig_mouse <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "REACTOME")
  
  # Create list: pathway name → gene vector
  reactome_t2g <- msig_mouse[, c("gs_name", "gene_symbol")]   # <- data.frame with 2 columns
  
  pseudobulk_DE_results <- lapply(pseudobulk_DE_results, function(df) {
    if (!is.null(df) && nrow(df) > 0) {
      rownames(df) <- stringr::str_to_title(rownames(df))
    }
    return(df)
  })
}else{
  # Create list: pathway name → gene vector
  reactome_t2g <- go_human_RM[, c("gs_name", "gene_symbol")]   # <- data.frame with 2 columns
}



gsea_results <- list()

for (clust in names(pseudobulk_DE_results)) {
  de_table <- pseudobulk_DE_results[[clust]]
  
  if (is.null(de_table) || nrow(de_table) == 0 || !"logFC" %in% colnames(de_table)) next
  
  # Prepare ranked gene list
  gene_list <- de_table$logFC
  names(gene_list) <- rownames(de_table)
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_list <- gene_list[!is.na(gene_list)]
  
  # Run GSEA using gene symbols + Reactome gene sets
  gsea <- tryCatch({
    GSEA(geneList = gene_list,
         TERM2GENE = reactome_t2g,
         pvalueCutoff = 0.05,
         verbose = FALSE)
  }, error = function(e) NULL)
  
  gsea_results[[clust]] <- gsea
}

library(enrichplot)

if (Species == "Mouse"){
  output_dir <- file.path(git_dir, "Integrated/Mouse/Plots/PseudoBulk_GSEA_Dotplots")
}else{
  output_dir <- file.path(git_dir, "Integrated/Plots/PseudoBulk_GSEA_Dotplots")
}
if (!(dir.exists(output_dir))){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}


for (clust in names(gsea_results)) {
  gsea_obj <- gsea_results[[clust]]
  if (is.null(gsea_obj) || nrow(gsea_obj) == 0) next
  
  p <- dotplot(gsea_obj, showCategory = 20, split = ".sign") +
    ggtitle(paste("Cluster", clust, "- Reactome GSEA")) +
    theme_minimal()
  if (Species == "Mouse"){
    out_file <- file.path(output_dir, paste0("MouseOnly_Cluster_", clust, "_Reactome_GSEA_dotplot.png"))
  }else{
    out_file <- file.path(output_dir, paste0("HumanOnly_Cluster_", clust, "_Reactome_GSEA_dotplot.png"))
  }
  png(out_file, width = 10, height = 8, units = "in", res = 300)
  print(p)
  dev.off()
}

### fibroblast/Epithelial bulk

### we've already removed epithelial cluster above so all cells here will be fibroblast in mouse
### epithelial in human

# Get raw counts
counts <- GetAssayData(VU40T.combined, slot = "counts")

# Metadata
meta <- VU40T.combined@meta.data

# Make sure sample_id and condition are present
stopifnot(all(c("sample", "condition") %in% colnames(meta)))

# Aggregate counts by sample
pseudobulk_counts <- sapply(unique(meta$sample), function(sid) {
  cell_ids <- rownames(meta)[meta$sample == sid]
  if (length(cell_ids) == 1) {
    counts[, cell_ids]
  } else {
    Matrix::rowSums(counts[, cell_ids])
  }
})

# Convert to matrix
pseudobulk_counts <- as.matrix(pseudobulk_counts)

# Sample-level metadata
sample_metadata <- meta %>%
  distinct(sample, condition) %>%
  filter(sample %in% colnames(pseudobulk_counts))

# Ensure sample order matches
sample_metadata <- sample_metadata[match(colnames(pseudobulk_counts), sample_metadata$sample), ]
# Replace dashes with dots (or underscores)
sample_metadata$condition <- make.names(sample_metadata$condition)


# Create DGE object
dge <- DGEList(counts = pseudobulk_counts, group = sample_metadata$condition)
dge <- calcNormFactors(dge)

# Design matrix
design <- model.matrix(~ 0 + sample_metadata$condition)
colnames(design) <- gsub("sample_metadata\\$condition", "", colnames(design))

# Estimate dispersion and fit
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

# Test LPS-P vs LPS-N
contrast <- makeContrasts(LPS.P - LPS.N, levels = design)
qlf <- glmQLFTest(fit, contrast = contrast)

# Top table
fibro_DE <- topTags(qlf, n = Inf)$table
if (Species == "Mouse"){
  write.csv(fibro_DE, file = file.path(git_dir, "Integrated/Mouse/MouseOnly_Fibroblast_bulk_DEs.csv"))
}else{
  write.csv(fibro_DE, file = file.path(git_dir, "Integrated/HumanOnly_Epithelial_bulk_DEs.csv"))
  
}

### now make volcano plot
if (Species == "Mouse"){
  output_dir <- file.path(git_dir, "Integrated/Mouse/Plots/Bulk_Analysis")
  out_file <- file.path(output_dir, "MouseOnly_Bulk_Fibroblast_volcano_FDR_0.1.png")
  celltype <- "Fibroblast" 
}else{
  output_dir <- file.path(git_dir, "Integrated/Plots/Bulk_Analysis")
  out_file <- file.path(output_dir, "HumanOnly_Bulk_Epithelial_volcano_FDR_0.1.png")
  celltype <- "Epithelial"
}
if (!(dir.exists(output_dir))){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}



# Create the pic
png(out_file, width = 8, height = 6, units = "in", res = 300)

select_labels <- head(rownames(fibro_DE[order(fibro_DE$FDR), ]), 10)
select_labels <- select_labels[!is.na(select_labels) & select_labels != ""]

print(EnhancedVolcano::EnhancedVolcano(
  fibro_DE,
  lab = rownames(fibro_DE),
  x = 'logFC',
  y = 'FDR',
  pCutoff = 0.1,
  FCcutoff = 1,
  title = paste0(celltype, ": LPS-P vs LPS-N"),
  pointSize = 2.0,
  labSize = 3,
  selectLab = select_labels,
  legendLabels = c('NS', 'Log2 FC', 'FDR', 'FDR & Log2 FC'),
  legendPosition = 'right',
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  max.overlaps = Inf	
))

dev.off()

### Bulk gsea

if (Species == "Mouse"){
  # Mouse Reactome gene sets from MSigDB
  msig_mouse <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "REACTOME")
  
  # Create list: pathway name → gene vector
  reactome_t2g <- msig_mouse[, c("gs_name", "gene_symbol")]   # <- data.frame with 2 columns
  
  rownames(fibro_DE) <- stringr::str_to_title(rownames(fibro_DE))
  
}


# Prepare ranked gene list
gene_list <- fibro_DE$logFC
names(gene_list) <- rownames(fibro_DE)
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!is.na(gene_list)]

# Run GSEA using gene symbols + Reactome gene sets
gsea <- tryCatch({
  GSEA(geneList = gene_list,
       TERM2GENE = reactome_t2g,
       pvalueCutoff = 0.1,
       verbose = FALSE)
}, error = function(e) NULL)


p <- dotplot(gsea, showCategory = 20, split = ".sign") +
  ggtitle(paste0(celltype, " LPS-P vs LPS-N - Reactome GSEA")) +
  theme_minimal()

out_file <- file.path(output_dir, paste0(Species,"Only", celltype, "_Reactome_GSEA_dotplot.png"))
png(out_file, width = 10, height = 12, units = "in", res = 300)
print(p)
dev.off()

### export epithelial data to h5ad building blocks for native pySCENIC
if (all(VU40T.combined$species == "Human")){
  library(Matrix)
  out <- file.path(proj_dir, "scRNAseqAnalysis/human_epi_mtx")
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  write.table(as.matrix(Matrix::t(GetAssayData(object = epi, slot = "counts"))), 
              file.path(out,'human_epi_counts.csv'), 
              sep = ',', row.names = T, col.names = T, quote = F)
}