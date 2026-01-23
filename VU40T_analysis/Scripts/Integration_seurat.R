Sys.setenv(BIOCFILECACHE_DIR = tempfile())

###libs
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)
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
input <- opt$input

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

human_cache <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")
mouse_cache <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", mirror = "useast")


### funcs


make_marker_dotplots <- function(
    seurat_obj,
    genelist,
    outPrefix,
    export = TRUE,
    plot_dir = ".",
    human_to_mouse_map = NULL,          # pass a cached mapping if you have it
    cluster_field = "seurat_clusters",  # which metadata field to group by
    resolution_tag = "0.6",             # used in filename
    max_width_in = 50                  # cap image width
) {
  # --- deps ---
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Package 'stringr' is required.")
  if (!requireNamespace("tools", quietly = TRUE)) stop("Package 'tools' is required.")
  
  # Offline ortholog mapping deps (replaces biomaRt)
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) stop("Package 'AnnotationDbi' is required for offline mapping.")
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("Package 'org.Hs.eg.db' is required for offline mapping.")
  if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) stop("Package 'org.Mm.eg.db' is required for offline mapping.")
  
  # --- output dir ---
  out_dir <- file.path(plot_dir, "MarkerDotplots")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # --- name handling (don’t mutate the input frame’s names) ---
  orig_names <- colnames(genelist)
  
  fix_second_token <- function(x) {
    parts <- stringr::str_split(x, "_")[[1]]
    if (length(parts) >= 2 && parts[2] != toupper(parts[2])) {
      parts[2] <- stringr::str_to_title(parts[2])
    }
    paste(parts, collapse = "_")
  }
  fixed_names   <- vapply(orig_names, fix_second_token, FUN.VALUE = character(1))
  display_names <- gsub("[._]", " ", fixed_names)
  safe_names    <- gsub("[^A-Za-z0-9]+", "_", display_names)
  
  # --- species flags (vector-safe) ---
  meta_species <- tryCatch(seurat_obj$species, error = function(e) NULL)
  is_mouse <- !is.null(meta_species) && all(meta_species == "Mouse", na.rm = TRUE)
  is_human <- !is.null(meta_species) && all(meta_species == "Human", na.rm = TRUE)
  
  mouse_features_upper <- is_mouse && all(rownames(seurat_obj) == toupper(rownames(seurat_obj)))
  
  # --- helper: offline human -> mouse mapping via HomoloGene IDs ---
  get_h2m_map_offline <- function(genes) {
    genes <- unique(genes[!is.na(genes) & nzchar(genes)])
    
    hs_db <- getNamespace("org.Hs.eg.db")$org.Hs.eg.db
    mm_db <- getNamespace("org.Mm.eg.db")$org.Mm.eg.db
    
    hs_cols <- AnnotationDbi::columns(hs_db)
    mm_cols <- AnnotationDbi::columns(mm_db)
    hs_kts  <- AnnotationDbi::keytypes(hs_db)
    mm_kts  <- AnnotationDbi::keytypes(mm_db)
    
    # Prefer HOMOLOGENE if available
    use_homologene <- ("HOMOLOGENE" %in% hs_cols) && ("HOMOLOGENE" %in% hs_kts) &&
      ("HOMOLOGENE" %in% mm_kts) && ("SYMBOL" %in% mm_cols)
    
    if (!use_homologene) {
      warning("Offline ortholog mapping via HOMOLOGENE not available in org.*.eg.db on this system. Proceeding without conversion.")
      return(NULL)
    }
    
    # Human SYMBOL -> HOMOLOGENE
    hs_map <- AnnotationDbi::select(
      hs_db,
      keys    = genes,
      keytype = "SYMBOL",
      columns = c("SYMBOL", "HOMOLOGENE")
    )
    hs_map <- hs_map[!is.na(hs_map$HOMOLOGENE) & nzchar(hs_map$HOMOLOGENE), , drop = FALSE]
    if (nrow(hs_map) == 0) return(NULL)
    
    # HOMOLOGENE -> Mouse SYMBOL
    mm_map <- AnnotationDbi::select(
      mm_db,
      keys    = unique(hs_map$HOMOLOGENE),
      keytype = "HOMOLOGENE",
      columns = c("HOMOLOGENE", "SYMBOL")
    )
    mm_map <- mm_map[!is.na(mm_map$SYMBOL) & nzchar(mm_map$SYMBOL), , drop = FALSE]
    if (nrow(mm_map) == 0) return(NULL)
    
    # Join human->homologene to mouse->symbol
    joined <- merge(
      hs_map[, c("SYMBOL", "HOMOLOGENE"), drop = FALSE],
      mm_map[, c("HOMOLOGENE", "SYMBOL"), drop = FALSE],
      by = "HOMOLOGENE",
      suffixes = c(".human", ".mouse")
    )
    
    # Build a named character vector: names = human SYMBOL, values = mouse SYMBOL
    # If multiple mouse symbols per human, keep the first (deterministic alphabetical)
    joined <- joined[order(joined$SYMBOL.human, joined$SYMBOL.mouse), , drop = FALSE]
    first_mouse <- tapply(joined$SYMBOL.mouse, joined$SYMBOL.human, function(x) x[1])
    
    out <- unname(first_mouse)
    names(out) <- names(first_mouse)
    out
  }
  
  # --- build mapping once if needed ---
  if (is_mouse && is.null(human_to_mouse_map)) {
    all_genes <- unique(unlist(genelist))
    if (length(all_genes) > 0) {
      human_to_mouse_map <- tryCatch(get_h2m_map_offline(all_genes), error = function(e) {
        warning(sprintf("Offline mapping failed: %s. Proceeding without conversion.", conditionMessage(e)))
        NULL
      })
    }
  }
  
  # If mouse object stores features in UPPERCASE, align mapping keys/values to uppercase
  if (is_mouse && !is.null(human_to_mouse_map) && mouse_features_upper) {
    human_to_mouse_map <- human_to_mouse_map[!is.na(names(human_to_mouse_map)) & nzchar(names(human_to_mouse_map))]
    names(human_to_mouse_map) <- toupper(names(human_to_mouse_map))
    human_to_mouse_map <- toupper(unname(human_to_mouse_map)) |> setNames(toupper(names(human_to_mouse_map)))
  }
  
  # --- plotting loop ---
  for (i in seq_along(genelist)) {
    genes_i <- unique(genelist[[i]])
    genes_i <- genes_i[!is.na(genes_i) & nzchar(genes_i)]
    if (length(genes_i) == 0) next
    
    # If mouse object and we have a mapping, convert human→mouse
    if (is_mouse && !is.null(human_to_mouse_map)) {
      if (mouse_features_upper) genes_i <- toupper(genes_i)
      genes_i <- unname(human_to_mouse_map[genes_i])
      genes_i <- unique(genes_i[!is.na(genes_i) & nzchar(genes_i)])
    }
    
    # Keep only genes present in the dataset
    genes_i <- intersect(genes_i, rownames(seurat_obj))
    if (length(genes_i) == 0) next
    
    # Build plot
    p <- Seurat::DotPlot(seurat_obj, features = genes_i, group.by = cluster_field) +
      Seurat::RotatedAxis() +
      ggplot2::ggtitle(display_names[i]) +
      ggplot2::scale_color_gradient(low = "lightgrey", high = "red") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9)
      ) +
      ggplot2::labs(x = "Genes", y = "Clusters")
    
    # Title-case gene labels for mouse display (cosmetic)
    if (is_mouse) {
      p <- p + ggplot2::scale_x_discrete(labels = function(x) stringr::str_to_title(x))
    }
    
    # Export or print
    if (export) {
      n_genes <- length(genes_i)
      plot_width <- if (n_genes <= 6) 4 else min(max(n_genes / 3, 6), max_width_in)
      
      # filename pieces
      species_tag <- if (is_human) "HumanONLY" else "Mouse"
      fname <- sprintf(
        "%s_Dotplot_clusterResolution%s_%s_%s_VU40T_combined.png",
        species_tag, resolution_tag, outPrefix, safe_names[i]
      )
      
      png(file = file.path(out_dir, fname),
          width = plot_width, height = 4, units = "in", res = 300)
      print(p)
      dev.off()
    } else {
      print(p)
    }
  }
}


if (is.null(input)){
  preIntegrationSeuratList <- readRDS(file.path(cache_dir, paste0("VU40T_singlets_only_sep_samples_",species,".RDS")))
} else {
  preIntegrationSeuratList <- readRDS(input)
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

# preIntegrationSeuratList <- readRDS("annotated_VU40T_sepsamples_singlets_only_mouse.rds")
# 
# preIntegrationSeuratList <- readRDS("R/VU40T_singlets_only_sep_samples.RDS")


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


### pre-integration QC

# Compute per-sample pseudobulk vector (gene means across cells)
pb_mat <- sapply(names(preIntegrationSeuratList), function(nm) {
  obj <- preIntegrationSeuratList[[nm]]
  
  DefaultAssay(obj) <- "RNA"
  
  x <- GetAssayData(obj, slot = "data")
  
  # gene-wise mean across cells (works for sparse matrices too)
  if (inherits(x, "dgCMatrix")) {
    Matrix::rowMeans(x)
  } else {
    rowMeans(as.matrix(x))
  }
})

# Ensure matrix shape: genes x samples
pb_mat <- as.matrix(pb_mat)

#pearson correlation
cor_mat <- cor(pb_mat, method = "pearson", use = "pairwise.complete.obs")

# Color mapping
col_fun <- circlize::colorRamp2(
  c(0.9,0.95, 1),
  c("#2166AC", "#F7F7F7", "#B2182B")
)

# ComplexHeatmap plot
ht <- Heatmap(
  cor_mat,
  name = "r",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  rect_gp = grid::gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(at = c(0.9, 0.95, 1)#,
                              # title_gp  = grid::gpar(fontsize = 8),
                              # labels_gp = grid::gpar(fontsize = 7)
                              ),
  # cell_fun = function(j, i, x, y, w, h, fill) {
  #   grid::grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = grid::gpar(fontsize = 10))
  # },
  column_title = "Intersample Whole Transcriptome Correlation",
  row_title = NULL
)

png(filename = file.path(plots_dir, 
                         paste0("IntersampleWTranscriptomeCorHeatMap_", species, "_only.png")),
    height = 4 ,width = 6, units = "in", res = 300)
draw(ht)
dev.off()

## diagnostic pca

merged <- merge(
  x = preIntegrationSeuratList[[1]],
  y = preIntegrationSeuratList[-1],
  add.cell.ids = names(preIntegrationSeuratList),
  project = "merged_preintegration"
)

merged$condition <- ifelse(grepl("LPS-N", merged$sample), "LPS-N", "LPS-P")
merged$passage <- gsub("_.*", "", merged$sample)


DefaultAssay(merged) <- "RNA"
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged, features = VariableFeatures(merged))
merged <- RunPCA(merged, features = VariableFeatures(merged), npcs = 50)


# change 'lps_status' to your metadata column name if needed
png(
  filename = file.path(plots_dir, "diagnostic_PCA_by_LPS_status.png"),
  width = 6,
  height = 5,
  units = "in",
  res = 300
)
DimPlot(merged, reduction = "pca", group.by = "condition", pt.size = 0.2) +
  ggtitle("Diagnostic PCA: colored by LPS treatment status")
dev.off()
png(
  filename = file.path(plots_dir, "diagnostic_PCA_by_sample.png"),
  width = 6,
  height = 5,
  units = "in",
  res = 300
)

DimPlot(merged, reduction = "pca", group.by = "sample", pt.size = 0.2) +
  ggtitle("Diagnostic PCA: colored by LPS treatment status")

dev.off()

rm(merged)
qc_df <- bind_rows(lapply(names(preIntegrationSeuratList), function(s) {
  
  obj <- preIntegrationSeuratList[[s]]
  
  data.frame(
    sample = s,
    nFeature_RNA = obj$nFeature_RNA,
    nCount_RNA   = obj$nCount_RNA,
    subsets_Mt_percent = obj$subsets_Mt_percent
  )
}))

qc_long <- qc_df %>%
  tidyr::pivot_longer(
    cols = c(nFeature_RNA, nCount_RNA, subsets_Mt_percent),
    names_to = "metric",
    values_to = "value"
  )

## pre integration QC plots
png(file = file.path(plots_dir, 
                     paste0(species,"Only_preIntegration_QC_violins.png")),
    width = 8, height = 4, units = "in", res = 300)

p <- ggplot(qc_long, aes(x = sample, y = value, fill = sample)) +
    geom_violin(trim = TRUE) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    theme_bw() +
    theme(legend.position = "none") +
    ylab("Value") +
    xlab("Sample") +
    ggtitle("QC metrics across samples (pre-integration)")
print(p)
dev.off()

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
if (species == "mouse"){
  optimal_pc <- 20  ## 15 human, 16 mouse - as per jackstraw plot when using all data together. 20 for mouse seperately, 10 for human seperately    
} else {
  
  optimal_pc <- 10  ## 15 human, 16 mouse - as per jackstraw plot when using all data together. 20 for mouse seperately, 10 for human seperately
}

    

# optimal_pc <- 15  ## 15 human, 16 mouse - as per jackstraw plot when using all data together. 20 for mouse seperately, 10 for human seperately    


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

if (species == "mouse"){
  png(file = file.path(plots, "/Integrated/Mouse/Plots/MouseOnly_20_PCs_ClustTree_VU40T_stabilityscores_integrated.png"), width = 16, height = 12, units = "in", res = 300)
  p1 <- clustree(VU40T.combined, prefix = "RNA_snn_res.", node_colour = "sc3_stability") ## best resolution lies between 0.6-0.8
  print(p1)
  dev.off()
  
  png(file = file.path(git_dir, "/Integrated/Mouse/Plots/MouseOnly_20_PCs_ClustTree_VU40T_integrated.png"), width = 16, height = 12, units = "in", res = 300)
  p1 <- clustree(VU40T.combined, prefix = "RNA_snn_res.") ## best resolution lies between 0.6-0.8
  print(p1)
  dev.off()
}



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
# VU40T.combined <- RunTSNE(VU40T.combined, reduction = "pca", dims = 1:optimal_pc, verbose = F, seed.use = 666)
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

# saveRDS(VU40T.combined, file = file.path(cache_dir, "VU40T_combined_joined_1_res_humanOnly.rds"))
if (species == "human"){
  optimum_res <- 0.3
  saveRDS(VU40T.combined, file = file.path(cache_dir, "VU40T_combined_joined_0.3_res_humanOnly.rds"))
} else {
  optimum_res <-0.6
  saveRDS(VU40T.combined, file = file.path(cache_dir, "VU40T_combined_joined_0.6_res_mouseOnly.rds"))
}


##mouse
# VU40T.combined<- readRDS("VU40T_combined_joined_1.0_res_mouse.rds")
##human
# VU40T.combined<- readRDS("~/R/VU40T_joined_res0.5.rds")
##VU40T.combined<- readRDS(file.path(cache_dir, "VU40T_combined_joined_1_res_humanOnly.rds")) ## overclustered
## from crash
# if (species == "human"){
  # VU40T.combined <- readRDS(file.path(cache_dir, "VU40T_combined_joined_0.3_res_humanOnly.rds"))
# } else {
  # VU40T.combined <- readRDS(file.path(cache_dir, "VU40T_combined_joined_0.6_res_mouseOnly.rds"))
# }
##find markers
##mouse only

### pre and post Integration QC

merged_pre <- merge(
  preIntegrationSeuratList[[1]],
  y = preIntegrationSeuratList[-1],
  add.cell.ids = names(preIntegrationSeuratList),
  project = "pre_integration"
)

DefaultAssay(merged_pre) <- "RNA"
merged_pre <- NormalizeData(merged_pre)
merged_pre <- FindVariableFeatures(merged_pre)
merged_pre <- ScaleData(merged_pre)
merged_pre <- RunPCA(merged_pre)
merged_pre <- RunUMAP(merged_pre, dims = 1:optimal_pc)

p1 <- DimPlot(merged_pre, reduction = "umap", group.by = "sample") +
  ggtitle("Before integration (RNA)")

DefaultAssay(VU40T.combined) <- "integrated"
p2 <- DimPlot(VU40T.combined, reduction = "umap", group.by = "sample") +
  ggtitle("After integration (CCA)")

png(file = file.path(plots_dir, 
                     paste0(species, 
                            "only_",optimal_pc,"PCs_UMAPs_pre_and_post_integration.png")),
    width = 10, height = 5, units = "in", res = 300)

print(p1+p2)
dev.off()

## implementing lisi as a metric here for batch mixing

emb_int <- Embeddings(VU40T.combined, reduction = "pca")[, 1:10, drop = FALSE]
meta_int <- VU40T.combined@meta.data

lisi_int <- lisi::compute_lisi(
  X = emb_int,
  meta_data = meta_int,
  label_colnames = c("sample")   # or "orig.ident" if that's your column
)



lisi_scores <- lisi::compute_lisi(emb_int, 
                                  meta_int, 
                                  label_colnames = c("sample", "seurat_clusters"))
summary(lisi_scores)

DefaultAssay(merged_pre) <- "RNA"

emb_pre <- Embeddings(merged_pre, reduction = "pca")[, 1:10, drop = FALSE]
meta_pre <- merged_pre@meta.data

lisi_pre <- lisi::compute_lisi(
  X = emb_pre,
  meta_data = meta_pre,
  label_colnames = c("sample")
)


df_lisi <- rbind(
  data.frame(stage = "Before integration (RNA)", iLISI = lisi_pre$sample),
  data.frame(stage = "After integration (CCA)",  iLISI = lisi_int$sample)
)

df_lisi$stage <- factor(df_lisi$stage,  
                           levels = c("Before integration (RNA)", "After integration (CCA)"))

wil <- wilcox.test(iLISI ~ stage, data = df_lisi, exact = FALSE)


p <- ggplot(df_lisi, aes(x = stage, y = iLISI)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    ylab("iLISI (sample mixing)") +
    xlab(NULL)

p <- p + ggpubr::stat_compare_means(
  method = "wilcox.test",
  comparisons = list(c("Before integration (RNA)", "After integration (CCA)")),
  label = "p.format"
)

print(p)
ggsave(filename = file.path(plots_dir, 
                            paste0(species,"only_lisi_BeforeandAfter_wilcoxon.png")),
        width = 4, height = 5, units = "in")
# compute medians
meds <- tapply(df_lisi$iLISI, df_lisi$stage, median)

effect_table <- data.frame(
  stage = names(meds),
  median_iLISI = as.numeric(meds)
)

# add delta (after - before)
effect_table$delta_from_before <- effect_table$median_iLISI -
  effect_table$median_iLISI[effect_table$stage == "Before integration (RNA)"]

write.csv(effect_table, file = file.path(res_dir, "iSLI_effectsize.csv"), row.names = F)


rm(preIntegrationSeuratList)

rm(merged_pre)

DefaultAssay(VU40T.combined) <- "RNA"

if (all(VU40T.combined$species == "Mouse")){
  Species <- "Mouse"
} else{
  Species <- "Human"
}

markers <- FindAllMarkers(VU40T.combined,
                          only.pos = T,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)




write.csv(markers, file = file.path(
  res_dir, 
  paste0(species,"Only",optimal_pc,"_PCs_ClusterMarkers_res_",optimum_res,"_VU40T_combined.csv")), 
  row.names = FALSE)

## convert back to camelcase if mouse.
if (all(VU40T.combined$species == "Mouse")){
  markers$gene <- sapply(markers$gene, stringr::str_to_title)
}




p1 <- DimPlot(VU40T.combined, reduction = "umap", label = T) 

if (all(VU40T.combined$species == "Mouse")){
  png(file = file.path(plots_dir, "MouseOnly_20_PCs_VU40T_UMAP_Resolution_0.6.png"), width = 8, height = 5, units = "in", res = 300)
}else{
  png(file = file.path(plots_dir, "HumanOnly_10_PCs_VU40T_UMAP_Resolution_0.3.png"), width = 8, height = 5, units = "in", res = 300)
}

print(p1)
dev.off()

# annotation validation - not run
# library(HCATonsilData)
# 
# 
# library("BiocFileCache")
# 
# Sys.setenv(BFC_CACHE=file.path(proj_dir, ".cache"))
# # Example: Set the cache to a specific directory
# bfc = BiocFileCache(cache = file.path(proj_dir, ".cache"), ask = FALSE)
# bfccache(bfc)
# 
# bfccache(bfc)
# 
# tools::R_user_dir("BiocFileCache", which="cache")
# 
# if (all(VU40T.combined$species == "Human")){
#   sce_HCATonsil_epithelial <- HCATonsilData(assayType = "RNA", cellType = "epithelial")
#   # Get counts and metadata manually
#   counts_epi <- counts(sce_HCATonsil_epithelial)
#   counts_epi <- as(counts_epi, "dgCMatrix")
#   logcounts_epi <- logcounts(sce_HCATonsil_epithelial)
#   logcounts_epi <- as(logcounts_epi, "dgCMatrix")
#   metadata_epi <- as.data.frame(colData(sce_HCATonsil_epithelial))
#   seu_epi <- CreateSeuratObject(counts = counts_epi, meta.data = metadata_epi)
#   seu_epi <- SetAssayData(seu_epi, slot = "data", new.data = logcounts_epi)
#   seu_ref <- seu_epi
#   rm(list = c(seu_epi, counts_epi, logcounts_epi, metadata_epi))
# }
# 
# if (all(VU40T.combined$species == "Mouse")){
#   sce_HCATonsil_mesenchymal <- HCATonsilData(assayType = "RNA", cellType = "FDC")
#   counts_fdc <- counts(sce_HCATonsil_mesenchymal)
#   logcounts_fdc <- logcounts(sce_HCATonsil_mesenchymal)
#   counts_fdc <- as(counts_fdc, "dgCMatrix")
#   logcounts_fdc <- as(logcounts_fdc, "dgCMatrix")
#   metadata_fdc <- as.data.frame(colData(sce_HCATonsil_mesenchymal))
#   seu_fdc <- CreateSeuratObject(counts = counts_fdc, meta.data = metadata_fdc)
#   seu_fdc <- SetAssayData(seu_fdc, slot = "data", new.data = logcounts_fdc)
#   seu_ref <- seu_fdc
#   rm(list = c(seu_fdc, counts_fdc, logcounts_fdc, metadata_fdc))
# }
# 
# seu_ref <- NormalizeData(seu_ref)
# seu_ref <- FindVariableFeatures(seu_ref)
# seu_ref <- ScaleData(seu_ref)
# seu_ref <- RunPCA(seu_ref)
# 
# anchors <- FindTransferAnchors(
#   reference = seu_ref,
#   query = VU40T.combined,
#   dims = 1:30,
#   reference.assay = "RNA",
#   query.assay = "RNA"
# )
# 
# predictions <- TransferData(
#   anchorset = anchors,
#   refdata = seu_fdc$annotation_20230508,  # or ref_seurat$celltype if that's what it's called
#   dims = 1:30
# )
# VU40T.combined <- AddMetaData(object = VU40T.combined, metadata = predictions)
# 
# p1 <- DimPlot(VU40T.combined, group.by = "predicted.id", split.by = "condition" , label = TRUE, label.size = 3)
# p2 <- DimPlot(VU40T.combined, group.by = "seurat_clusters", split.by = "condition")
# if(all(VU40T.combined$species == "Mouse")){
#   png(file = paste0(git_dir,"/Integrated/Mouse/Plots/Azimuth_annotated_UMAP_VU40T_combined.png"), width = 12, height = 10, units = "in", res = 300)
# }else{
#   png(file = paste0(git_dir,"/Integrated/Plots/Azimuth_annotated_UMAP_VU40T_combined.png"), width = 12, height = 10, units = "in", res = 300)
# }
# print(p1 + p2)
# dev.off()
# 
# ## for mouse
# if (all(VU40T.combined$species == "Mouse")) {
#   colnames(VU40T.combined) = sapply(colnames(VU40T.combined), stringr::str_to_title)
#   ref <- MouseRNAseqData(cell.ont = "nonna")
# }else{
#   ref <- readRDS("~/R/HPCA_reference.rds") 
# }
# 
# 
# cluster_labels <- VU40T.combined$seurat_clusters
# bp <- BiocParallel::MulticoreParam(workers = 8) 
# # Run SingleR annotation
# sce <- Seurat::as.SingleCellExperiment(VU40T.combined)
# 
# ### doesn't work with human-aligned integrated. kinda works with mouse
# ## doesn't work with mouse only
# pred <- SingleR::SingleR(
#   test = sce,
#   ref = ref,
#   labels = ref$label.main,
#   clusters =  cluster_labels,
#   assay.type.test = 1,
#   BPPARAM=bp
# )
# 
# 
# 
# # Apply cluster annotations to each cell
# VU40T.combined[["SingleR_HTA_cluster_label"]] <- pred$labels[
#   match(VU40T.combined$seurat_clusters, rownames(pred))
# ]
# 
# png(file = file.path(git_dir, "Integrated/Plots/HumanOnly_10_PCs_SingleR_annotated_VU40T_combined_res_0.3_labelMain.png"), width = 10, height = 8, units = "in", res = 300)
# 
# # Visualize cell type annotation on UMAP
# p <- Seurat::DimPlot(VU40T.combined, group.by='SingleR_HTA_cluster_label', split.by = "condition")
# p2 <- DimPlot(VU40T.combined, reduction = "umap", split.by = "condition", group.by = "seurat_clusters", label = T) 
# print(p + p2)
# #print (p2)
# dev.off()


### marker dotplots
genelist <- read.csv("~/Markers_for_dotplots_2_set1.csv", header = T) ## both human and mouse markers

if(species == "mouse"){
  VU40T.combined <- subset(VU40T.combined, idents = 13, invert = TRUE)
}



make_marker_dotplots(VU40T.combined, 
                     genelist, 
                     export = T, 
                     outPrefix = paste0(species, "_onlySET1"), 
                     plot_dir = plots_dir,
                     resolution_tag = optimum_res)

## genelist 2 -> human only
if (all(VU40T.combined$species == "Human")){
  genelist <- read.csv("~/markers_secondset_human_ep.csv", header = T)
  make_marker_dotplots(VU40T.combined, 
                       genelist, 
                       export = T, 
                       outPrefix = paste0(species, "_onlySET2_Epithelial_res",optimum_res), 
                       plot_dir = plots_dir)
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
    plots_dir,
    paste0("HeatMap_epiMarkers_clusterResolution",optimum_res,"_VU40T_combined.png")),
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
    plots_dir,
    paste0("Marker_Overlay_UMAPs_EpiMarkers_Viridis_clusterResolution",optimum_res,"_VU40T_combined.png")),
    width = 14, height = 12, units = "in", res = 300)
  p <- FeaturePlot(VU40T.combined, 
                   features = umap_markers, 
                   label = F,# label.size = 3, repel = T,
  ) & 
    scale_color_viridis_c()
  print(p)
  dev.off()
  ### new human dotplots
  genelist <- readxl::read_excel(file.path(proj_dir,"scRNAseqAnalysis/Markers_for_dotplots_2.xlsx"), sheet = 8)
  make_marker_dotplots(VU40T.combined, 
                       genelist, 
                       export = T, 
                       outPrefix = paste0(species, "_onlySET3_Epithelial",optimum_res), 
                       plot_dir = plots_dir)
  
  ### and violin plots
  genelist <- readxl::read_excel(file.path(proj_dir,"scRNAseqAnalysis/Markers_for_dotplots_2.xlsx"), sheet = 9,col_names = F)
  
  
  for (col in seq_len(ncol(genelist))) {
    png(file = file.path(
      plots_dir,
      paste0("ViolinPlot_EpiMarkers_set3_clusterResolution",optimum_res,"_VU40T_combined",col,".png")),
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
    dev.off()
  }
  ### new human dotplots -> HNSCC and EMT
  genelist <- readxl::read_excel(file.path(proj_dir,"scRNAseqAnalysis/Markers_for_dotplots_3.xlsx"), sheet = 11, skip = 1)
  make_marker_dotplots(VU40T.combined, genelist, export = T, outPrefix = paste0(species, "_Fig4_Epithelial",optimum_res), plot_dir = plots_dir)
  
}
## genelist 3 -> mouse and human epithelial, fibroblast and EMT markers
genelist <- as.data.frame(read.csv("~/Markers_for_dotplots_2_ep_2nd_set.csv", header = T, skip = 1))
genelist <- as.data.frame(lapply(genelist, toupper))
genelist[] <- lapply(genelist, function(x) gsub("CDH2", "", x))

if (all(VU40T.combined$species == "Mouse")){
  genelist <- genelist[, c(4,5)] ## both human and mouse markers w/out emt
  colnames(genelist) <- gsub("_[HM]", "", colnames(genelist))
  genelist <- apply(genelist, stringr::str_to_title)
  
}else{
  genelist <- genelist[, c(1:2)] ## both human and mouse markers, no need for EMT here
  colnames(genelist) <- gsub("_[HM]", "", colnames(genelist))
}
colnames(genelist) <- gsub("_", " ", colnames(genelist))
colnames(genelist) <- gsub("markers", "Markers", colnames(genelist))
## convert to named vectors
genevectors <- lapply(genelist, function(x) unique(x[x != "" & !is.na(x)]))

if(species == "mouse"){
  genevectors <- lapply(genevectors, stringr::str_to_title)
}


png(file = file.path(
  plots_dir,
  paste0("CombinedDotplot_epi_fibro_Markers_clusterResolution",optimum_res,"_VU40T_combined.png")),
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

if (all(VU40T.combined$species == "Mouse")){

  
  ## make stacked barplot of num fibroblasts per clust
  # clust_cell_df <- as.data.frame(table(VU40T.combined$seurat_clusters))
  # colnames(clust_cell_df) <- c("Cluster", "No. Cells")
  # write.csv(clust_cell_df, file = file.path(git_dir, "Integrated/Mouse/MouseOnly_fibroblast_cells_per_clust.csv"))
  # clust_cell_df$`Cell Type` <- "Fibroblast"
  # ## plot barplot
  # 
  # png(file = file.path(
  #   plots_dir,
  #   "Integrated/Mouse/Plots/MouseONLY_Stacked_Bar_numCells_per_cluster_Resolution",optimum_res,"_VU40T_combined.png"),
  #   width = 4, height = 6, units = "in", res = 300)
  # 
  # p <- clust_cell_df %>% 
  #   ggplot(aes(x = `Cell Type`, y = `No. Cells`, fill = Cluster)) + 
  #   geom_bar(position="stack", stat="identity") + 
  #   theme_minimal() + xlab("")
  # print(p)
  # dev.off()
  
  # png(file = file.path(
  #   git_dir,
  #   "Integrated/Mouse/Plots/MouseONLY_Stacked_Bar_percent_Cells_per_cluster_Resolution0.6_VU40T_combined.png"),
  #   width = 4, height = 6, units = "in", res = 300)
  # p <- clust_cell_df %>% 
  #   ggplot(aes(x = `Cell Type`, y = `No. Cells`, fill = Cluster)) + 
  #   geom_bar(position="fill", stat="identity") + 
  #   theme_minimal() +
  #   ylab("Proportion of Cells") + xlab("") + 
  #   scale_y_continuous(labels = scales::percent)
  # print(p)
  # dev.off()
  # genelist <- read.csv("~/Fibroblasts_dotplots_markers.csv", header = F)
  # names(genelist) <- "Fibroblast Markers"
  # 
  # png(file = file.path(
  #   git_dir,
  #   "Integrated/Mouse/Plots/MouseONLY_MarkerDotplots/Dotplot_Fibro_Markers_clusterResolution0.6_VU40T_combined.png"),
  #   width = 12, height = 4, units = "in", res = 300)
  # 
  # p <- DotPlot(
  #   VU40T.combined,
  #   features = genelist,
  #   group.by = "seurat_clusters",
  # ) +
  #   scale_x_discrete(labels = function(x) stringr::str_to_title(x)) +
  #   RotatedAxis() +
  #   scale_color_gradient(low = "lightgrey", high = "red") +
  #   theme(
  #     axis.text.x = element_text(angle = 45, hjust = 1),
  #     axis.text.y = element_text(face = "italic")  # Optional: italics for gene names
  #   ) + labs(y = "Clusters")  # 👈 Relabel y-axis
  # print(p)
  # dev.off()
  
  ## fibroblast heatmap
  # heatmap_markers <- read.csv("~/fibroblast_heatmap.csv", header = T)
  fibro2_marker_ls <- list()
  for (i in seq(1,6)){
    fibro2_marker_ls[[i]] <- readxl::read_xlsx(file.path(proj_dir ,"/scRNAseqAnalysis/Fibroblasts_dotplots.xlsx"), sheet = i)
  }
  heatmap_markers <- fibro2_marker_ls[[2]]
  ### for heatmap
  # heatmap_markers <- genelist[, c(1:5)]
  colnames(heatmap_markers) <- stringr::str_replace(colnames(heatmap_markers), "_", " ")
  ## for fibro list 2
  colnames(heatmap_markers) <- stringr::str_replace(colnames(heatmap_markers), "Oxid phosph", "Oxidative phosphorylation")
  colnames(heatmap_markers) <- stringr::str_to_title(colnames(heatmap_markers))
  colnames(heatmap_markers) <- stringr::str_replace(colnames(heatmap_markers), "Ecm", "ECM")
  colnames(heatmap_markers) <- stringr::str_replace(colnames(heatmap_markers), "To", "to")
  heatmap_markers <- as.data.frame(apply(heatmap_markers, 2 , stringr::str_to_title))
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
    circlize::rand_color(length(lev)),
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
    plots_dir,
    paste0("Heatmap_FibroMarkers_set2_clusterResolution",optimum_res,"_VU40T_combined.png")),
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
  ViolinGenes <- stringr::str_to_title(ViolinGenes)
  png(file = file.path(
    plots_dir,
    paste0("ViolinPlot_FibroMarkers_set2_clusterResolution",optimum_res,"_VU40T_combined.png")),
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
    p[[i]][["labels"]][["x"]] <- ""
    p[[i]][["labels"]][["y"]] <- ""
  } 
  
  p_stack <- wrap_plots(p, ncol = 2, nrow = 5,  guides = "collect")
  
  print(p_stack)
  dev.off()
  
  ### v2 violinplots
  
  
  ViolinGenes <- as.data.frame(fibro2_marker_ls[[5]])
  ViolinGenes <- as.data.frame(apply(ViolinGenes, 2, stringr::str_to_title))
  
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
        plots_dir,
        paste0(
          "V2_ViolinPlot_FibroMarkers_",
          gsub("[^A-Za-z0-9._-]+", "_", panel_chr),
          "_clusterResolution", optimum_res, "_VU40T_combined.png"
        )
      ),
      width = 5, height = h_in, units = "in", res = 300
    )
    print(p_stack)
    dev.off()
  }
  
  
  
  ### overlay umaps
  
  umap_markers <- fibro2_marker_ls[[6]]
  colname <- colnames(umap_markers)[1]
  umap_markers <- c(colname, as.character(umap_markers[[1]]))
  umap_markers<- append(umap_markers, "ZBTB7B")
  umap_markers <- stringr::str_to_title(umap_markers)
  
  png(file = file.path(
    plots_dir,
    paste0("Marker_Overlay_UMAPs_FibroMarkers_set2_clusterResolution",optimum_res,"_VU40T_combined.png")),
    width = 15, height = 25, units = "in", res = 300)
  p <- FeaturePlot(VU40T.combined, 
                   features = umap_markers, 
                   label = F,# label.size = 3, repel = T,
  ) & 
    scale_color_viridis_c()
  
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
    plots_dir,
    paste0("CombinedDotplot_CSC_EMT_Markers_clusterResolution",optimum_res,"_VU40T_combined.png")),
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
    plots_dir,
    paste0("ViolinPlot_CSC_EPI_Markers_clusterResolution",optimum_res,"_VU40T_combined.png")),
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
    plots_dir,
    paste0("ViolinPlot_Basal_Markers_clusterResolution",optimum_res,"_VU40T_combined.png")),
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









# #### epithelial cluster 6 inspection
# 
# DimPlot(VU40T.combined, label = T)
# 
# 
# ### look at QC fail
# 
# if (species == "human"){
#   mtx_dir <- file.path(proj_dir, "BaseSpace/LPS_VU40T_QC_and_counts/cellranger/mtx_conversions")
# } else {
#   mtx_dir <- file.path(proj_dir, "BaseSpace/LPS_VU40T_QC_and_counts_Mouse/cellranger/mtx_conversions")
# }
# 
# Samples <- list.files(path = mtx_dir, 
#                       pattern = "^P",
#                       recursive = F
# ) 
# 
# 
# 
# ##check we only have 4 samples
# ## get only cellbender seurat_obj paths
# seurat_objs <- list.files(path = mtx_dir, 
#                           pattern = "^P.*cellbender.*.seurat.rds$",
#                           full.names = T,
#                           recursive = T
# )
# 
# 
# 
# SeuratList <- list()
# for(i in seurat_objs){
#   SeuratList[[i]] <- LoadSeuratRds(i)
#   Idents(SeuratList[[i]]) <- SeuratList[[i]]$sample
# }
# 
# 
# # Replace SeuratList names with simplified format like P18_LPS-N, P22_LPS-P, etc.
# names(SeuratList) <- gsub("^.*/(P[0-9]+)_LPS-([NP]).*$", "\\1_LPS-\\2", names(SeuratList))
# 
# ## Mito genes were isolated from hg38 GTF (chrM)
# 
# has_mito_prefix <- any(sapply(SeuratList, function(obj) {
#   any(grepl("^MT-", rownames(obj[["RNA"]]))) |
#     any(grepl("^mt-", rownames(obj[["RNA"]])))
# }))
# 
# if (!has_mito_prefix) {
#   mito_path <- "/rds/projects/g/gendood-3dmucosa/BaseSpace/LPS_VU40T_QC_and_counts/cellranger/mkref/cellranger_reference/genes/MTgenes.txt"
#   
#   if (file.exists(mito_path)) {
#     mito_genes <- read.delim(mito_path, sep = "\t", header = FALSE)
#     mito_genes <- c(t(mito_genes))
#   } else {
#     stop("❌ Mitochondrial gene list not found at: ", mito_path)
#   }
# }
# 
# 
# ## thresholds worked out from doing histogram of log-transformed count data
# thresholds_df <- data.frame(
#   sample = c("P18_LPS-N", "P18_LPS-P", "P22_LPS-N", "P22_LPS-P"),
#   cutoff = c(300, 350, 300, 300)
# )
# 
# seurat_filtered_list <- list() 
# all_qc_metadata <- list()
# for (sample in seq_along(SeuratList)) {
#   seurat_obj <- SeuratList[[sample]]
#   sce <- as.SingleCellExperiment(seurat_obj,assay = "RNA")
#   sample_id <- unique(seurat_obj$sample)
#   
#   # Detect mitochondrial genes
#   if (! has_mito_prefix) {
#     mito_flag <- rownames(sce) %in% mito_genes
#   } else if (any(grepl("^MT-", rownames(seurat_obj[["RNA"]]))) &
#              !any(grepl("^mt-", rownames(seurat_obj[["RNA"]])))
#   ) {
#     mito_genes_detected <- rownames(seurat_obj[["RNA"]])[grepl("^MT[-]?", rownames(seurat_obj[["RNA"]]))]
#     mito_flag <- rownames(sce) %in% mito_genes_detected
#   } else {
#     mito_genes_detected <- rownames(seurat_obj[["RNA"]])[grepl("^mt[-]?", rownames(seurat_obj[["RNA"]]))]
#     mito_flag <- rownames(sce) %in% mito_genes_detected
#   }
#   
#   # Calculate QC metrics
#   qc_metrics <- scuttle::perCellQCMetrics(sce, subsets = list(Mt = mito_flag))
#   
#   # Apply QC thresholds
#   qc_metrics$low_lib <- isOutlier(qc_metrics$sum, log = TRUE, type = "lower", nmads = 5)
#   qc_metrics$low_feats <- isOutlier(qc_metrics$detected, log = TRUE, type = "lower", nmads = 5) |
#     qc_metrics$sum < thresholds_df$cutoff[thresholds_df$sample == sample_id]
#   qc_metrics$high_mito <- qc_metrics$subsets_Mt_percent > 20
#   qc_metrics$qc_pass <- !(qc_metrics$low_feats | qc_metrics$high_mito)
#   
#   # Add metadata
#   qc_metrics <- as.data.frame(qc_metrics)
#   qc_metrics <- qc_metrics[colnames(seurat_obj), , drop = FALSE]
#   seurat_obj <- AddMetaData(seurat_obj, metadata = qc_metrics)
#   seurat_obj$qc_status <- ifelse(seurat_obj$qc_pass, "Pass", "Fail")
#   
#   # Subset Seurat object
#   meta_pass <- seurat_obj@meta.data %>% filter(qc_pass)
#   ### let's get qc fail this time
#   seurat_filtered <- subset(seurat_obj, cells = setdiff(colnames(seurat_obj), rownames(meta_pass)))
#   seurat_filtered_list[[sample]] <- seurat_filtered
#   
#   # Store for combined plots
#   qc_metrics$sample_id <- sample_id
#   qc_metrics$nFeature_RNA <- seurat_obj$nFeature_RNA[rownames(qc_metrics)]
#   qc_metrics$percent_mt <- qc_metrics$subsets_Mt_percent
#   qc_metrics$qc_status <- seurat_obj$qc_status[rownames(qc_metrics)]
#   all_qc_metadata[[sample]] <- qc_metrics
#   
#   # Print summary
#   message("\nQC summary for ", sample_id, "\n")
#   print(qc_metrics %>% summarise(
#     total_cells = n(),
#     kept = sum(qc_pass),
#     percent_kept = mean(qc_pass) * 100,
#     low_feats = sum(low_feats),
#     high_mito = sum(high_mito),
#     low_lib = sum(low_lib)
#   ))
# }
# 
# ## reassign sample name to slices of seurat obj list:
# names(seurat_filtered_list) <- names(SeuratList)
# 
# 
# # Combine all metadata
# qc_df <- do.call(rbind, all_qc_metadata)
# 
# # Plot combined histogram: nFeature_RNA
# p1 <- ggplot(qc_df, aes(x = detected, fill = qc_pass)) +
#   geom_histogram(alpha = 0.5, bins = 100, position = "identity") +
#   scale_x_log10() +
#   theme_minimal() +
#   facet_wrap(~sample_id, scales = "free_y") +
#   labs(title = "nFeature_RNA: Pass vs Fail by Sample",
#        x = "log10(Number of detected genes)", y = "Cell count")
# print(p1)
# 
# ## now let's look at these qc fail cells
# qc_fail_cell_counts_df <- data.frame(
#   sample = names(SeuratList),
#   n_cells = sapply(SeuratList, ncol)
# )
# 
# for (sample in names(SeuratList)) {
#   
#   # Get the Seurat object for the current sample
#   seurat_obj <- SeuratList[[sample]]
#   
#   # Filter species_df for this sample
#   df_sample <- subset(species_df, sample == sample)
#   
#   # Make sure barcodes match
#   df_sample <- df_sample[!duplicated(df_sample$barcode), ]
#   rownames(df_sample) <- df_sample$barcode
#   
#   # Match the species to the Seurat object's barcodes
#   matched_species <- df_sample[colnames(seurat_obj), "species_label"]
#   
#   # Add to metadata
#   seurat_obj$species <- matched_species
#   
#   # Store back in list
#   SeuratList[[sample]] <- seurat_obj
# }
# 
# for (i in seq_along(SeuratList)) {
#   
#   seurat_obj <- SeuratList[[i]]
#   
#   # Keep only cells from dominant species
#   cells_to_keep <- colnames(seurat_obj)[seurat_obj$species == "Human"]
#   
#   # Subset Seurat object
#   seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
#   
#   # Store it back
#   SeuratList[[i]] <- seurat_obj
# }
# 
# humanOnly_qc_fail_cell_counts_df <- data.frame(
#   sample = names(SeuratList),
#   n_cells = sapply(SeuratList, ncol)
# )
# qc_fail_cell_counts_df
# humanOnly_qc_fail_cell_counts_df
# 
# qcFailed_list <- SeuratList
# ## now merge them all into 1 object
# 
# # Ensure the list is named (this is crucial!)
# if (is.null(names(qcFailed_list)) || any(names(qcFailed_list) == "")) {
#   names(qcFailed_list) <- paste0("S", seq_along(qcFailed_list))  # or your real sample IDs
# }
# 
# # (A) Set sample metadata and prefix cell barcodes BEFORE merging
# for (nm in names(qcFailed_list)) {
#   qcFailed_list[[nm]]$sample <- nm
#   qcFailed_list[[nm]] <- RenameCells(
#     qcFailed_list[[nm]],
#     new.names = paste0(nm, "_", Cells(qcFailed_list[[nm]]))  # forces prefix
#   )
# }
# 
# # (B) Merge all at once (preferred pattern in Seurat)
# qry_fail <- merge(
#   x  = qcFailed_list[[1]],
#   y  = qcFailed_list[2:length(qcFailed_list)],
#   add.cell.ids = names(qcFailed_list)  # harmless here since we've already prefixed
# )
# 
# qry_fail$qc_status <- "QC_fail"
# Idents(qry_fail) <- "QC_fail"
# 
# # sanity check
# head(colnames(qry_fail))         # should look like "S1_ATGCCAGCAACGTTCC-1"
# table(qry_fail$sample)
# 
# ## just to make sure
# DefaultAssay(VU40T.combined) <- "RNA"
# 
# DefaultAssay(qry_fail) <- "RNA"
# qry_fail <- NormalizeData(qry_fail, verbose = FALSE)
# qry_fail <- FindVariableFeatures(qry_fail, verbose = FALSE)
# 
# npcs_ref <- ncol(Embeddings(VU40T.combined, "pca"))
# dims_use <- 1:10  # e.g., 1:10
# 
# ## store model for ref
# VU40T.combined <- RunUMAP(
#   VU40T.combined,
#   reduction = "pca",
#   dims = 1:10,
#   return.model = TRUE        # <-- critical
#   # umap.method = "uwot"     # default is fine; model is saved with uwot too
#   # metric = "cosine"
# )
# 
# VU40T.combined <- RunUMAP(
#   VU40T.combined,
#   reduction = "pca",
#   dims = 1:10,         # same number of PCs you used before!
#   return.model = TRUE  # <--- this is critical
# )
# 
# # 2) Find anchors and map QC-fails onto ref UMAP
# anchors <- FindTransferAnchors(
#   reference = VU40T.combined,
#   query = qry_fail,
#   dims = 1:10,
#   normalization.method = "LogNormalize",
#   reference.reduction = "pca"
# )
# 
# qry_fail <- MapQuery(
#   anchorset = anchors,
#   query = qry_fail,
#   reference = VU40T.combined,
#   refdata = list(ref_cluster = Idents(VU40T.combined)),   # optional: predicted cluster labels
#   reference.reduction = "pca",
#   reduction.model = "umap"                      # project into ref’s UMAP space
# )
# 
# # ensure qry_fail has PCA first
# qry_fail <- RunPCA(qry_fail, npcs = 1:10)
# 
# qry_fail <- ProjectUMAP(
#   object      = qry_fail,
#   umap.model  = ref[["umap"]]@misc$model,
#   reduction   = "pca"            # use the PCA you just computed
# )
# 
# 
# qry_fail$qc_status <- "QC_fail"
# 
# 
# merged <- merge(VU40T.combined, qry_fail)
# 
# # highlight vector
# fail_cells <- Cells(merged)[merged$qc_status == "QC_fail"]
# 
# # clusters as colors + black overlay for QC-fails
# DimPlot(
#   merged, reduction = "umap",
#   group.by = "seurat_clusters",
#   label = TRUE, repel = TRUE, shuffle = TRUE
# ) + DimPlot(
#   merged, reduction = "umap",
#   cells = fail_cells,
#   cols = "black", pt.size = 1.6
# )

# Keep ref clusters for pass cells; give fail cells their own Ident
Idents(merged) <- ifelse(merged$qc_status == "QC_fail",
                         "QC_fail",
                         as.character(Idents(merged)))

# Plots
DimPlot(merged, reduction = "umap", group.by = "qc_status")
DimPlot(merged, reduction = "umap", label = TRUE, repel = TRUE)

# Optional: per-sample view (if you stored sample)
if ("sample" %in% colnames(merged@meta.data)) {
  DimPlot(merged, reduction = "umap", group.by = "sample", shuffle = TRUE)
}

p_ref  <- DimPlot(VU40T.combined, reduction="umap",
                  group.by="seurat_clusters", label=TRUE, repel=TRUE, pt.size=0.4, cols=NULL)

p_fail <- DimPlot(qry_fail, reduction="ref.umap",
                  cells=Cells(qry_fail), cols= viridisLite::viridis(n=9), pt.size=0.4)
png(filename = file.path(plots_dir, "human_onlyQC_fail_and_processedSeurat.png"), width = 12, height = 8, units = "in", res = 300)

print(p_ref + p_fail)   # patchwork overlay
dev.off()
