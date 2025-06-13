# Load packages
library(Seurat)
library(Matrix)
library(tidyverse)

# Helper function: compute species % per barcode
compute_species_df <- function(human_counts, mouse_counts) {
  all_barcodes <- union(names(human_counts), names(mouse_counts))
  
  df <- data.frame(
    barcode = all_barcodes,
    human_counts = human_counts[all_barcodes],
    mouse_counts = mouse_counts[all_barcodes]
  )
  
  df[is.na(df)] <- 0
  df$total_counts <- df$human_counts + df$mouse_counts
  df$percent_human <- df$human_counts / df$total_counts * 100
  df$percent_mouse <- df$mouse_counts / df$total_counts * 100
  
  return(df)
}

convert_ensembl_to_symbol <- function(seurat_obj, species = c("human", "mouse"), assay = "RNA") {
  # Load biomaRt if needed
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install("biomaRt")
  }
  library(biomaRt)
  
  # Match species and set Ensembl dataset
  species <- match.arg(species)
  dataset <- switch(
    species,
    human = "hsapiens_gene_ensembl",
    mouse = "mmusculus_gene_ensembl"
  )
  symbol_attr <- switch(
    species,
    human = "hgnc_symbol",
    mouse = "mgi_symbol"
  )
  
  # Connect to Ensembl
  message("Connecting to Ensembl for ", species, "...")
  ensembl <- useEnsembl(biomart = "genes", dataset = dataset)
  
  # Extract Ensembl IDs
  ens_ids <- rownames(seurat_obj[[assay]])
  
  # Get gene symbol mapping
  gene_map <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", symbol_attr),
    values = ens_ids,
    mart = ensembl
  )
  
  # Clean map
  gene_map <- gene_map[gene_map[[symbol_attr]] != "", ]
  gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]
  common_ids <- intersect(rownames(seurat_obj[[assay]]), gene_map$ensembl_gene_id)
  
  # Replace Ensembl IDs with symbols
  new_names <- toupper(gene_map[[symbol_attr]][match(common_ids, gene_map$ensembl_gene_id)])
  rownames(seurat_obj[[assay]])[match(common_ids, rownames(seurat_obj[[assay]]))] <- new_names
  
  message("Replaced ", length(common_ids), " Ensembl IDs with gene symbols.")
  return(seurat_obj)
}


## func to get aggregate rows w/ duplicated mouse gene names
aggregate_by_gene_fast <- function(mat) {
  # Convert sparse matrix to triplet format (i,j,x)
  m_coo <- summary(mat)
  genes <- rownames(mat)
  barcodes <- colnames(mat)
  
  # Add gene names to triplet
  m_coo$gene <- genes[m_coo$i]
  m_coo$barcode <- barcodes[m_coo$j]
  
  # Convert to data.table
  dt <- as.data.table(m_coo)
  
  # Sum duplicated gene-barcode combinations
  dt_agg <- dt[, .(x = sum(x)), by = .(gene, barcode)]
  
  # Convert back to sparse matrix
  gene_levels <- unique(dt_agg$gene)
  barcode_levels <- unique(dt_agg$barcode)
  
  sparse_mat <- Matrix::sparseMatrix(
    i = match(dt_agg$gene, gene_levels),
    j = match(dt_agg$barcode, barcode_levels),
    x = dt_agg$x,
    dimnames = list(gene_levels, barcode_levels)
  )
  
  return(sparse_mat)
}

# For each sample:
species_list <- list()
seurat_list <- list()

git_dir <- "~/OralMucosa/VU40T_analysis/"
proj_dir <- "/rds/projects/g/gendood-3dmucosa/"
plots_dir <- file.path(git_dir, "Mixed_species_plots/")
if (!(dir.exists(plots_dir))) {
  dir.create(plots_dir)
}
## wd
mtx_dir_mouse <- file.path(proj_dir, "BaseSpace/LPS_VU40T_QC_and_counts_Mouse/cellranger/mtx_conversions")
mtx_dir_human <- file.path(proj_dir, "BaseSpace/LPS_VU40T_QC_and_counts/cellranger/mtx_conversions")
Samples <- list.files(path = mtx_dir_human, 
                      pattern = "^P",
                      recursive = F
) ##check we have all samples
## get only cellbender seurat_obj paths
seurat_objs_human <- list.files(path = mtx_dir_human, 
                          pattern = "^P.*cellbender.*.seurat.rds$",
                          full.names = T,
                          recursive = T
)

seurat_objs_mouse <- list.files(path = mtx_dir_mouse, 
                                pattern = "^P.*cellbender.*.seurat.rds$",
                                full.names = T,
                                recursive = T
)


for (i in Samples) {
  # Load CellBender matrices
  seurat_human <- readRDS(paste0(mtx_dir_human, "/", i,"/" , i, "_cellbender_filter_matrix.seurat.rds"))
  seurat_mouse <- readRDS(paste0(mtx_dir_mouse, "/", i, "/", i, "_cellbender_filter_matrix.seurat.rds"))
  
  ##get total counts per barcode
  human_counts_total <- Matrix::colSums(GetAssayData(seurat_human, slot = "counts"))
  mouse_counts_total <- Matrix::colSums(GetAssayData(seurat_mouse, slot = "counts"))
  
  # Compute species score df
  species_df <- compute_species_df(human_counts_total, mouse_counts_total)
  species_df$sample <- i
  species_df$species_label <- case_when(
    species_df$percent_human >= 70 ~ "Human",
    species_df$percent_mouse >= 70 ~ "Mouse",
    TRUE ~ "Mixed"
  )
  # Store species df and Seurat object for later
  species_list[[i]] <- species_df
  
  # For visualization, pick one object to hold metadata (e.g. human Seurat object)
  barcodes_in_seurat <- colnames(seurat_human)
  species_df_sub <- species_df[species_df$barcode %in% barcodes_in_seurat, ]
  species_df_sub <- species_df_sub[match(barcodes_in_seurat, species_df_sub$barcode), ]
  
  ## set quite a liberal threshold for species mixing 
  ## due to expected ~85% sequence homology between mus musculus and homo sapiens
  seurat_human$percent_human <- species_df_sub$percent_human
  seurat_human$percent_mouse <- species_df_sub$percent_mouse
  seurat_human$species_label <- case_when(
    seurat_human$percent_human >= 70 ~ "Human",
    seurat_human$percent_mouse >= 70 ~ "Mouse",
    TRUE ~ "Mixed"
  )
  
  # Add sample ID
  seurat_human$sample <- paste0("Sample", i)
  
  seurat_list[[i]] <- seurat_human
}

## inspect species distribution
ggplot(species_list[[1]][species_list[[1]]$species_label == "Mixed", ]) + geom_histogram(aes(x = percent_human))

## output this for further analysis on seperate alignments

species_df_all <- do.call(rbind, species_list)
write.csv(species_df_all, file = file.path(
  git_dir, "VU40T_species_assignment_by_reads.csv"
)
            )

# Combine all samples into one Seurat object
combined <- merge(seurat_list[[1]], y = seurat_list[2:4])

# Visualization — UMAP (you can do more!)
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:20)

# Plot species scores

png(filename = file.path(
  plots_dir, "human_mouse_seurat_percent_species_UMAP_all_cells.png"), 
  width = 12, 
  height = 8, 
  units = "in", 
  res =  300)

p <- FeaturePlot(combined, features = c("percent_human", "percent_mouse"))
print(p)
dev.off()

# Plot species label

png(filename = file.path(
  plots_dir, "human_mouse_seurat_species_assignment_preprocessing_UMAP.png"), 
  width = 12, 
  height = 8, 
  units = "in", 
  res =  300)

p <- DimPlot(combined, group.by = "species_label", label = TRUE) + ggtitle("Species assignment per cell")
print(p)
dev.off()


combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

# Plot species label on UMAP
DimPlot(combined, group.by = "species_label",  label = TRUE) + ggtitle("Species assignment on UMAP")


# ggplot(combined@meta.data %>% filter(species_label == "Mouse" | species_label == "Mixed"), 
#        aes(x = percent_human, fill = species_label)) +
#   geom_density(alpha = 0.5) +
#   theme_classic() +
#   ggtitle("Percent human in Mouse + Mixed cells")

combined$condition <- ifelse(grepl("LPS-N", combined$sample), "LPS-N", "LPS-P")
combined$passage <- gsub("_.*", "", combined$sample)

joined_layers <- JoinLayers(combined)




markers <- FindAllMarkers(object = joined_layers, min.pct = 0.25, only.pos = T, logfc.threshold = 0.25)

sig_markers <- markers[markers$p_val_adj <= 0.05, ]


# Plot species assignment

png(filename = file.path(
  plots_dir, "human_mouse_seurat_clusters_preprocessing_UMAP.png"), 
    width = 12, 
    height = 8, 
    units = "in", 
    res =  300)

p <- DimPlot(combined, group.by = "seurat_clusters", split.by = "condition", label = TRUE)
print(p)
dev.off()


VlnPlot(combined, features = "percent_human", group.by = "seurat_clusters") # or group.by = "species_label_refined"


FeaturePlot(combined, features = c("EPCAM", "KRT8", "KRT18"))


library(msigdbr)
library(dplyr)

# human_fibro <- msigdbr(species = "Homo sapiens") %>%
#   filter(grepl("fibroblast", tolower(gs_name)))
# 
# # For mouse fibroblast markers
# mouse_fibro <- msigdbr(species = "Mus musculus") %>%
#   filter(grepl("fibroblast", tolower(gs_name)))
# 
# mouse_genes <- unique(mouse_fibro$gene_symbol)
# 
# human_epi <- msigdbr(species = "Homo sapiens") %>%
#   filter(grepl("epithelial", tolower(gs_name)))
# 
# # For mouse fibroblast markers
# mouse_epi <- msigdbr(species = "Mus musculus") %>%
#   filter(grepl("epithelial", tolower(gs_name)))
# # 
# # human_genes <- unique(mouse_fibro$gene_symbol)
# # mouse_genes <- unique(mouse_fibro$gene_symbol)
# 
# human_genes <- unique(human_epi$gene_symbol)
# mouse_genes <- unique(mouse_epi$gene_symbol)
# 
# # Make gene symbols comparable
# human_uc_epi <- toupper(human_genes)
# mouse_uc_epi <- toupper(mouse_genes)
# 
# shared_fibroblast_markers <- intersect(mouse_uc, human_uc)
# human_only_blacklist <- setdiff(mouse_uc, human_uc)
# 
# shared_epi_markers <- intersect(mouse_uc, human_uc)
# mouse_only_epi_blacklist <- setdiff(mouse_uc, human_uc)
# 
# cat("✅ Shared fibroblast markers (Human & Mouse):\n")
# print(shared_fibroblast_markers)
# 
# cat("\n🚫 Mouse-only markers (blacklist):\n")
# print(mouse_only_epi_blacklist)
# 
# markers_list <- read.csv("~/OralMucosa/VU40T_analysis/Integrated/VU40T_significant_markers_resolution_0.5.csv")
# 
# filtered_markers <- markers %>%
#   filter(gene %in% mouse_only_blacklist) ## didn't filter any fibroblast markers
# 
# 
# filtered_markers <- markers %>%
#   filter(!gene %in% mouse_only_epi_blacklist) ## filtered 9 mouse epithelial markers
# filtered_markers <- markers %>%filter(gene %in% mouse_only_epi_blacklist) ## filtered 9 mouse epithelial markers
# 
# VU40T_0.5res

