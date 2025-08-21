### pseudotime monocle3 w epithelial data

library(Seurat)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
library(parallely)

# --- determine number of cores ---
available <- parallelly::availableCores()
cores_to_use <- max(1, available - 1)
message("🧠 Using ", cores_to_use, " core(s) for graph_test().")


##functs
make_marker_dotplots <- function(seurat_obj, genelist, outPrefix, export = T, pseudo_order = FALSE){
  colnames(genelist) <- gsub("[._]", " ", colnames(genelist))
  ## revert to whitespaces if whitespaces in original colname or change underscores to whitespace

  if (pseudo_order == TRUE) {
    ## Order clusters in pseudotime:
    cluster_order <- seurat_obj@meta.data %>%
      dplyr::group_by(seurat_clusters) %>%
      dplyr::summarise(avg_pseudotime = mean(pseudotime, na.rm = TRUE)) %>%
      dplyr::arrange(avg_pseudotime) %>%
      dplyr::pull(seurat_clusters)
    
    # Relevel the factor according to pseudotime
    seurat_obj$seurat_clusters <- factor(seurat_obj$seurat_clusters, levels = cluster_order)
    
    # Set output prefix for pseudotime plots
      marker_plot_dir_prefix <- "Pseudotime_humanOnly/HumanONLY_Pseudotime_"
    
  } else {
    # Not using pseudotime ordering, set standard output prefix
      marker_plot_dir_prefix <- "MarkerDotplots/HumanONLY/Res0.3/HumanONLY_"
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
    if (colnames(genelist)[i] != "Maturation Trajectory TFs"){
      # Generate dot plot
      p <- DotPlot(seurat_obj, features = gene_set, group.by = "seurat_clusters") +
        RotatedAxis() +
        ggplot2::ggtitle(colnames(genelist)[i]) +
        scale_color_gradient(low = "lightgrey", high = "red") +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1)) + 
        labs(x = "Genes", y = "Clusters")  # Relabel axes
    } else {
      p <- DotPlot(seurat_obj, features = gene_set, group.by = "seurat_clusters") +
        coord_flip() +
        ggplot2::ggtitle(colnames(genelist)[i]) +
        scale_color_gradient(low = "lightgrey", high = "red") +
        labs(x = "Genes", y = "Clusters")  # Relabel axes
    }

    # Export or print
    if (export) {
      
      num_genes <- length(na.omit(gene_set[[i]]))
      if (num_genes <= 6) {
        plot_width <- 6  # Small gene sets: fixed small width
      } else {
        plot_width <- min(max(2 * num_genes + 4, 10), 50)
      }
      colnames(genelist) <- gsub(" ", "_", colnames(genelist))
        if(colnames(genelist)[i] != "Maturation Trajectory TFs"){
          png(file = file.path(
            git_dir,
            paste0("Integrated/Plots/", marker_plot_dir_prefix ,"Dotplot_clusterResolution0.3_",
                   outPrefix, "_", colnames(genelist)[i], "_VU40T_combined.png")),
            width = plot_width, height = 4, units = "in", res = 300)
        } else { 
          ## axis flipped for maturation trajectory TFs so dims should be flipped
          png(file = file.path(
            git_dir,
            paste0("Integrated/Plots/", marker_plot_dir_prefix ,"Dotplot_clusterResolution0.3_",
                   outPrefix, "_", colnames(genelist)[i], "_VU40T_combined.png")),
            height = plot_width, width = 4, units = "in", res = 300)
        }


      colnames(genelist) <- gsub("[._]", " ", colnames(genelist))
      print(p)
      dev.off()
    } else {
      print(p)
    }
  }
}


#### make pseudotime ridgeplots

PseudotimeRidgePlot <- function(seurat_obj, 
                                features, 
                                pseudotime_col = "pseudotime", 
                                cluster_col = "seurat_clusters", 
                                ...) {
  # Check if required columns exist
  if (!(pseudotime_col %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Column", pseudotime_col, "not found in Seurat metadata."))
  }
  if (!(cluster_col %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Column", cluster_col, "not found in Seurat metadata."))
  }
  
  # Get cluster order by average pseudotime
  cluster_order <- seurat_obj@meta.data %>%
    group_by(!!sym(cluster_col)) %>%
    summarize(mean_pt = mean(!!sym(pseudotime_col), na.rm = TRUE)) %>%
    arrange(mean_pt) %>%
    pull(!!sym(cluster_col))
  
  # Apply ordering to cluster factor
  seurat_obj@meta.data[[cluster_col]] <- factor(
    seurat_obj@meta.data[[cluster_col]],
    levels = cluster_order
  )
  
  # Generate ridgeplot
  RidgePlot(seurat_obj,
            features = features,
            group.by = cluster_col,
            ...)
}



git_dir <- "~/OralMucosa/VU40T_analysis"


## to avoid hitting that 20gb quota on my home dir,
## --- directories ---


proj_dir <- "/rds/projects/g/gendood-3dmucosa/"
cache_dir <- file.path(proj_dir, "rds_cache")
plot_dir <- file.path(git_dir, "Integrated/Plots/Pseudotime_humanOnly")
tbls_dir <- file.path(git_dir, "Integrated/Pseudotime_humanOnly")

## dir check
if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(tbls_dir)) dir.create(tbls_dir, recursive = TRUE)


## load in data
message("📦 Loading Seurat object...")
rds_path <- file.path(cache_dir, "VU40T_combined_joined_0.3_res_humanOnly.rds")
stopifnot(file.exists(rds_path))
epis <- readRDS(rds_path)

## --- convert to Monocle3 CDS ---
message("🔄 Converting Seurat object to CDS...")
cds <- as.cell_data_set(epis)
## get genenames
fData(cds)$gene_short_name <- rownames(fData(cds))

recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
# recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
list.cluster <- epis@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- epis@reductions$umap@cell.embeddings



## --- clustering and trajectory ---
message("🧬 Clustering cells and learning trajectory...")
cds <- learn_graph(cds, use_partition = F)

png(filename = file.path(plot_dir, "Epi_pseudoTime_trajectory_non_ordered.png"),width = 8, height = 4, units = "in", res = 300)
p1 <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)
print(p1)
dev.off()

### order in pseudo time

message("⏳ Ordering cells in pseudotime...")
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == 7]))

p1 <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)

p2 <- plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F, group_label_size = 5)
png(filename = file.path(plot_dir, "Epi_pseudoTime_trajectory_ordered.png"),width = 8, height = 4, units = "in", res = 300)
print(p1+p2)
dev.off()

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

# ggplot(data.pseudo, aes(monocle3_pseudotime, seurat_clusters, fill = seurat_clusters)) + geom_boxplot()

png(filename = file.path(plot_dir, "Epi_pseudoTime_trajectory_boxplot_ordered.png"),width = 8, height = 4, units = "in", res = 300)
p1 <- ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) + 
    geom_boxplot() + 
    xlab("Monocle3 Pseudotime") + 
    ylab("Seurat Clusters (Reordered by Pseudotime)")
print(p1)
dev.off()

### save back to seurat obj
epis$pseudotime <- pseudotime(cds)

### gene expression change as func of pseudotime

## --- graph_test with error handling ---
## don't run interactively on RStudio
message("📊 Running graph_test...")
deg <- tryCatch({
  graph_test(cds, neighbor_graph = "principal_graph", cores = cores_to_use)
}, error = function(e) {
  message("❌ graph_test() failed: ", e$message)
  return(NULL)
})

if (is.null(deg)) {
  warning("⚠️ graph_test() failed. Skipping DEG saving.")
} else {
  passqc_deg <- deg %>% arrange(q_value) %>% filter(status == "OK")
  write.csv(passqc_deg, file = file.path(tbls_dir, "epithelial_PASS_QC_pseudotime_degs.csv"), row.names = F)
  saveRDS(epis, file = file.path(cache_dir, "VU40T_epithelial_withpseudotime.RDS"))
  message("✅ graph_test results saved successfully.")
}


epis <- readRDS(file.path(cache_dir, "VU40T_epithelial_withpseudotime.RDS"))
pseudotime_marker_list <- as.vector(t(read.csv(file.path(git_dir, "Integrated/pseudotime_genes_umap.csv"),header = F)))

## remove empty strs from the transposition
pseudotime_marker_list <- pseudotime_marker_list[nzchar(pseudotime_marker_list)]

png(filename = file.path(plot_dir, "Epi_pseudoTime_trajectory_umaps_TFs_and_KRTs.png"),width = 10, height = 10, units = "in", res = 300)
p1 <- FeaturePlot(epis, features = pseudotime_marker_list, label = T )
print(p1)
dev.off()
FeaturePlot(epis, features = c("SOX4", "RUNX1", "ELF1", "GRHL2", "PITX1", "TP63", "CEBPB" ),label = T)

RidgePlot(epis, features = c("SOX4", "RUNX1", "ELF1", "GRHL2", "PITX1", "TP63", "CEBPB" ), sort = T)

## make ridgeplots

# Assume pseudotime is stored in metadata

# Calculate average pseudotime per cluster
cluster_order <- epis@meta.data %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarize(mean_pt = mean(pseudotime, na.rm = TRUE)) %>%
  dplyr::arrange(mean_pt) %>%
  dplyr::pull(seurat_clusters)


png(filename = file.path(plot_dir, "Epi_pseudoTime_trajectory_Ridgeplots_TFs_and_KRTs.png"),width = 10, height = 16, units = "in", res = 300)
p1 <- PseudotimeRidgePlot(epis, features = pseudotime_marker_list, combine = T )
print(p1)
dev.off()


FeaturePlot(epis, features = "pseudotime", label = T) + patchwork::plot_annotation(title = "Epithelial Pseudotime") 

# write.csv(table(epis$seurat_clusters), file = file.path(git_dir, "num_cells_human_only_epi_clust.csv"))

epi_clust <-data.frame(table(epis$seurat_clusters))
p <- ggplot(epi_clust, aes(x = Var1, y = Freq)) +
      geom_bar(stat = "identity") +
      labs(x = "Cluster", y = "Cell Count", title = "Cell Counts per Epithelial Cluster") +
      theme_minimal()
png(file.path(git_dir, "Integrated/Plots/number_cells_per_clust_humanOnlyEpi.png"), width = 8, height = 4, units = "in", res = 300)
print(p)
dev.off()

### trajectory dotplots

marker_list_dp <- read.csv(file.path(git_dir, "Integrated/TF_trajectory_dotplot_markers.csv"), encoding = "UTF-8", check.names = FALSE, stringsAsFactors = FALSE)


make_marker_dotplots(epis, marker_list_dp, export = T, outPrefix = "TF_trajectory_markers", pseudo_order = TRUE)
marker_list_dp <- read.csv("~/marker_list_3.csv")
make_marker_dotplots(epis, marker_list_dp, export = T, outPrefix = "IL_WNT_COX_OXstress_trajectory_markers")