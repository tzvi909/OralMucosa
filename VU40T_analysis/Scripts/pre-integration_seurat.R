#### libs

library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)
library(SingleCellExperiment)


git_dir <- "~/OralMucosa/VU40T_analysis"

###functs
optimize_kparam <- function(seurat_obj, 
                            dims = 1:20, 
                            k_values = c(10, 20, 30, 50), 
                            resolution = 0.1, 
                            assay = "RNA", 
                            slot = "data",
                            do_plot = TRUE) {
  
  require(Seurat)
  require(cluster)
  require(dplyr)
  require(ggplot2)
  
  results <- list()
  metrics <- data.frame(k.param = integer(), 
                        n_clusters = integer(), 
                        avg_silhouette = numeric())
  
  for (k in k_values) {
    message("Testing k.param = ", k)
    
    # Copy to avoid modifying original object
    obj_copy <- seurat_obj
    
    # Run neighborhood graph and clustering
    obj_copy <- FindNeighbors(obj_copy, dims = dims, k.param = k, verbose = FALSE)
    obj_copy <- FindClusters(obj_copy, resolution = resolution, verbose = FALSE)
    
    clusters <- obj_copy$seurat_clusters
    results[[paste0("k", k)]] <- clusters
    
    # Calculate silhouette width
    pca <- Embeddings(obj_copy, "pca")[, dims]
    sil <- cluster::silhouette(as.integer(clusters), dist(pca))
    avg_sil <- mean(sil[, 3])
    
    # Store metrics
    metrics <- rbind(metrics, data.frame(k.param = k, 
                                         n_clusters = length(unique(clusters)), 
                                         avg_silhouette = avg_sil))
  }
  
  if (do_plot) {
    p1 <- ggplot(metrics, aes(x = k.param, y = n_clusters)) + 
      geom_line() + geom_point() +
      labs(title = "Number of Clusters vs k.param", y = "Number of Clusters", x = "k.param")
    
    p2 <- ggplot(metrics, aes(x = k.param, y = avg_silhouette)) + 
      geom_line() + geom_point() +
      labs(title = "Average Silhouette vs k.param", y = "Avg. Silhouette Width", x = "k.param")
    
    print(p1)
    print(p2)
  }
  
  return(list(clusterings = results, metrics = metrics))
}





optimize_resolution_by_silhouette <- function(seurat_obj, 
                                              resolutions = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0),
                                              dims = 1:min_pc,
                                              verbose = TRUE) {
  require(Seurat)
  require(cluster)
  require(dplyr)
  require(ggplot2)
  
  # Make a copy to avoid overwriting the original object
  sobj <- seurat_obj
  
  # Build SNN graph
  sobj <- FindNeighbors(sobj, dims = dims, verbose = verbose)
  
  # Run clustering at multiple resolutions
  for (res in resolutions) {
    sobj <- FindClusters(
      sobj, resolution = res, verbose = verbose,
      algorithm = 1,               # optional, can be changed
      group.singletons = FALSE    # optional
    )
    # Rename the active identity column
    cluster_col <- paste0("res_", res)
    sobj[[cluster_col]] <- sobj@meta.data$seurat_clusters
  }
  
  # Function to compute silhouette from clustering column
  compute_silhouette <- function(obj, cluster_col) {
    pca_mat <- Embeddings(obj, "pca")[, dims]
    clusters <- obj[[cluster_col]][,1]
    if (length(unique(clusters)) < 2) return(NA)
    sil <- cluster::silhouette(as.integer(as.factor(clusters)), dist(pca_mat))
    mean(sil[, 3])
  }
  
  # Compute silhouette score for each resolution
  sil_scores <- sapply(resolutions, function(res) {
    col_name <- paste0("RNA_snn_res.", res)
    compute_silhouette(sobj, col_name)
  })
  
  sil_df <- data.frame(resolution = resolutions, silhouette = sil_scores)
  
  # Generate plot
  p <- ggplot(sil_df, aes(x = resolution, y = silhouette)) +
    geom_line() + geom_point(size = 2) +
    theme_minimal() +
    labs(title = "Silhouette Score vs Resolution",
         x = "Resolution", y = "Avg Silhouette Width")
  
  print(p)
  
  # Choose best resolution
  best_res <- sil_df$resolution[which.max(sil_df$silhouette)]
  if (verbose) {
    cat("✅ Best resolution by silhouette:", best_res,
        "(score =", round(max(sil_df$silhouette), 3), ")\n")
  }
  
  list(
    seurat = sobj,
    silhouette_df = sil_df,
    best_resolution = best_res,
    silhouette_plot = p
  )
}



get_optimal_pc_count <- function(seurat_obj, variance_cutoff = 90, percent_drop = 0.1) {
  stdv <- seurat_obj[["pca"]]@stdev
  sum_stdv <- sum(stdv)
  percent_stdv <- (stdv / sum_stdv) * 100
  cumulative <- cumsum(percent_stdv)
  
  # First condition: capture >90% variance but avoid very flat PCs
  co1 <- which(cumulative > variance_cutoff & percent_stdv < 5)[1]
  
  # Second condition: dropoff in variance explained > 0.1%
  diffs <- percent_stdv[1:(length(percent_stdv) - 1)] - percent_stdv[2:length(percent_stdv)]
  co2 <- sort(which(diffs > percent_drop), decreasing = TRUE)[1] + 1
  
  # Fallback if one of them is NA
  if (is.na(co1)) co1 <- Inf
  if (is.na(co2)) co2 <- Inf
  
  min_pc <- min(co1, co2)
  
  return(min_pc)
}



###main
seurat_singlets_list <- readRDS("VU40T_singlets_only_sep_samples_Mouse.RDS")



pc_summary <- data.frame(sample = character(), optimal_pcs = integer(), stringsAsFactors = FALSE)

for (sample_name in names(seurat_singlets_list)) {
  seurat_obj <- seurat_singlets_list[[sample_name]]
  optimal_pc <- get_optimal_pc_count(seurat_obj)
  
  pc_summary <- rbind(pc_summary, data.frame(
    sample = sample_name,
    optimal_pcs = optimal_pc
  ))
}


# for (i in seq_along(seurat_singlets_list)){
#   sample_name <- names(seurat_singlets_list)[i]
#   min_pc <- pc_summary$optimal_pcs[pc_summary$sample == sample_name]
#   k_param_res[[i]] <- optimize_kparam_by_silhouette(seurat_singlets_list[[i]], k_values = c(0.1, 0.25, 0.5, 0.75, 1.0), dims = 1:min_pc)
# }

k_params <- list()
k_vals <- c(10, 20, 30, 40, 50)
resolution_use <- 0.1

kparam_res <- lapply(names(seurat_singlets_list), function(sample_name) {
  message("Running on: ", sample_name)
  
  seurat_obj <- seurat_singlets_list[[sample_name]]
  
  # Get optimal PCs for this sample
  dims_use <- pc_summary$optimal_pcs[pc_summary$sample == sample_name][[1]]
  
  optimize_kparam(seurat_obj, 
                  dims = dims_use, 
                  k_values = k_vals, 
                  resolution = resolution_use,
                  do_plot = TRUE)
})
names(kparam_res) <- names(seurat_singlets_list)
all_K_metrics <- do.call(rbind, lapply(names(kparam_res), function(sample_name) {
  df <- kparam_res[[sample_name]]$metrics
  df$sample <- sample_name
  return(df)
}))

best_k_by_sample <- all_K_metrics %>%
  group_by(sample) %>%
  slice_max(order_by = avg_silhouette, n = 1, with_ties = FALSE) %>%
  ungroup()


ggplot(all_K_metrics, aes(x = k.param, y = avg_silhouette, color = sample)) +
  geom_line() + geom_point() +
  labs(title = "Silhouette Score by k.param", y = "Avg. Silhouette", x = "k.param")


# # Name each result with the corresponding sample
# names(results_list) <- names(seurat_singlets_list)

best_global_k <- all_K_metrics %>%
  group_by(k.param) %>%
  summarise(mean_silhouette = mean(avg_silhouette)) %>%
  slice_max(mean_silhouette, n = 1)

## can leave k param as default for human and mouse





resolution_res <- list()
for (i in seq_along(seurat_singlets_list)){
  sample_name <- names(seurat_singlets_list)[i]
  min_pc <- pc_summary$optimal_pcs[pc_summary$sample == sample_name]
  resolution_res[[i]] <- optimize_resolution_by_silhouette(seurat_singlets_list[[i]], resolutions = c(0.1, 0.25, 0.5, 0.75, 1.0), dims = 1:min_pc)
}



## best k param is default and best resolution is 0.1

## rerun clustering with optimised dims

for (i in seq_along(seurat_singlets_list)){
  seurat_singlets_list[[i]] <- FindNeighbors(seurat_singlets_list[[i]], dims=1:min_pc)
  seurat_singlets_list[[i]] <- FindClusters(seurat_singlets_list[[i]], resolution=0.1)
}


for (i in seq_along(seurat_singlets_list)){
  sample_name <- names(seurat_singlets_list)[i]
  
  plot1 <- DimPlot(seurat_singlets_list[[i]], reduction = "umap", label = T)
  
  plot2 <- DimPlot(seurat_singlets_list[[i]], reduction = "tsne", label = T)
  
  png(file = paste0(git_dir, "/Seperate_samples/Mouse/Plots/SingletsOnly_UMAP_and_t-sne_", sample_name, ".png"), width = 10, height = 5, units = "in", res = 300, )
  
  print(
    plot1 + plot2 +
      patchwork::plot_annotation(title = paste0("Doublet-Filtered UMAP and t-SNE: ", sample_name))
  )
  
  dev.off()
}

# optional quick annotation
bp <- BiocParallel::MulticoreParam(workers = 8)
ref <- readRDS("./R/HPCA_reference.rds")

library(SingleR)
for (i in seq_along(seurat_singlets_list)) {
  seurat_obj <- seurat_singlets_list[[i]]
  sample_name <- names(seurat_singlets_list)[i]

  # outfile <- file.path(git_dir, "Seperate_samples/SingleR_annotations", paste0(sample_name, "_SingleR_asSepCells.rds"))

  # if (file.exists(outfile)) {
  #   message(sample_name, " already annotated. Skipping.")
  #   next
  # }

  message("Annotating ", sample_name)

  # Convert to SingleCellExperiment
  sce <- as.SingleCellExperiment(seurat_obj)

  # Get cluster labels
  cluster_labels <- seurat_obj$seurat_clusters

  # Run SingleR annotation
  pred <- SingleR::SingleR(
    test = sce,
    ref = ref,
    labels = ref$label.fine,
    clusters =  cluster_labels,
    #assay.type.test = 1,
    BPPARAM=bp
  )

  # Save prediction
  # saveRDS(pred, outfile)

  # Apply cluster annotations to each cell
  seurat_obj[["SingleR_cluster_label"]] <- pred$labels[
    match(seurat_obj$seurat_clusters, rownames(pred))
  ]

  ## per cell annotations

  # # Check alignment (rownames(pred) should match colnames(seurat_obj))
  # if (!all(rownames(pred) == colnames(seurat_obj))) {
  #   stop("Cell names in SingleR result do not match Seurat object.")
  # }
  #
  #
  # # (Optional) assign pruned labels or other metadata
  # seurat_obj$SingleR_pruned <- pred$pruned.labels
  #
  # # # Apply cluster annotations to each cell
  # seurat_obj[["SingleR_label"]] <- pred$labels


  # Store back into list
  seurat_singlets_list[[i]] <- seurat_obj
}

##checkpoint
saveRDS(seurat_singlets_list, "annotated_VU40T_sepsamples_singlets_only_mouse.rds")

for (i in seq_along(seurat_singlets_list)) {
  sample_name <- names(seurat_singlets_list)[i]
  seurat_obj <- seurat_singlets_list[[i]]

  p1 <- DimPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE, repel = TRUE)

  p2 <- DimPlot(seurat_obj, group.by = "SingleR_cluster_label", label = F, repel = TRUE)
  png(file = file.path(git_dir, paste0("/Seperate_samples/Mouse/Plots/Seurat_clusters_vs_SingleR_UMAP_ascluster_", sample_name, ".png")), 
      width = 16, height = 5, units = "in", res = 300)
  print(p1 + p2 + patchwork::plot_annotation(title = sample_name))
  dev.off()
}
# output_dir <- file.path(git_dir, "Seperate_samples/Plots")

SeuratList_markers <- list()
for (i in seq_along(seurat_singlets_list)){
  SeuratList_markers[[i]] <- FindAllMarkers(seurat_singlets_list[[i]],
                                            only.pos = T,
                                            min.pct = 0.25,
                                            logfc.threshold = 0.25)
}

# Set names to match SeuratList for easy referencing


names(SeuratList_markers) <- names(seurat_singlets_list)


for (i in seq_along(SeuratList_markers)) {
  # Extract sample name and corresponding marker table
  sample_name <- names(SeuratList_markers)[i]
  markers <- SeuratList_markers[[i]]
  seurat_obj <- seurat_singlets_list[[sample_name]]
  
  # Split markers by cluster and get top 3 per cluster
  topMarkers <- split(markers, markers$cluster)
  top10 <- lapply(topMarkers, head, n = 10)
  top10 <- do.call("rbind", top10)
  
  # Display or save FeaturePlot
  print(paste("Top 10 markers for", sample_name))
  # print(top10)
  
  
  # Generate FeaturePlots
  plots <- lapply(top10$gene, function(gene) {
    FeaturePlot(seurat_obj, features = gene) + ggtitle(gene)
  })
  
  # Combine with patchwork, arrange in 10 columns
  combined_plot <- patchwork::wrap_plots(plots, ncol = 10)
  
  # Save to PNG
  output_file <- file.path(git_dir, paste0("Seperate_samples/Plots/top10_markers_", sample_name, ".png"))
  
  # Adjust image dimensions
  n_rows <- ceiling(length(plots) / 10)
  width_px <- 10 * 300   # ~300 px per plot
  height_px <- n_rows * 350
  
  png(filename = output_file,
      width = width_px, height = height_px, res = 150)
  print(combined_plot)
  dev.off()
}

for (i in seq_along(seurat_singlets_list)) {
  sample_name <- names(SeuratList_markers)[i]
  markers <- SeuratList_markers[[i]]
  seurat_obj <- seurat_singlets_list[[sample_name]]
  
  
  # Split markers by cluster and get top 3 per cluster
  topMarkers <- split(markers, markers$cluster)
  top10 <- lapply(topMarkers, head, n = 10)
  top10 <- do.call("rbind", top10)
  
  # # Display or save FeaturePlot
  # print(paste("Top 10 markers for", comb.names[i]))
  # # print(top10)
  # 
  
  # Generate FeaturePlots
  plots <- lapply(top10$gene, function(gene) {
    FeaturePlot(seurat_obj, features = gene) + ggplot2::ggtitle(gene)
  })
  
  # Combine with patchwork, arrange in 10 columns
  combined_plot <- patchwork::wrap_plots(plots, ncol = 10)
  
  # Save to PNG
  output_file <- file.path(git_dir, paste0("/Seperate_samples/Mouse/Plots/top10Markers_", sample_name, ".png"))
  
  # Adjust image dimensions
  n_rows <- ceiling(length(plots) / 10)
  width_px <- 10 * 300   # ~300 px per plot
  height_px <- n_rows * 350
  
  png(filename = output_file,
      width = width_px, height = height_px, res = 150)
  print(combined_plot)
  dev.off()
}



for (i in seq_along(SeuratList_markers)) {
  sample_name <- names(SeuratList_markers)[i]
  markers <- SeuratList_markers[[i]][SeuratList_markers[[i]]$p_val_adj <= 0.05, ]
  
  write.csv(markers, file = file.path(git_dir, paste0("/Seperate_samples/Mouse/ClusterMarkers_", sample_name, ".csv")), row.names = FALSE)
}


