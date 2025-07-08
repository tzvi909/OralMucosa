### pseudotime monocle3 w epithelial data

library(Seurat)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)


git_dir <- "~/OralMucosa/VU40T_analysis"


## to avoid hitting that 20gb quota on my home dir,
proj_dir <- "/rds/projects/g/gendood-3dmucosa/"
cache_dir <- file.path(proj_dir, "rds_cache")
plot_dir <- file.path(git_dir, "Integrated/Plots/Pseudotime_humanOnly")
## check for cache dir
if(!(dir.exists(cache_dir))){
  dir.create(cache_dir, recursive = T)
}
if(!(dir.exists(plot_dir))){
  dir.create(plot_dir, recursive = T)
}

## load in data

epis <- readRDS(file.path(cache_dir, "VU40T_combined_joined_0.3_res_humanOnly.rds"))

cds <- as.cell_data_set(epis)
## get genenames
fData(cds)$gene_short_name <- rownames(fData(cds))

recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
list.cluster <- epis@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- epis@reductions$umap@cell.embeddings


cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
                                 group_label_size = 5) + theme(legend.position = "right")
cluster.before.traj


## actual pseudotime bit

cds <- learn_graph(cds, use_partition = F)
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)


### order in pseudo time

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

ggplot(data.pseudo, aes(monocle3_pseudotime, seurat_clusters, fill = seurat_clusters)) + geom_boxplot()

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot()


### gene expression change as func of pseudotime

deg <- graph_test(cds, neighbor_graph = "principal_graph")

# Get top 5 significant genes with 'OK' status
top_genes <- deg %>%
  arrange(q_value) %>%
  filter(status == "OK") %>%
  head(5) %>%
  pull(gene_short_name)  # or use `row.names(.)` if gene names are rownames

# Plot them
FeaturePlot(epis, features = top_genes, label = T)

FeaturePlot(epis, features = c("SOX4", "RUNX1", "ELF1", "GRHL2", "PITX1", "TP63", "CEBPB" ),label = T)

epis$pseudotime <- pseudotime(cds)
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

