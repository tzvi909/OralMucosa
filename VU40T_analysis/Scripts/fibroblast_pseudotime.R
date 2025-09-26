### pseudotime monocle3 w fibroblast data

library(Seurat)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
library(biomaRt)
library(parallelly)
library(argparse)


# ── Get Args ─────────────────────────────────────────────────────

parser <- ArgumentParser(description='Run Pseudotime Analysis on Mouse Fibroblast Cells')
parser$add_argument('-r', '--root_cluster', default = 9, type="integer",
                    help='set root cluster num for analysis: \n
                    9 for proliferating cells, 0 for undifferentiated cells')

argv <- parser$parse_args()
# ── Core detection ───────────────────────────────────────────────
cores_to_use <- max(1, parallelly::availableCores() - 1)
message("🧠 Using ", cores_to_use, " core(s) for graph_test().")

# ── Paths ────────────────────────────────────────────────────────
git_dir <- "~/OralMucosa/VU40T_analysis"
proj_dir <- "/rds/projects/g/gendood-3dmucosa/"
cache_dir <- file.path(proj_dir, "rds_cache")
plot_dir <- file.path(git_dir, "Integrated/Mouse/Plots/Pseudotime_MouseOnly")

dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load Data ─────────────────────────────────────────────────────
VU40T.combined <- readRDS(file.path(cache_dir, "VU40T_combined_joined_0.6_res_mouseOnly.rds"))

# Get rid of contaminating Epithelial cluster 13
# though it has no bearing on the pseudotime analysis as it is completely seperate from all other clusters

VU40T.combined <- subset(VU40T.combined, idents = "13", invert = T)


cds <- as.cell_data_set(VU40T.combined)
fData(cds)$gene_short_name <- rownames(fData(cds))

# ── Add clustering and UMAP info ─────────────────────────────────
recreate.partitions <- factor(rep(1, ncol(cds)), levels = 1)
names(recreate.partitions) <- colnames(cds)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
cds@clusters@listData[["UMAP"]][["clusters"]] <- VU40T.combined@active.ident
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- VU40T.combined@reductions$umap@cell.embeddings

# ── Trajectory: Learn graph ──────────────────────────────────────
root_cluster <- argv$root_cluster
cds <- learn_graph(cds, use_partition = FALSE)

# ── Plot trajectory (before ordering) ────────────────────────────
png(file.path(plot_dir, 
                paste0("Fibro_pseudoTime_root_clust_", root_cluster , "_trajectory_non_ordered.png")
                ), width = 8, height = 4, units = "in", res = 300)
p1 <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = FALSE,
                 label_branch_points = F, label_roots = TRUE, label_leaves = FALSE,
                 group_label_size = 5)
print(p1)
dev.off()

# ── Root cluster auto-detection ──────────────────────────────────
# Get all Seurat clusters

# Order pseudotime using those cells
message("Root cluster set to: ", root_cluster, "for this pseudotime analysis")
root_cells <- colnames(cds[, clusters(cds) == root_cluster])


cds <- order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)


# ── Plot ordered trajectory ──────────────────────────────────────
p1 <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = F, label_roots = F, label_leaves = F,
           group_label_size = 5, label_principal_points = F, trajectory_graph_segment_size = 2)

p2 <- plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F, group_label_size = 5, 
           label_principal_points = F, trajectory_graph_segment_size = 2)

png(file.path(plot_dir, 
    paste0("Fibro_pseudoTime_root_clust_", root_cluster ,"_trajectory_ordered.png")
    ), width = 8, height = 4, units = "in", res = 300)
print(p1 + p2)
dev.off()

# ── Boxplot of pseudotime vs cluster ─────────────────────────────
cds$monocle3_pseudotime <- pseudotime(cds)
VU40T.combined$pseudotime <- cds$monocle3_pseudotime
data.pseudo <- as.data.frame(colData(cds))

png(file.path(plot_dir, 
    paste0("Fibro_pseudoTime_trajectory_root_clust_", root_cluster ,"_boxplot_ordered.png")
    ), width = 8, height = 4, units = "in", res = 300)
p_box <- ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) +
  geom_boxplot() +
  xlab("Monocle3 Pseudotime") +
  ylab("Seurat Clusters (Reordered by Pseudotime)") +
  theme_minimal()
print(p_box)
dev.off()

# ── DEG Testing ──────────────────────────────────────────────────
## Do not run in RStudio interactively - it will take ~ 1hr
message("📊 Running graph_test...")
deg <- tryCatch({
  graph_test(cds, neighbor_graph = "principal_graph", cores = cores_to_use)
}, error = function(e) {
  message("❌ graph_test() failed: ", e$message)
  return(NULL)
})

if (!is.null(deg)) {
  passqc_deg <- deg %>% arrange(q_value) %>% filter(status == "OK")
  write.csv(passqc_deg, file = file.path(plot_dir, paste0("fibroblast_root_clust_", root_cluster ,"_PASS_QC_pseudotime_degs.csv")), row.names=FALSE)
  message("✅ DEG results saved.")
} else {
  warning("⚠️ graph_test() failed — skipping DEG save.")
}

# ── Save final Seurat object ─────────────────────────────────────
saveRDS(VU40T.combined, file = file.path(cache_dir, paste0("VU40T_Fibroblast_withpseudotime_rootclust_", root_cluster, ".RDS")))
message("💾 Seurat object with pseudotime saved.")

# # Get top 5 significant genes with 'OK' status
# top_genes <- deg %>%
#   arrange(q_value) %>%
#   filter(status == "OK") %>%
#   head(5) %>%
#   pull(gene_short_name)  # or use `row.names(.)` if gene names are rownames
# 
# # Plot them
# FeaturePlot(epis, features = top_genes, label = T)
# 
# pseudotime_marker_list <- as.vector(t(read.csv(file.path(git_dir, "Integrated/pseudotime_genes_umap.csv"),header = F)))
# 
# ## remove empty strs from the transposition
# pseudotime_marker_list <- pseudotime_marker_list[nzchar(pseudotime_marker_list)]
# 
# png(filename = file.path(plot_dir, "Epi_pseudoTime_trajectory_umaps_TFs_and_KRTs.png"),width = 10, height = 10, units = "in", res = 300)
# p1 <- FeaturePlot(epis, features = pseudotime_marker_list, label = T )
# print(p1)
# dev.off()
# FeaturePlot(epis, features = c("SOX4", "RUNX1", "ELF1", "GRHL2", "PITX1", "TP63", "CEBPB" ),label = T)
# 
# 
# FeaturePlot(epis, features = "pseudotime", label = T) + patchwork::plot_annotation(title = "Epithelial Pseudotime") 
# 
# # write.csv(table(epis$seurat_clusters), file = file.path(git_dir, "num_cells_human_only_epi_clust.csv"))
# 
# epi_clust <-data.frame(table(epis$seurat_clusters))
# p <- ggplot(epi_clust, aes(x = Var1, y = Freq)) +
#       geom_bar(stat = "identity") +
#       labs(x = "Cluster", y = "Cell Count", title = "Cell Counts per Epithelial Cluster") +
#       theme_minimal()
# png(file.path(git_dir, "Integrated/Plots/number_cells_per_clust_humanOnlyEpi.png"), width = 8, height = 4, units = "in", res = 300)
# print(p)
# dev.off()
# 
# ### trajectory dotplots
# # 
# # marker_list_dp <- read.csv(file.path(git_dir, "Integrated/TF_trajectory_dotplot_markers.csv"))
# # 
# # 
# # make_marker_dotplots(epis, marker_list_dp, export = T, outPrefix = "TF_trajectory_markers", pseudo_order = TRUE)
# # marker_list_dp <- read.csv("~/marker_list_3.csv")
# # make_marker_dotplots(epis, marker_list_dp, export = T, outPrefix = "IL_WNT_COX_OXstress_trajectory_markers")
