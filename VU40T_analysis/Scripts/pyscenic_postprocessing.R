library(Seurat)
library(ggplot2)
library(dplyr)

#### dirs setup again




## to avoid hitting that 20gb quota on my home dir,
proj_dir <- "/rds/projects/g/gendood-3dmucosa/"
analysis_dir <- file.path(proj_dir, "scRNAseqAnalysis/")
git_dir <- file.path(analysis_dir, "OralMucosa/VU40T_analysis")
cache_dir <- file.path(proj_dir, "rds_cache")

## check for dirs
if(!(dir.exists(cache_dir))){
  dir.create(cache_dir, recursive = T)
}

if(!(dir.exists(analysis_dir))){
  dir.create(analysis_dir, recursive = T)
}

if(!(dir.exists(git_dir))){
  dir.create(git_dir, recursive = T)
}
### Epi Seurat Obj
epi <- readRDS(file.path(cache_dir, "VU40T_combined_joined_0.3_res_humanOnly.rds"))

### pyscenic aucell mtx
auc_fname <- file.path(analysis_dir, "pyscenic_res/human_epi_auc_mtx.csv")

auc <- data.table::fread(auc_fname,
                check.names = FALSE, index = "Cell")
data.table::setnames(auc, 1, "Cell")
auc <- as.data.frame(auc)
rownames(auc) <- auc[[1]]
auc[[1]] <- NULL
### align auc to seurat obj

common <- intersect(rownames(auc), colnames(epi))
epi <- subset(epi, cells = common)
mat <- as(as.matrix(t(auc[common, , drop = FALSE])), "dgCMatrix")  

epi[["SCENIC"]] <- CreateAssayObject(
  counts   = mat
)
DefaultAssay(epi) <- "SCENIC"

epi <- RunUMAP(VU40T.combined, reduction = "SCENIC", dims = 1:optimal_pc, verbose = F,  seed.use = 666)
DimPlot(epi, reduction = "SCENIC")

head(rownames(epi[["SCENIC"]]))       # regulon names
# VlnPlot(epi, features = head(rownames(epi[["SCENIC"]]), 3), alpha = 0)
# FeaturePlot(epi, features = c("TP63(+)", "SOX2(+)"))      # adjust names to yours

## binarise regulons
m <- GetAssayData(epi, slot = "data")  # or 'counts' (same here)
# per-regulon threshold at, e.g., 95th percentile
thr_mat <- as.data.frame(matrixStats::rowQuantiles(as.matrix(m), probs = seq(0,1,0.01), na.rm = TRUE))

## 1) Name the quantile columns for easy picking
colnames(thr_mat) <- paste0("q", sprintf("%02d", 0:100))  # q00, q01, ..., q100

## 2) Choose a percentile and extract a per-row threshold vector
p <- 95  # e.g., 95th percentile
thr_vec <- thr_mat[[ paste0("q", sprintf("%02d", p)) ]]   # numeric, length = nrow(m)
names(thr_vec) <- rownames(m)

## 3) Build a mask of values ≥ per-row threshold (sparse-friendly)
mask <- Matrix::Matrix( t( t(m) >= thr_vec ), sparse = TRUE )  # same shape as m; TRUE where keep

## (a) Filter values (zero-out below threshold), keep matrix numeric
m_above <- m * mask   # element-wise multiply; below-threshold → 0

## (b) Keep only rows that have at least one cell ≥ threshold
rows_keep <- Matrix::rowSums(mask) > 0
m_rows_kept <- m[rows_keep, , drop = FALSE]


nz_q <- function(x, p = 0.95) { x <- x[x > 0]; if (length(x)) unname(quantile(x, p)) else Inf }
thr_vec <- vapply(seq_len(nrow(m)), \(i) nz_q(as.numeric(m[i, ]), 0.95), numeric(1))
mask <- Matrix( t( t(m) >= thr_vec ), sparse = TRUE )nz_q <- function(x, p = 0.95) { x <- x[x > 0]; if (length(x)) unname(quantile(x, p)) else Inf }
thr_vec <- vapply(seq_len(nrow(m)), \(i) nz_q(as.numeric(m[i, ]), 0.95), numeric(1))
mask <- Matrix::Matrix( t( t(m) >= thr_vec ), sparse = TRUE )

thr <- setNames(as.numeric(thr), rownames(m))         # named numeric vector

# Perform quantile normalization
qn <- function(.data){
  data_sort <- apply(.data, 2, sort)
  row_means <- rowMeans(data_sort)
  data_sort <- matrix(row_means, 
                      nrow = nrow(data_sort), 
                      ncol = ncol(data_sort), 
                      byrow = TRUE
  )
  index_rank <- apply(.data, 2, order)
  normalized_data <- matrix(nrow = nrow(.data), ncol = ncol(.data))
  for(i in 1:ncol(.data)){
    normalized_data[,i] <- data_sort[index_rank[,i], i]
  }
  return(normalized_data)
}

normalized_data <- qn(m)

df <- data.frame(
  barcode = names(epi$seurat_clusters),
  cluster = as.character(epi$seurat_clusters),   # keep as character; or leave as factor if you want
  stringsAsFactors = FALSE
)

write.csv(df, file = file.path(git_dir, "Integrated/epi_seurat_barcode_annotation.csv"))














bin <- sweep(m, 1, thr, FUN = ">=") * 1L              # 0/1 matrix
bin <- Matrix::Matrix(bin, sparse = TRUE)
# store as a second layer
epi[["SCENIC"]] <- SetAssayData(epi[["SCENIC"]], slot = "counts", new.data = m)
epi[["SCENIC"]] <- SetAssayData(epi[["SCENIC"]], slot = "data",   new.data = m)
# binarized as a separate assay if you like:
epi[["SCENIC_bin"]] <- CreateAssayObject(counts = bin)
DefaultAssay(epi) <- "SCENIC_bin"
epi <- ScaleData(
  epi,
  assay   = "SCENIC_bin",
  features = rownames(epi[["SCENIC_bin"]]),
  do.center = TRUE,
  do.scale  = TRUE,
  scale.max = 3,      # optional: clip extreme z-scores for prettier heatmaps
  verbose   = FALSE
)

prev <- Matrix::rowMeans(bin)             # fraction ON per regulon
keep <- names(sort(prev, decreasing = TRUE))
DefaultAssay(epi) <- "SCENIC_bin"
mat <- GetAssayData(epi, assay = "SCENIC_bin", slot = "scale.data")[keep, , drop = FALSE]
ord <- rownames(mat)[hclust(dist(mat))$order]  # hierarchical order of regulons

p <- DoHeatmap(
  epi, features = keep, assay = "SCENIC_bin",
  slot = "scale.data", disp.min = -2, disp.max = 2, raster = TRUE
) + scale_fill_gradient2(low="navy", mid="white", high="firebrick3",
                         midpoint=0, limits=c(-2,2), na.value="grey90")
print(p) 

DefaultAssay(epi) <- "SCENIC_bin"
mat <- GetAssayData(epi, assay = "SCENIC_bin", slot = "scale.data")[keep, , drop = FALSE]

p<- pheatmap::pheatmap(
  mat,
  cluster_rows = TRUE,    # y-axis dendrogram ✅
  cluster_cols = FALSE,   # set TRUE if you also want columns clustered (can be heavy)
  show_colnames = FALSE,
  color = colorRampPalette(c("navy","white","firebrick3"))(255),
  breaks = seq(-2, 2, length.out = 256)
)

avg <- AverageExpression(epi, assays = "SCENIC_bin", slot = "scale.data")$SCENIC_bin
pheatmap::pheatmap(avg[keep, ], cluster_rows = F, cluster_cols = TRUE)

options(repr.plot.width = 8, repr.plot.height = 20)
print(p)
# dotplot (size = % ON, color = mean of 0/1)
# DotPlot(epi, features = keep, assay = "SCENIC_bin") + RotatedAxis()

# UMAP for a few regulons (binary)
FeaturePlot(epi, features = c("TP63(+)", "SOX2(+)"),
            assay = "SCENIC_bin", cols = c("grey85","red"))

heatmap(epi[["SCENIC_bin"]]$counts)
### filter by regulon importance

regulon_fname <- file.path(analysis_dir, "pyscenic_res/human_epi_regulons.csv")
regulon_df <- read.csv(regulon_fname, check.names = F)
regulon_df[1]



df2 <- regulon_df %>%
  mutate(weight_pct = parse_number(weight)) %>%  # "53.2%" -> 53.2
  filter(weight_pct > 50) %>%                    # > 50%
  arrange(desc(weight_pct), .by_group = FALSE)

head(df2)
