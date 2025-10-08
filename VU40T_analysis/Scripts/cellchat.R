## --- Libs ---

library(dplyr)
library(Seurat)
library(patchwork)
library(biomaRt)
library(CellChat)
library(parallelly)
library(future)
library(grid)

#---functs---

### custom heatmap function to reorder the cols to get epi first before fibroblasts
reordered_netAnalysis_signalingRole_heatmap <- function(object, signaling = NULL, pattern = c("outgoing", "incoming","all"), slot.name = "netP",
                                              color.use = NULL, color.heatmap = "BuGn",
                                              title = NULL, width = 10, height = 8, font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE){
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[,i] <- centr[[i]]$outdeg
    incoming[,i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  } else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  } else if (pattern == "all") {
    mat <- t(outgoing+ incoming)
    legend.name <- "Overall"
  }
  
  # reorder to split the classes properly
  epi_levels <- c("0-E", "1-E", "2-E", "3-E", "4-E", "5-E", "6-E", "7-E")
  fib_levels <- c("0-Fib", "1-Fib", "2-Fib", "3-Fib", "4-Fib", "5-Fib", "6-Fib",
                  "7-Fib", "8-Fib", "9-Fib", "10-Fib", "11-Fib", "12-Fib")   
  
  new_order <- c(epi_levels, fib_levels)
  mat <- mat[, new_order, drop = FALSE]
  
  
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  } else {
    title <- paste0(paste0(legend.name, " signaling patterns"), " - ",title)
  }
  
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
  mat[mat == 0] <- NA
  
  
  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)
  
  df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  names(color.use) <- colnames(mat)
  col_annotation <- ComplexHeatmap::HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha2 = ComplexHeatmap::HeatmapAnnotation(Strength = ComplexHeatmap::anno_barplot(colSums(mat.ori), 
                                                                                  border = FALSE,
                                                                                  gp = gpar(fill = color.use, col=color.use)), 
                                                                                  show_annotation_name = FALSE)
  
  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum)*1.1, max(pSum)*1.5, length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }
  
  ha1 = ComplexHeatmap::rowAnnotation(Strength = ComplexHeatmap::anno_barplot(pSum, border = FALSE), show_annotation_name = FALSE)
  
  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  } else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), round(max(mat, na.rm = T), digits = 1))
  }
  ht1 = ComplexHeatmap::Heatmap(mat, col = color.heatmap.use, na_col = "white", name = "Relative strength",
                  bottom_annotation = col_annotation, top_annotation = ha2, right_annotation = ha1,
                  cluster_rows = cluster.rows,cluster_columns = cluster.rows,
                  row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                  width = unit(width, "cm"), height = unit(height, "cm"),
                  column_title = title,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                              border = NA, at = legend.break,
                                              legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
  )
  #  draw(ht1)
  return(ht1)
}


# --- determine number of cores ---
available <- parallelly::availableCores()
cores_to_use <- max(1, available - 1)
message("🧠 Using ", cores_to_use, " core(s) for cellchat.")

## to avoid hitting that 20gb quota on my home dir,
## --- directories ---


proj_dir <- "/rds/projects/g/gendood-3dmucosa/"
analysis_dir <- file.path(proj_dir, "scRNAseqAnalysis")
git_dir <- file.path(analysis_dir, "OralMucosa/VU40T_analysis")
cache_dir <- file.path(proj_dir, "rds_cache")
tbls_dir <- file.path(git_dir, "Integrated/CellChat_results")
plot_dir <- file.path(tbls_dir, "Plots")

## dir check

dir_chk_list <- c(git_dir, cache_dir,tbls_dir, plot_dir )

for (path in dir_chk_list){
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}


epis <- readRDS(file.path(cache_dir, "VU40T_combined_joined_0.3_res_humanOnly.rds"))

fibros <-  readRDS(file.path(cache_dir, "VU40T_combined_joined_0.6_res_mouseOnly.rds"))

## omit epithelial cluster in fibros to avoid overcounting cells

# 2) Keep all idents except "13"
fibros@graphs <- list()
fibros <- subset(fibros, seurat_clusters != "13")

### quick violin plots for fig 6 of TGFA and TGFBR3 before merging obj 
### showing 1-way communication of the L-R that is dependent on epi contact
### w fibroblasts in the model

genestoplot <- c("TGFA", "TGFBR3")

p <- lapply(seq_along(genestoplot), function(i) {
  g <- genestoplot[i]
  gp <- VlnPlot(
    epis,
    features = g,
    group.by = "seurat_clusters",
    pt.size  = 0,
    combine  = TRUE
  ) + 
    NoLegend() +
    labs(title = g, x = "",
         y = "") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14, , vjust = 0.4),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 8, angle = 0),
      axis.title.x = element_blank(),
      axis.text.x  = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14),
      axis.ticks.x = element_blank(),
      plot.margin = margin(2, 6, 2, 5)
    )
})


p2 <- lapply(seq_along(genestoplot), function(i) {
  g <- genestoplot[i]
  VlnPlot(
    fibros,
    features = g,
    group.by = "seurat_clusters",
    pt.size  = 0,
    combine  = TRUE) + 
    NoLegend() +
    labs(title = stringr::str_to_sentence(tolower(g)), 
         x = "", y = "") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14, vjust = 0.3),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 8, angle = 0),
      axis.title.x = element_blank(),
      axis.text.x  = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14),
      axis.ticks.x = element_blank(),
      plot.margin = margin(2, 6, 2, 5)
    )
})
plots_list <- append(p, p2)
p_stack <- wrap_plots(
  plots_list,
  ncol = 2,
  guides = "keep"
) +
  plot_layout(
    guides = "collect"
  ) +
  plot_annotation(
    title = "TGFA-TGFBR3 L-R Epithelial and Fibroblast Expression"
  ) &
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.margin = margin(20, 5, 20, 5),
    panel.spacing = unit(2, "lines")
  )
print(p_stack)

png(file.path(plot_dir, "VlnPlot_TGFA-TGFBR3_L-R.png"),
    width = 10, height = 10, units = "in", res = 300)
print(p_stack)
dev.off()


epis$seurat_clusters <- sapply(epis$seurat_clusters, function(x) paste0(x, "-human"))
fibros$seurat_clusters <- sapply(fibros$seurat_clusters, function(x) paste0(x, "-mouse"))

## convert fibro gene names to human homologues

# Connect to Ensembl human mart (homolog info is already stored here)
mmart <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl",
  host = "https://www.ensembl.org"
)

mouse_genes <- rownames(fibros)

# ✅ These attributes are from the *same attribute page* (homologs)
conv <- getBM(
  attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"),
  filters    = "external_gene_name",
  values     = unique(mouse_genes),
  mart       = mmart
)

## mapping dict
mouse_to_human <- setNames(
  toupper(conv$hsapiens_homolog_associated_gene_name),
  conv$external_gene_name
)
mouse_to_human <- mouse_to_human[mouse_to_human != "" & !is.na(mouse_to_human)]

# 5) Create the new symbol vector for your rows (fallback to original mouse symbol if no homolog)
new_symbols <- toupper(ifelse(mouse_genes %in% names(mouse_to_human),
                              mouse_to_human[mouse_genes],
                              mouse_genes))

## remove any poss duplicates by aggregating counts w new symbols
cnt <- LayerData(fibros, assay = DefaultAssay(fibros), layer = "counts")

coo <- summary(cnt)                 # columns: i, j, x  (1-based)

# 2) Map each nonzero row to its new gene symbol, drop NAs/empty
genes_for_nz <- new_symbols[coo$i]
keep <- !is.na(genes_for_nz) & genes_for_nz != ""
coo <- coo[keep, , drop = FALSE]
genes_for_nz <- genes_for_nz[keep]

# 3) Choose gene order (levels). Use first-appearance order; change if you want a custom order.
gene_levels <- unique(genes_for_nz)

# 4) Build aggregated sparse matrix (duplicates auto-summed)
agg_cnt <- Matrix::sparseMatrix(
  i = match(genes_for_nz, gene_levels),
  j = coo$j,
  x = coo$x,
  dims = c(length(gene_levels), ncol(cnt)),
  dimnames = list(gene_levels, colnames(cnt))
)


# If you also have data/scale.data and want to fully rebuild, recompute after creating the assay.
new_assay <- CreateAssayObject(counts = agg_cnt)
## sanity check
if (all(colnames(new_assay) == colnames(fibros))){
  message(
    "Mouse fibroblast genes mapped to human homologues successfully.\n
Ignore next warning."
  )
} else {
  warning("cell ID's have been altered in gene count aggregation after mouse -> human homologue mapping")
}

no_agg_genes <- length(rownames(fibros)) - length(rownames(new_assay))
message(paste0("number of genes aggregated after homologue renaming = ", no_agg_genes))
## ignore any warning here if sanity check successful
fibros[[DefaultAssay(fibros)]] <- new_assay

### merge objs back with each other

DefaultAssay(epis)   <- "RNA"
DefaultAssay(fibros) <- "RNA"

EpiFib <- merge(
  epis, y = fibros,
  add.cell.ids = c("EPI","FIB"),
  project = "EpiFib"
)

DefaultAssay(EpiFib) <- "RNA"
EpiFib <- JoinLayers(EpiFib, assay = "RNA")
EpiFib <- NormalizeData(EpiFib)
# EpiFib$celltype <- ifelse(grepl("^EPI_", colnames(EpiFib)), "Epithelia", "Fibroblasts")

data.input <- LayerData(EpiFib, assay = "RNA", layer = "data")
stopifnot(nrow(data.input) > 0, ncol(data.input) > 0)
group_col <- "seurat_clusters"



names(EpiFib@meta.data)[names(EpiFib@meta.data) == "sample"] <- "samples"
EpiFib@meta.data$samples <- factor(EpiFib@meta.data$samples)
meta <- EpiFib@meta.data[colnames(EpiFib), , drop = FALSE]  # align rows to cell

# make grouping a factor and drop NA rows
meta[[group_col]] <- droplevels(factor(meta[[group_col]]))
meta <- meta[!is.na(meta[[group_col]]), , drop = FALSE]
stopifnot(all(rownames(meta) == colnames(data.input)))
# align columns
data.input <- data.input[, rownames(meta)]

cellchat  <- createCellChat(object = data.input, 
                            meta = meta, 
                            group.by = group_col)

cellchat@DB <- CellChatDB.human

cellchat <- subsetData(cellchat)
cat("data.signaling dims:", dim(cellchat@data.signaling), "\n")


cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 0) We'll use the raw signaling matrix (you already saw project is empty)
data_mat <- cellchat@data.signaling
genes_in_data <- rownames(data_mat)

# helper: does 'g' refer to a complex in the DB?
is_complex <- function(g) g %in% rownames(cellchat@DB$complex)

# helper: are ALL subunits for 'g' present in the data?
complex_ok <- function(g) {
  if (!is_complex(g)) return(g %in% genes_in_data)
  subs <- unlist(cellchat@DB$complex[g, grep("^subunit", names(cellchat@DB$complex)), drop = FALSE],
                 use.names = FALSE)
  subs <- subs[subs != ""]
  length(subs) > 0 && all(subs %in% genes_in_data)
}

# 1) filter LR table to pairs whose ligand+receptor (and their subunits) all exist
lr0 <- cellchat@LR$LRsig
ok_lig <- vapply(lr0$ligand,   complex_ok, logical(1))
ok_rec <- vapply(lr0$receptor, complex_ok, logical(1))
lr_keep <- lr0[ok_lig & ok_rec, , drop = FALSE]

if (nrow(lr_keep) == 0L) {
  stop("After checking complex subunits, no LR pairs remain. 
       Likely many complex subunits are missing from your gene set.")
}

# (optional) see what was dropped
dropped <- lr0[!(ok_lig & ok_rec), c("ligand","receptor"), drop = FALSE]
head(dropped)

# 2) replace LR list and re-run the standard steps
cellchat@LR$LRsig <- lr_keep

# ensure idents are named & aligned to matrix columns
if (is.null(names(cellchat@idents))) names(cellchat@idents) <- colnames(data_mat)
cellchat@data.signaling <- data_mat[, names(cellchat@idents), drop = FALSE]
cellchat@idents <- droplevels(cellchat@idents)

# avoid sparse->data.frame issues in aggregate()
if (!inherits(cellchat@data.signaling, "matrix")) {
  cellchat@data.signaling <- as.matrix(cellchat@data.signaling)
}
future::plan(strategy = "multisession", workers = cores_to_use)
# 3) compute
## takes 10 mins on rstudio session with 8 cores used
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, type = "triMean")
# # then:
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# make sure idents are named by cell barcodes
if (is.null(names(cellchat@idents))) {
  names(cellchat@idents) <- colnames(cellchat@data.signaling)
}

# rename levels: "...-human" -> "...-epithelial", "...-mouse" -> "...-fibroblast"
lv <- levels(cellchat@idents)
lv_new <- sub("-human$", "-E", lv)
lv_new <- sub("-mouse$", "-Fib", lv_new)
levels(cellchat@idents) <- lv_new

## and rownames in counts and weights for visualisation

fix_names <- function(x) {
  stringr::str_replace_all(x, c("-human$" = "-E", "-mouse$" = "-Fib"))
}

for (slot in c("count", "weight", "prob", "pval")) {
  m <- cellchat@net[[slot]]
  rownames(m) <- fix_names(rownames(m))
  colnames(m) <- fix_names(colnames(m))
  cellchat@net[[slot]] <- m
}

rownames(cellchat@netP$prob) <- fix_names(rownames(cellchat@netP$prob))
colnames(cellchat@netP$prob) <- fix_names(colnames(cellchat@netP$prob))


# Compute the network centrality scores - takes 2 secs
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

saveRDS(cellchat, 
        file = file.path(cache_dir, "VU40T_cellchat_res.rds")
)

cellchat <- readRDS(file.path(cache_dir, "VU40T_cellchat_res.rds"))
## visualisation


levels(cellchat@idents)

# reorder to split the classes properly
epi_levels <- c("0-E", "1-E", "2-E", "3-E", "4-E", "5-E", "6-E", "7-E")   # <-- edit to your names
fib_levels <- c("0-Fib", "1-Fib", "2-Fib", "3-Fib", "4-Fib", "5-Fib", "6-Fib",
                "7-Fib", "8-Fib", "9-Fib", "10-Fib", "11-Fib", "12-Fib")   # <-- edit to your names

new_order <- c(epi_levels, fib_levels)





# quick check
table(cellchat@idents)
# vertex weights aligned to matrix order
library(magick)

# vertex weights aligned to matrix order
gs <- table(cellchat@idents)
gs <- as.numeric(gs[match(rownames(cellchat@net$count), names(gs))])

f1 <- file.path(plot_dir, "VU40T_num_interactions.png")
f2 <- file.path(plot_dir, "VU40T_weight_interactions.png")
fout <- file.path(plot_dir, "VU40T_num_and_weight_interactions.png")

# 1) draw number-of-interactions
png(f1, width = 2000, height = 2000, res = 300)
par(mar = c(1,1,2,1), xpd = TRUE)
netVisual_circle(cellchat@net$count,
                 vertex.weight = gs, weight.scale = TRUE,
                 label.edge = FALSE, title.name = "Number of interactions")
dev.off()

# 2) draw interaction strength
png(f2, width = 2000, height = 2000, res = 300)
par(mar = c(1,1,2,1), xpd = TRUE)
netVisual_circle(cellchat@net$weight,
                 vertex.weight = gs, weight.scale = TRUE,
                 label.edge = FALSE, title.name = "Interaction strength")
dev.off()

# 3) stitch them side-by-side
img1 <- image_read(f1)
img2 <- image_read(f2)
combo <- image_append(c(img1, img2))  # horizontally
image_write(combo, path = fout)

message("Saved: ", fout)


# inputs
mat <- cellchat@net$weight  # or cellchat@net$count
# vertex weights aligned to matrix order
groupSize <- table(cellchat@idents)
groupSize <- as.numeric(groupSize[match(rownames(mat), names(groupSize))])

# optional: order panels epi first, then fibro; tweak patterns to your labels
grp <- rownames(mat)
num <- suppressWarnings(as.integer(sub("^([0-9]+).*", "\\1", grp)))
is_epi <- grepl("epi|E|epithelial$", grp, ignore.case = TRUE)
ord <- order(!is_epi, num, grp)     # epi (TRUE) first, then fibro, numeric within each
mat <- mat[ord, , drop = FALSE]


# ---- per-panel rendering ----
tmp_dir <- file.path(tempdir(), "cellchat_circles_21")
dir.create(tmp_dir, showWarnings = FALSE)
png_files <- character(nrow(mat))
ewmax <- max(mat, na.rm = TRUE)      # consistent edge scaling



for (i in seq_len(nrow(mat))) {
  f_i <- file.path(tmp_dir, sprintf("circle_%02d_%s.png", i, make.names(rownames(mat)[i])))
  png_files[i] <- f_i
  png(f_i, width = 1600, height = 1600, res = 300)
  par(mar = c(1,1,2,1), xpd = TRUE)
  # isolate sender i
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(
    mat2,
    vertex.weight   = groupSize,
    weight.scale    = TRUE,
    edge.weight.max = ewmax,
    label.edge      = FALSE,
    title.name      = rownames(mat)[i]
  )
  dev.off()
}

# ---- stitch into a grid ----
fout <- file.path(plot_dir, "VU40T_weight_interactions_from_each_grp.png")
imgs <- image_read(png_files)

# pick a neat grid for 21 panels (e.g., 7 cols × 3 rows)
ncol_grid <- 7
n <- length(imgs)
if (n %% ncol_grid != 0) {
  pad <- ncol_grid - (n %% ncol_grid)
  blank <- image_blank(width = image_info(imgs[1])$width,
                       height = image_info(imgs[1])$height,
                       color = "white")
  imgs <- c(imgs, rep(list(blank), pad))
}

rows <- split(imgs, rep(seq_len(length(imgs) / ncol_grid), each = ncol_grid))
row_imgs <- lapply(rows, function(x) image_append(image_join(x)))     # horizontal rows
final <- Reduce(function(a, b) image_append(image_join(c(a, b)), stack = TRUE), row_imgs)

image_write(final, path = fout)
message("Saved: ", fout)

# ### do pathway plots now
pathways.show <- c("SPP1") 
# 
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")








# pathways.show <- c("LAMININ") 
# 
# netAnalysis_contribution(cellchat, signaling = pathways.show)
# dev.off()
# 
# pathways.show <- c("EGF") 
# 
# netAnalysis_contribution(cellchat, signaling = pathways.show)
# dev.off()
# 
# pathways.show <- c("PDGF") 
# 
# netAnalysis_contribution(cellchat, signaling = pathways.show)
# dev.off()
# Hierarchy plot

# Heatmap
# par(mfrow=c(1,1))
# netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
# dev.off()



# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
# fib_idents = grep("Fib", levels(cellchat@idents))
# # vertex.receiver = seq(1,length(target_idents)) # Fibroblast clusters will be set to left most Target nodes, Epithelial clusters will be used on the right hand side
# epi_idents = grep("E", levels(cellchat@idents))

# Combine and sort lexicographically
lexicographic_order <- sort(c(levels(cellchat@idents)))
## use indices for lexicographic epi_idents before reordering which didn't work
epi_idents_old <- grep("E", lexicographic_order)

object <- cellchat ## to avoid object errors
setwd(plot_dir)  ## netVisual exports all plots to the working directory
for (i in 1:length(pathways.show.all)) {
  
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(object = object,
            signaling = pathways.show.all[i],
            vertex.receiver = epi_idents_old,
            layout = "hierarchy" ,
            out.format = "pdf")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 8, height = 5, units = 'in', dpi = 300)
}
dev.off() ## force X11 device to close after this if it does not




## identify signalling roles


# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups

# while (dev.cur() != 1 && names(dev.cur()) != "RStudioGD") dev.off()
for (i in 1:length(pathways.show.all)) {
  pdf(paste0(pathways.show.all[i], "_signallingRole_network.pdf"), width = 8, height = 2.5)
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show.all[i],
                                    width = 8, height = 2.5, font.size = 10,
                                    cluster.rows = F,
                                    cluster.cols = F)
  dev.off()
}

signalling_groups <- list(
  Group1 = c("LAMININ", "COLLAGEN", "SPP1", "FN1"),
  Group2 = c("ncWNT", "EGF", "NRG", "MIF"),
  Group3 = c("FGF", "PDGF", "VEGF", "ANGPTL", "CypA"),
  Group4 = c("NOTCH", "ADGRL", "THBS", "SLITRK", "MK", "TENASCIN"),
  Group5 = c("TRAIL", "EPHA", "EPHB", "APP", "SLIT", "VISFATIN")
)

## extract data from cellchat object used for 
patterns <- c("incoming", "outgoing")

object <- cellchat
slot.name <- "netP"
centr <- slot(object, slot.name)$centr
outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
dimnames(outgoing) <- list(levels(object@idents), names(centr))
dimnames(incoming) <- dimnames(outgoing)
for (i in 1:length(centr)) {
  outgoing[,i] <- centr[[i]]$outdeg
  incoming[,i] <- centr[[i]]$indeg
}

for (pattern in patterns){
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  } else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  } else if (pattern == "all") {
    mat <- t(outgoing+ incoming)
    legend.name <- "Overall"
  }
  mat <- mat[, new_order, drop = FALSE]
  
  # if (is.null(title)) {
  #   title <- paste0(legend.name, " signaling patterns")
  # } else {
  #   title <- paste0(paste0(legend.name, " signaling patterns"), " - ",title)
  # }
  
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
  mat[mat == 0] <- NA
  print(mat)
  # write.csv(mat, file.path(tbls_dir, paste0(pattern, "_communication_centrality_scores.csv")))

}

#   
#   
#   if (is.null(color.use)) {
#     color.use <- scPalette(length(colnames(mat)))
#   }
#   color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)
#   
#   df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
#   names(color.use) <- colnames(mat)
#   col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
#                                       show_legend = FALSE, show_annotation_name = FALSE,
#                                       simple_anno_size = grid::unit(0.2, "cm"))
#   ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
#   
#   pSum <- rowSums(mat.ori)
#   pSum.original <- pSum
#   pSum <- -1/log(pSum)
#   pSum[is.na(pSum)] <- 0
#   idx1 <- which(is.infinite(pSum) | pSum < 0)
#   if (length(idx1) > 0) {
#     values.assign <- seq(max(pSum)*1.1, max(pSum)*1.5, length.out = length(idx1))
#     position <- sort(pSum.original[idx1], index.return = TRUE)$ix
#     pSum[idx1] <- values.assign[match(1:length(idx1), position)]
#   }
#   
#   ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), show_annotation_name = FALSE)
#   
#   if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
#     legend.break <- max(mat, na.rm = T)
#   } else {
#     legend.break <- c(round(min(mat, na.rm = T), digits = 1), round(max(mat, na.rm = T), digits = 1))
#   }
#   ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = "Relative strength",
#                 bottom_annotation = col_annotation, top_annotation = ha2, right_annotation = ha1,
#                 cluster_rows = cluster.rows,cluster_columns = cluster.rows,
#                 row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
#                 width = unit(width, "cm"), height = unit(height, "cm"),
#                 column_title = title,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
#                 heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
#                                             border = NA, at = legend.break,
#                                             legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
#   )
#   #  draw(ht1)
#   return(ht1)
# }
# }


for (group in seq_along(signalling_groups)) {
  group_name <- names(signalling_groups)[group]
  group_signals <- signalling_groups[[group]]
  
  ht1 <- reordered_netAnalysis_signalingRole_heatmap(cellchat,
                                           pattern = "outgoing",
                                           height = 8,
                                           width = 13,
                                           signaling = group_signals,cluster.cols = T)
  
  ht2 <- reordered_netAnalysis_signalingRole_heatmap(cellchat,
                                           pattern = "incoming",
                                           height = 8,
                                           width = 12,
                                           signaling = group_signals, cluster.cols = T)

  png(file.path(plot_dir, paste0("incoming_and_outgoing_signalingRole_heatmap_", group_name, ".png")),
      width = 14, height = 5, units = "in", res = 300)
  
  print(ht1 + ht2) 
  dev.off()
}

ht1 <- reordered_netAnalysis_signalingRole_heatmap(cellchat,
                                                   pattern = "outgoing",
                                                   height = 13,
                                                   width = 12)

ht2 <- reordered_netAnalysis_signalingRole_heatmap(cellchat,
                                                   pattern = "incoming",
                                                   height = 13,
                                                   width = 12)

png(file.path(plot_dir, "incoming_and_outgoing_signalingRole_heatmap_allPathways.png"),
    width = 13, height = 7, units = "in", res = 300)

print(ht1 + ht2) 
dev.off()

## pattern analysis

library(NMF)
library(ggalluvial)

#outgoing pattern analysis

png("outgoingSelectK_cophenetic_silhouette_scores.png", 
    width = 14, height = 8, units = "in", res = 300)
selectK(cellchat, pattern = "outgoing")
dev.off()
nPatterns = 3 ## cophenetic drops dramatically after this, silhouette drops at 2 though not as relatively sharply
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# river plot
png("outgoing_patterns_river.png", 
    width = 12, height = 6, units = "in", res = 300)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()
# dot plot looks confusing af
# netAnalysis_dot(cellchat, pattern = "outgoing")

#incoming pattern analysis

png("incomingSelectK_cophenetic_silhouette_scores.png", 
    width = 14, height = 8, units = "in", res = 300)
selectK(cellchat, pattern = "incoming")
dev.off()

nPatterns = 2 # both begin to drop quite a lot after that
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
png("incoming_patterns_river.png", 
    width = 12, height = 6, units = "in", res = 300)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()



set1 <- c("LAMININ", "COLLAGEN", "SPP1", "FN1", "ncWNT", "EGF", "NRG", "MIF")
set2 <- c("NOTCH", "ADGRL", "SLITRK")


### for pathways we care about from NMF analysis
pdf("Set1_CirclePlots_CellChat.pdf", width = 6, height = 6)

for (pw in set1) {
  if (!(pw %in% cellchat@netP$pathways)) {
    warning(paste("❌ Pathway", pw, "not found."))
    next
  }
  
  CellChat::netVisual_aggregate(
    object = cellchat,
    signaling = pw,
    layout = "circle",
    remove.isolate = TRUE
  )
  title(main = pw, cex.main = 1.5, font.main = 2)
}

dev.off()

pdf("Set2_CirclePlots_CellChat.pdf", width = 6, height = 6)

for (pw in set2) {
  if (!(pw %in% cellchat@netP$pathways)) {
    warning(paste("❌ Pathway", pw, "not found."))
    next
  }
  
  CellChat::netVisual_aggregate(
    object = cellchat,
    signaling = pw,
    layout = "circle",
    remove.isolate = TRUE
  )
  title(main = pw, cex.main = 1.5, font.main = 2)
}

dev.off()


# Load required packages only if not already loaded
if (!require("UpSetR")) install.packages("UpSetR")


library(UpSetR)
library(dplyr)

# ligand and receptor pairs
receptors <- c(
  "ITGA2_ITGB1", "ITGA3_ITGB1", "ITGA6_ITGB1", "ITGAV_ITGB8",
  "ITGA6_ITGB4", "ITGAV_ITGB1", "ITGAV_ITGB6", 
  "CD44", "SDC1", "SDC4"
)

# Pathways of interest
pathways_of_interest <- c("LAMININ", "COLLAGEN", "FN1", "SPP1")

# 2. Subset communication dataframe from CellChat
comm_df <- subsetCommunication(cellchat) 

# Filter the communication data frame to those pathways
comm_subset <- comm_df %>%
  filter(pathway_name %in% pathways_of_interest)

# Extract receptor list per pathway
receptors_by_pathway <- comm_subset %>%
  group_by(pathway_name) %>%
  summarise(receptors = list(unique(receptor))) %>%
  tibble::deframe()  # converts to named list

df <- tibble::tibble(
  pathway = rep(names(receptors_by_pathway), times = sapply(receptors_by_pathway, length)),
  receptor = unlist(receptors_by_pathway)
)

binary_mat <- df %>%
  mutate(value = 1) %>%
  tidyr::pivot_wider(names_from = pathway, values_from = value, values_fill = 0)

png(file.path(plot_dir, "UpsetPlot_receptors_in_SPP1_laminin_col_fn1.png"),
    width = 8, height  = 6, units = "in", res = 300)
p<-ComplexUpset::upset(binary_mat, intersect = names(receptors_by_pathway), name = "Receptors by Pathway")
print(p)
dev.off()


# # Filter and prepare the data
# dot_data <- comm_df %>%
#   filter(pathway_name %in% pathways_of_interest) %>%
#   filter(receptor %in% receptors) %>%
#   mutate(
#     safe_pval = pmax(pval, 1e-300)  # Prevent division by zero
#   ) %>%
#   group_by(pathway_name, interaction_name_2) %>%
#   summarise(
#     mean_prob = mean(prob, na.rm = TRUE),
#     mean_pval = mean(safe_pval, na.rm = TRUE),
#     .groups = "drop"
#   )

# # Ensure factors are ordered
# dot_data$pathway_name <- factor(dot_data$pathway_name, levels = pathways_of_interest)
# dot_data$interaction_name_2 <- factor(dot_data$interaction_name_2, levels = rev(unique(dot_data$interaction_name_2)))
# 
# # Plot

# png(file.path(plot_dir, "Dotplot_L-R_pairs_SPP1_laminin_col_fn1.png"),
#     width = 4, height  = 16, units = "in", res = 300)
# p <- ggplot(dot_data, aes(x = pathway_name, y = interaction_name_2)) +
#   geom_point(aes(size = 1 / mean_pval, color = mean_prob)) +
#   viridis::scale_color_viridis(option = "viridis", name = "Comm. Prob") +
#   scale_size_continuous(
#     name = "p-value",
#     breaks = 1 / c(1e-2, 1e-5, 1e-10, 1e-50, 1e-100),
#     labels = c("1e-2", "1e-5", "1e-10", "1e-50", "1e-100"),
#     range = c(2, 6)
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.title = element_blank(),
#     legend.position = "right",
#     panel.grid = element_blank()
#   ) +
#   ggtitle("L-R Pairs: Collagen,\nFN1, Laminin, SPP1")
# print(p)
# dev.off()

receptor_counts <- comm_subset %>%
  group_by(pathway_name, receptor) %>%
  summarise(lr_count = n_distinct(interaction_name_2), .groups = "drop")

all_combinations <- expand.grid(
  receptor = unique(comm_subset$receptor),
  pathway_name = pathways_of_interest
)

receptor_counts_complete <- all_combinations %>%
  left_join(receptor_counts, by = c("receptor", "pathway_name")) %>%
  mutate(lr_count = ifelse(is.na(lr_count), 0, lr_count))

receptor_counts_complete$receptor <- stringr::str_replace_all(receptor_counts_complete$receptor ,"_","/")
receptor_counts_complete$pathway_name <- factor(receptor_counts_complete$pathway_name, levels = c("LAMININ", "COLLAGEN", "FN1", "SPP1"))

png(file.path(plot_dir, "Heatmap_freq_Common_Sig_receptors_SPP1_laminin_col_fn1.png"),
    width = 9, height  = 8, units = "in", res = 600)

p<- ggplot(receptor_counts_complete, aes(x = pathway_name, y = receptor, fill = lr_count)) +
    geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal(base_size = 12) +
    labs(title = "Frequency of L-R Pairs per Receptor common to each Pathway",
         x = "Pathway", y = "Receptor", fill = "L-R pair count") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14),
      panel.grid = element_blank()
    )
print(p)
dev.off()

# Create binary matrix
binary_matrix <- sapply(pathways_of_interest, function(pw) {
  lr_in_pathway <- comm_subset %>%
    filter(pathway_name == pw) %>%
    pull(receptor)
  
  as.integer(receptors %in% lr_in_pathway)
})

## set order
colnames(binary_matrix) <- factor(colnames(binary_matrix), levels = c("LAMININ", "COLLAGEN", "FN1", "SPP1"))


# Set rownames for clarity
rownames(binary_matrix) <- receptors
rownames(binary_matrix) <- stringr::str_replace_all(rownames(binary_matrix) ,"_","/")
# View
print(binary_matrix)


df_long <- reshape2::melt(binary_matrix)
colnames(df_long) <- c("Receptor", "Pathway", "Involved")

png(file.path(plot_dir, "Heatmap_Common_Sig_receptors_SPP1_laminin_col_fn1.png"),
    width = 5, height  = 7, units = "in", res = 300)
# Plot heatmap
p<-ggplot(df_long, aes(x = Pathway, y = Receptor, fill = factor(Involved))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "steelblue")) +
  theme_minimal(base_size = 12) +
  labs(title = "Common Receptor Involvement\nin Signalling Pathways",
       fill = "Involved"
       ) + 
  NoLegend()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(hjust = 0.5),
    plot.title = element_text(hjust = 0.5)
  )
print(p)
dev.off()


## look at functional similarity
# 
# cellchat <- computeNetSimilarity(cellchat, type = "functional")
# cellchat <- netEmbedding(cellchat, type = "functional")
# cellchat <- netClustering(cellchat, type = "functional")
# # Visualization in 2D-space
# netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
# # netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

# cellchat <- computeNetSimilarity(cellchat, type = "structural")
# reticulate::use_virtualenv()
# cellchat <- netEmbedding(cellchat, type = "structural")
# cellchat <- netClustering(cellchat, type = "structural")
# # Visualization in 2D-space
# netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
# netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)