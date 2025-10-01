## --- Libs ---

library(dplyr)
library(Seurat)
library(patchwork)
library(biomaRt)
library(CellChat)
library(parallelly)
library(future)

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

saveRDS(cellchat, 
        file = file.path(cache_dir, "VU40T_cellchat_res.rds")
)

cellchat <- readRDS(file.path(cache_dir, "VU40T_cellchat_res.rds"))
## visualisation
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

levels(cellchat@idents)

# reorder to split the classes properly
epi_levels <- c("0-E", "1-E", "2-E", "3-E", "4-E", "5-E", "6-E", "7-E")   # <-- edit to your names
fib_levels <- c("0-Fib", "1-Fib", "2-Fib", "3-Fib", "4-Fib", "5-Fib", "6-Fib",
                "7-Fib", "8-Fib", "9-Fib", "10-Fib", "11-Fib", "12-Fib")   # <-- edit to your names

new_order <- c(epi_levels, fib_levels)
cellchat@idents <- factor(cellchat@idents, levels = new_order)



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
# pathways.show <- c("SPP1") 
# 
# netAnalysis_contribution(cellchat, signaling = pathways.show)
# dev.off()
# 
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
fib_idents = grep("Fib", levels(cellchat@idents))
# vertex.receiver = seq(1,length(target_idents)) # Fibroblast clusters will be set to left most Target nodes, Epithelial clusters will be used on the right hand side
epi_idents = grep("E", levels(cellchat@idents))
object <- cellchat ## to avoid object errors
setwd(plot_dir)  ## netVisual exports all plots to the working directory
for (i in 1:length(pathways.show.all)) {
  
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(object = object,
            signaling = pathways.show.all[i],
            vertex.receiver = epi_idents,
            layout = "hierarchy" ,
            out.format = "pdf")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 8, height = 5, units = 'in', dpi = 300)
}
dev.off() ## force X11 device to close after this if it does not

## identify signalling roles

# Compute the network centrality scores - takes 2 secs
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
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

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", 
                                         height = 12,
                                         width = 13)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", 
                                         height = 12,
                                         width = 12)
png("incoming_and_outgoing_signallingRole_heatmap.png", 
    width = 14, height = 8, units = "in", res = 300)
ht1 + ht2
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