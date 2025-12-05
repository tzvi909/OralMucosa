### libs

library(dplyr)
library(Seurat)
library(scater)
library(scuttle)
library(ggplot2)
library(tibble)
library(SingleCellExperiment)
library(patchwork)
library(optparse)

### opts

option_list <- list(
  make_option(c("-s", "--species"),
              type = "character",
              default = "human",
              help = "species to use [default = %default]. Options: human or mouse",
              metavar = "character"),
  make_option(c("-i", "--input"),
              type = "character",
              help = "input RDS dir to load",
              metavar = "character"),
  make_option(c("-n", "--n_cores"),
              type = "integer",
              default = 1,
              help = "Number of cores to use to run  [default = %default]",
              metavar = "character"),
  make_option(c("-c", "--cachedir"),
              type = "character",
              default = "rds_cache",
              help = "Path to save intermediate RDS caches [default = %default]",
              metavar = "character"),
  make_option(c("-o", "--outdir"),
              type = "character",
              default = ".",
              help = "Path to save results [default = %default]",
              metavar = "character"),
  make_option(c("-m", "--mito_path"),
              type = "character",
              help = "Path to load in MT-gene list for hg38",
              metavar = "character")
)


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
print(opt)

n_cores <- opt$n_cores
cache_dir <- opt$cachedir
mtx_dir <- opt$input
outdir <- opt$outdir
mito_path <- opt$mito_path

# Store selected species
species <- tolower(opt$species)

# Validate species input
if (!species %in% c("human", "mouse")) {
  stop("Invalid species. Use 'human' or 'mouse'.")
}


### dir setup

out_prefix_dir <- file.path(outdir, paste0("Seperate_samples/", stringr::str_to_title(species)))

proj_dir <- "/rds/projects/g/gendood-3dmucosa/"
analysis_dir <- file.path(proj_dir, "scRNAseqAnalysis/")
git_dir <- file.path(analysis_dir, "OralMucosa/VU40T_analysis")
out_prefix_dir <- file.path(git_dir, paste0("Seperate_samples/", stringr::str_to_title(species)))
plotQC_dir <- file.path(out_prefix_dir, "Plots/QC")
cache_dir <- file.path(proj_dir, "rds_cache")
mito_path <- file.path(proj_dir, "BaseSpace/LPS_VU40T_QC_and_counts/cellranger/mkref/cellranger_reference/genes/MTgenes.txt")
## check for dirs recursively

chk_dir_list <- list(cache_dir, out_prefix_dir, plotQC_dir)

for (path in chk_dir_list){
  if(!(dir.exists(path))){
    dir.create(path, recursive = T)
  }
}

# if (species == "human"){
#   mtx_dir <- file.path(proj_dir, "BaseSpace/LPS_VU40T_QC_and_counts/cellranger/mtx_conversions")
# } else {
#   mtx_dir <- file.path(proj_dir, "BaseSpace/LPS_VU40T_QC_and_counts_Mouse/cellranger/mtx_conversions")
# }

Samples <- list.files(path = mtx_dir, 
                      pattern = "^P",
                      recursive = F
) 

cat("samples identified for analysis: ", Samples)


##check we only have 4 samples
## get only cellbender seurat_obj paths
seurat_objs <- list.files(path = mtx_dir, 
                          pattern = "^P.*cellbender.*.seurat.rds$",
                          full.names = T,
                          recursive = T
)



SeuratList <- list()
for(i in seurat_objs){
  SeuratList[[i]] <- LoadSeuratRds(i)
  Idents(SeuratList[[i]]) <- SeuratList[[i]]$sample
}

if (length(SeuratList) == 4){
  message("All samples loaded into Seurat from this paper")
} else {
  message("Partial QC and analysis started")
}

### ENS ID conversion to HGNC for mouse genome
if (species == "mouse"){
  library(biomaRt)
  
  # 1. Connect to Ensembl (mouse)
  ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
  for(i in seq_along(SeuratList)){
    seurat_obj <- SeuratList[[i]]
    # 2. Get Ensembl IDs from your Seurat object (assumes default assay is RNA)
    ens_ids <- rownames(seurat_obj[["RNA"]])
    
    # 3. Query BioMart for gene symbol mappings
    gene_map <- getBM(
      attributes = c("ensembl_gene_id", "mgi_symbol"),
      filters = "ensembl_gene_id",
      values = ens_ids,
      mart = ensembl
    )
    
    # 4. Clean
    gene_map <- gene_map[gene_map$mgi_symbol != "", ]
    gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]
    
    # 5. Replace rownames in Seurat object
    common_ids <- intersect(rownames(seurat_obj[["RNA"]]), gene_map$ensembl_gene_id)
    new_names <- gene_map$mgi_symbol[match(common_ids, gene_map$ensembl_gene_id)]
    new_names <- make.unique(new_names)
    
    rownames(seurat_obj[["RNA"]])[match(common_ids, rownames(seurat_obj[["RNA"]]))] <- new_names
    
    SeuratList[[i]] <- seurat_obj
    message("Mouse ENS IDs changed to HGNC sucessfully for sample ", unique(SeuratList[[i]]$sample) )
  }
}

# Replace SeuratList names with simplified format like P18_LPS-N, P22_LPS-P, etc.
names(SeuratList) <- gsub("^.*/(P[0-9]+)_LPS-([NP]).*$", "\\1_LPS-\\2", names(SeuratList))

## Mito genes were isolated from hg38 GTF (chrM)

has_mito_prefix <- any(sapply(SeuratList, function(obj) {
  any(grepl("^MT-", rownames(obj[["RNA"]]))) |
    any(grepl("^mt-", rownames(obj[["RNA"]])))
}))

if (!has_mito_prefix) {
  if (file.exists(mito_path)) {
    mito_genes <- read.delim(mito_path, sep = "\t", header = FALSE)
    mito_genes <- c(t(mito_genes))
  } else {
    stop("❌ Mitochondrial gene list not found at: ", mito_path)
  }
}


## thresholds worked out from doing histogram of log-transformed count data
thresholds_df <- data.frame(
  sample = c("P18_LPS-N", "P18_LPS-P", "P22_LPS-N", "P22_LPS-P"),
  cutoff = c(300, 350, 300, 300)
)

has_mito_prefix <- any(sapply(SeuratList, function(obj) {
  any(grepl("^MT-", rownames(obj[["RNA"]]))) |
    any(grepl("^mt-", rownames(obj[["RNA"]])))
}))

if (!has_mito_prefix) {
  if (file.exists(mito_path)) {
    mito_genes <- read.delim(mito_path, sep = "\t", header = FALSE)
    mito_genes <- c(t(mito_genes))
  } else {
    if (species == "human"){
      stop("❌ Mitochondrial gene list not found at: ", mito_path)
    }
  }
}

for(i in seq_along(SeuratList)){
  if (species == "human"){
    SeuratList[[i]][["percent.mt"]] <- PercentageFeatureSet(SeuratList[[i]], features = mito_genes)
  } else {
    SeuratList[[i]][["percent.mt"]] <- PercentageFeatureSet(SeuratList[[i]], pattern = "^mt-")
    
  }
}

seurat_filtered_list <- list() 
all_qc_metadata <- list()
for (sample in seq_along(SeuratList)) {
  seurat_obj <- SeuratList[[sample]]
  sce <- as.SingleCellExperiment(seurat_obj,assay = "RNA")
  sample_id <- unique(seurat_obj$sample)
  
  # Detect mitochondrial genes
  if (! has_mito_prefix) {
    mito_flag <- rownames(sce) %in% mito_genes
  } else if (any(grepl("^MT-", rownames(seurat_obj[["RNA"]]))) &
             !any(grepl("^mt-", rownames(seurat_obj[["RNA"]])))
  ) {
    mito_genes_detected <- rownames(seurat_obj[["RNA"]])[grepl("^MT[-]?", rownames(seurat_obj[["RNA"]]))]
    mito_flag <- rownames(sce) %in% mito_genes_detected
  } else {
    mito_genes_detected <- rownames(seurat_obj[["RNA"]])[grepl("^mt[-]?", rownames(seurat_obj[["RNA"]]))]
    mito_flag <- rownames(sce) %in% mito_genes_detected
  }

  # Calculate QC metrics
  qc_metrics <- perCellQCMetrics(sce, subsets = list(Mt = mito_flag))
  
  # Apply QC thresholds
  qc_metrics$low_lib <- isOutlier(qc_metrics$sum, log = TRUE, type = "lower", nmads = 5)
  qc_metrics$low_feats <- isOutlier(qc_metrics$detected, log = TRUE, type = "lower", nmads = 5) |
  qc_metrics$sum < thresholds_df$cutoff[thresholds_df$sample == sample_id]
  qc_metrics$high_mito <- qc_metrics$subsets_Mt_percent > 20
  qc_metrics$qc_pass <- !(qc_metrics$low_feats | qc_metrics$high_mito)
  
  # Add metadata
  qc_metrics <- as.data.frame(qc_metrics)
  qc_metrics <- qc_metrics[colnames(seurat_obj), , drop = FALSE]
  seurat_obj <- AddMetaData(seurat_obj, metadata = qc_metrics)
  seurat_obj$qc_status <- ifelse(seurat_obj$qc_pass, "Pass", "Fail")
  
  # Subset Seurat object
  meta_pass <- seurat_obj@meta.data %>% filter(qc_pass)
  seurat_filtered <- subset(seurat_obj, cells = rownames(meta_pass))
  seurat_filtered_list[[sample]] <- seurat_filtered
  
  # Store for combined plots
  qc_metrics$sample_id <- sample_id
  qc_metrics$nFeature_RNA <- seurat_obj$nFeature_RNA[rownames(qc_metrics)]
  qc_metrics$percent_mt <- qc_metrics$subsets_Mt_percent
  qc_metrics$qc_status <- seurat_obj$qc_status[rownames(qc_metrics)]
  all_qc_metadata[[sample]] <- qc_metrics
  
  # Print summary
  message("\nQC summary for ", sample_id, "\n")
  print(qc_metrics %>% summarise(
    total_cells = n(),
    kept = sum(qc_pass),
    percent_kept = mean(qc_pass) * 100,
    low_feats = sum(low_feats),
    high_mito = sum(high_mito),
    low_lib = sum(low_lib)
  ))
}

## reassign sample name to slices of seurat obj list:
names(seurat_filtered_list) <- names(SeuratList)


# Combine all metadata
qc_df <- do.call(rbind, all_qc_metadata)

# Plot combined histogram: nFeature_RNA
png(filename = file.path(plotQC_dir, "QC_pass_fail_combined_nReads_nFeatures.png"), width = 12, height = 8, res = 300, units = "in")
p1 <- ggplot(qc_df, aes(x = detected, fill = qc_pass)) +
  geom_histogram(alpha = 0.5, bins = 100, position = "identity") +
  scale_x_log10() +
  theme_minimal() +
  facet_wrap(~sample_id, scales = "free_y") +
  labs(title = "nFeature_RNA: Pass vs Fail by Sample",
       x = "log10(Number of detected genes)", y = "Cell count")
print(p1)
dev.off()

png(filename = file.path(plotQC_dir, "QC_combined_MTcounts_FeatureScatter.png"), width = 10, height = 7, res = 300, units = "in")
p2 <- ggplot(qc_df, aes(x = nFeature_RNA, y = percent_mt, color = qc_pass)) +
  geom_point(alpha = 0.3, size = 0.5) +
  facet_wrap(~sample_id) +
  theme_minimal() +
  labs(title = "nFeature_RNA vs Mito Percent (by Sample)",
       x = "Number of detected genes", y = "Percent mitochondrial reads")
print(p2)
dev.off()

for (sample_name in names(seurat_filtered_list)) {
  seurat_obj <- seurat_filtered_list[[sample_name]]
  
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, features=VariableFeatures(seurat_obj))
  
  # Store back in processed list
  seurat_filtered_list[[sample_name]] <- seurat_obj
  print(paste0("✅ Normalised and scaled: ", sample_name))
}

saveRDS(seurat_filtered_list, file = file.path(cache_dir, paste0("NormAndScaled_VU40T_Seurat_filtered_individual_samples_list_", species,".rds")))

message(paste0("Normalized and scaled RDS QC'd checkpoint saved to: ", 
             file.path(cache_dir, paste0("NormAndScaled_VU40T_Seurat_filtered_individual_samples_list_", species,".rds"))))
      
