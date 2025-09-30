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

### dir setup

proj_dir <- "/rds/projects/g/gendood-3dmucosa/"
analysis_dir <- file.path(proj_dir, "scRNAseqAnalysis/")
git_dir <- file.path(analysis_dir, "OralMucosa/VU40T_analysis")
out_prefix_dir <- file.path(git_dir, "OralMucosa/Seperate_samples")

cache_dir <- file.path(proj_dir, "rds_cache")

## check for dirs recursively

chk_dir_list <- list(analysis_dir, git_dir, cache_dir)

for (path in chk_dir_list){
  if(!(dir.exists(path))){
    dir.create(path, recursive = T)
  }
}

option_list <- list(
  make_option(c("-s", "--species"),
              type = "character",
              default = "human",
              help = "Species to use [default = %default]. Options: human or mouse",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Store selected species
species <- tolower(opt$species)

# Validate species input
if (!species %in% c("human", "mouse")) {
  stop("Invalid species. Use 'human' or 'mouse'.")
}

if (species == "human"){
  mtx_dir <- file.path(proj_dir, "BaseSpace/LPS_VU40T_QC_and_counts/cellranger/mtx_conversions")
  plotQC_dir <- file.path(out_prefix_dir, "Plots/QC")
} else {
  mtx_dir <- file.path(proj_dir, "BaseSpace/LPS_VU40T_QC_and_counts_Mouse/cellranger/mtx_conversions")
  plotQC_dir <- file.path(out_prefix_dir, "Mouse/Plots/QC")
}


Samples <- list.files(path = mtx_dir, 
                      pattern = "^P",
                      recursive = F
) ##check we only have 4 samples
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

# Replace SeuratList names with simplified format like P18_LPS-N, P22_LPS-P, etc.
names(SeuratList) <- gsub("^.*/(P[0-9]+)_LPS-([NP]).*$", "\\1_LPS-\\2", names(SeuratList))

## Mito genes were isolated from hg38 GTF (chrM)

if(species == "human"){
  mito_genes <- read.delim("/rds/projects/g/gendood-3dmucosa/BaseSpace/LPS_VU40T_QC_and_counts/cellranger/mkref/cellranger_reference/genes/MTgenes.txt", sep = "\t", header = F)
  
  ## convert to vector
  mito_genes <- c(t(mito_genes))
}
thresholds_df <- data.frame(
  sample = c("P18_LPS-N", "P18_LPS-P", "P22_LPS-N", "P22_LPS-P"),
  cutoff = c(300, 350, 300, 300)
)

seurat_filtered_list <- list() 
for (i in seq_along(SeuratList)) {
  seurat_obj <- SeuratList[[i]]
  sce <- as.SingleCellExperiment(seurat_obj)
  
  # Detect mitochondrial genes
  if (!any(grepl("^MT", rownames(seurat_obj[["RNA"]])))) {
    mito_flag <- rownames(sce) %in% mito_genes
  } else {
    mito_genes_detected <- rownames(seurat_obj[["RNA"]])[grepl("^MT[-]?", rownames(seurat_obj[["RNA"]]))]
    mito_flag <- rownames(sce) %in% mito_genes_detected
  }

  
  # Calculate QC metrics
  qc_metrics <- perCellQCMetrics(sce, subsets = list(Mt = mito_flag))
  
  # Apply QC thresholds
  qc_metrics$low_lib <- isOutlier(qc_metrics$sum, log = TRUE, type = "lower", nmads = 5)
  qc_metrics$low_feats <- isOutlier(qc_metrics$detected, log = TRUE, type = "lower", nmads = 5) |
    qc_metrics$sum < thresholds_df$cutoff[thresholds_df$sample == names(SeuratList)[i]]
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
  seurat_filtered_list[[sample_id]] <- seurat_filtered
  
  # Store for combined plots
  qc_metrics$sample_id <- sample_id
  qc_metrics$nFeature_RNA <- seurat_obj$nFeature_RNA[rownames(qc_metrics)]
  qc_metrics$percent_mt <- qc_metrics$subsets_Mt_percent
  qc_metrics$qc_status <- seurat_obj$qc_status[rownames(qc_metrics)]
  all_qc_metadata[[i]] <- qc_metrics
  
  # Print summary
  cat("\nQC summary for", sample_id, "\n")
  print(qc_metrics %>% summarise(
    total_cells = n(),
    kept = sum(qc_pass),
    percent_kept = mean(qc_pass) * 100,
    low_feats = sum(low_feats),
    high_mito = sum(high_mito),
    low_lib = sum(low_lib)
  ))
}


# Combine all metadata
qc_df <- do.call(rbind, all_qc_metadata)

# Plot combined histogram: nFeature_RNA
png(filename = file.path(plotQC_dir, "QC_pass_fail_combined_nReads_nFeatures.png"), width = 12, height = 8, res = 300)
p1 <- ggplot(qc_df, aes(x = detected, fill = qc_pass)) +
  geom_histogram(alpha = 0.5, bins = 100, position = "identity") +
  scale_x_log10() +
  theme_minimal() +
  facet_wrap(~sample_id, scales = "free_y") +
  labs(title = "nFeature_RNA: Pass vs Fail by Sample",
       x = "log10(Number of detected genes)", y = "Cell count")
print(p1)
dev.off()

png(filename = file.path(plotQC_dir, "QC_combined_FeatureScatter.png"), width = 10, height = 7, res = 300)
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
  
  # Standard Seurat preprocessing
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  
  # Store back in processed list
  seurat_filtered_list[[sample_name]] <- seurat_obj
  
  cat(paste0("✅ Normalized and scaled: ", sample_name, "\n"))
}

saveRDS(seurat_filtered_list, file = file.path(cache_dir, paste0("NormAndScaled_VU40T_Seurat_filtered_individual_samples_list_", species,".rds")))