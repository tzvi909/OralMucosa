### Doubletfinder R

library(DoubletFinder)
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)
library(SingleCellExperiment)
library(optparse)


### functs

# Find significant PCs
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

run_doubletfinder_custom <- function(seu_sample_subset, pc_summary, multiplet_rate = NULL, nCores = 1){
  # for debug
  #seu_sample_subset <- samp_split[[1]]
  # Print sample number
  
  sample_name = unique(seu_sample_subset[['sample']]$sample)
  print(paste0("Sample ", sample_name, '...........')) 
  
  if(is.null(multiplet_rate)){
    print('multiplet_rate not provided....... estimating multiplet rate from cells in dataset')
    
    # 10X multiplet rates table
    #https://rpubs.com/kenneditodd/doublet_finder_example
    multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                      'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
                                      'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))
    
    print(multiplet_rates_10x)
    
    multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(seu_sample_subset@meta.data)) %>% 
      dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
      dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
    
    print(paste('Setting multiplet rate to', multiplet_rate))
  }
  
  # pK identification (no ground-truth) 
  #introduces artificial doublets in varying props, merges with real data set and 
  # preprocesses the data + calculates the prop of artficial neighrest neighbours, 
  # provides a list of the proportion of artificial nearest neighbours for varying
  sample <- seu_sample_subset ## in this case as i've done the preprocessing above
  min_pc <- pc_summary$optimal_pcs[pc_summary$sample == sample_name]
  
  # combinations of the pN and pK
  sweep_list <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE, num.cores = nCores)   
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats) # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  optimal.pk <- bcmvn %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
  
  ## Homotypic doublet proportion estimate
  annotations <- sample@meta.data$seurat_clusters # use the clusters as the user-defined cell types
  homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets
  
  nExp.poi <- round(multiplet_rate * nrow(sample@meta.data)) # multiply by number of cells to get the number of expected multiplets
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
  
  # run DoubletFinder
  sample <- doubletFinder(seu = sample, 
                          PCs = 1:min_pc, 
                          pK = optimal.pk, # the neighborhood size used to compute the number of artificial nearest neighbours
                          nExp = nExp.poi.adj) # number of expected real doublets
  # change name of metadata column with Singlet/Doublet information
  colnames(sample@meta.data)[grepl('DF.classifications.*', colnames(sample@meta.data))] <- "doublet_finder"
  
  # Subset and save
  # head(sample@meta.data['doublet_finder'])
  # singlets <- subset(sample, doublet_finder == "Singlet") # extract only singlets
  # singlets$ident
  double_finder_res <- sample@meta.data['doublet_finder'] # get the metadata column with singlet, doublet info
  double_finder_res <- tibble::rownames_to_column(double_finder_res, "row_names") # add the cell IDs as new column to be able to merge correctly
  return(double_finder_res)
}

option_list <- list(
  make_option(c("-s", "--species"),
              type = "character",
              default = "human",
              help = "species to use [default = %default]. Options: human or mouse",
              metavar = "character"),
  make_option(c("-i", "--input"),
              type = "character",
              default = ,
              help = "input RDS dir to load",
              metavar = "character"),
  make_option(c("-n", "--n_cores"),
              type = "integer",
              default = 1,
              help = "Number of cores to use to run  [default = %default]",
              metavar = "character"),
  make_option(c("-c", "--cachedir"),
              type = "character",
              default = 1,
              help = "Path to save intermediate RDS caches [default = %default]",
              metavar = "character"),
  make_option(c("-o", "--outdir"),
              type = "character",
              default = 1,
              help = "Path to save results [default = %default]",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
print(opt)

n_cores <- opt$n_cores
cache_dir <- opt$cachedir
outdir <- opt$outdir
species <- opt$species
input <- opt$input

# Store selected species
species <- tolower(opt$species)

# Validate species input
if (!species %in% c("human", "mouse")) {
  stop("Invalid species. Use 'human' or 'mouse'.")
}

### dir setup

out_prefix_dir <- file.path(outdir, "Seperate_samples")

proj_dir <- "/rds/projects/g/gendood-3dmucosa/"
analysis_dir <- file.path(proj_dir, "scRNAseqAnalysis/")
git_dir <- file.path(analysis_dir, "OralMucosa/VU40T_analysis")
out_prefix_dir <- file.path(git_dir, "Seperate_samples")
cache_dir <- file.path(proj_dir, "rds_cache")


## check for dirs recursively

chk_dir_list <- list(out_prefix_dir, cache_dir)

for (path in chk_dir_list){
  if(!(dir.exists(path))){
    dir.create(path, recursive = T)
  }
}
if (is.null(input)){
  seurat_filtered_list <- readRDS(file.path(cache_dir, paste0("NormAndScaled_VU40T_Seurat_filtered_individual_samples_list_", species,".rds")))
  } else {
  seurat_filtered_list <- input
  ### sanity check if using opt$input arg to make sure species matches seurat list RDS
  if (any(sapply(seurat_filtered_list, function(obj) any(grepl("^ENSMUS", rownames(obj[["RNA"]])))))) {
    species <- "mouse"
  } else {
    species <- "human"
  }
}

plots_dir <- file.path(out_prefix_dir, paste0(stringr::str_to_title(species), "/Plots/DoubletFinder"))
res_dir <- file.path(out_prefix_dir, paste0(stringr::str_to_title(species), "/Doublet_detection"))

chk_dir_list <- list(plots_dir, res_dir)

for (path in chk_dir_list){
  if(!(dir.exists(path))){
    dir.create(path, recursive = T)
  }
}


pc_summary <- data.frame(sample = character(), optimal_pcs = integer(), stringsAsFactors = FALSE)

for (sample_name in names(seurat_filtered_list)) {
  seurat_obj <- seurat_filtered_list[[sample_name]]
  optimal_pc <- get_optimal_pc_count(seurat_obj)
  
  pc_summary <- rbind(pc_summary, data.frame(
    sample = sample_name,
    optimal_pcs = optimal_pc
  ))
}
print(pc_summary)
write.csv(pc_summary, file = file.path(res_dir, "pc_summary.csv"), row.names = F)

# dimUsed <- 30

for (i in seq_along(seurat_filtered_list)){
  sample_name = names(seurat_filtered_list)[i]
  seurat_filtered_list[[i]] <- RunUMAP(seurat_filtered_list[[i]], dims=1:pc_summary$optimal_pcs[pc_summary$sample == sample_name], seed.use=666, verbose=FALSE)
}

for (i in seq_along(seurat_filtered_list)){
  sample_name = names(seurat_filtered_list)[i]
  seurat_filtered_list[[i]] <- FindNeighbors(seurat_filtered_list[[i]], 
                                             dims = 1:pc_summary$optimal_pcs[pc_summary$sample == sample_name]
  )
  seurat_filtered_list[[i]] <- FindClusters(seurat_filtered_list[[i]], resolution=0.1)
}


for (i in seq_along(seurat_filtered_list)){
  obj <- seurat_filtered_list[[i]]
  sample_name <- names(seurat_filtered_list)[i]
  png(file = file.path(plots_dir, paste0("preDoubletFinder_UMAP_", sample_name, ".png"))
      , width = 10, height = 5, units = "in", res = 300)

  p <- DimPlot(obj, reduction = "umap", label = TRUE)
  print(p)
  
  dev.off()
  
  
}

##checkpoint

saveRDS(seurat_filtered_list, file = file.path(cache_dir, paste0("preDoubletFinder_VU40T_seuratList_",species,".rds")))


## DoubletFinder


##run if crash
seurat_filtered_list <- readRDS(file.path(cache_dir, paste0("preDoubletFinder_VU40T_seuratList_",species,".rds")))

pc_summary <- read.csv(file.path(res_dir, "pc_summary.csv"))




for (i in seq_along(seurat_filtered_list)) {
  sample_name <- names(seurat_filtered_list)[i]
  seurat_obj <- seurat_filtered_list[[i]]
  
  message("Running DoubletFinder on ", sample_name)
  
  # Run the custom function
  doublet_meta <- run_doubletfinder_custom(
    seurat_obj, 
    pc_summary = pc_summary, 
    nCores = n_cores ##run with 1 for Rstudio. change when run as standalone
  )
  
  # Merge result back into Seurat object
  doublet_meta <- tibble::column_to_rownames(doublet_meta, var = "row_names")
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, doublet_meta[rownames(seurat_obj@meta.data), , drop = FALSE])
  
  # Save back into list
  seurat_filtered_list[[i]] <- seurat_obj
}


## delete duplicate doubletfinder cols

seurat_filtered_list <- lapply(seurat_filtered_list, function(seu) {
  seu@meta.data <- seu@meta.data[, !duplicated(colnames(seu@meta.data), fromLast = TRUE)]
  return(seu)
})

saveRDS(seurat_filtered_list, file = file.path(cache_dir, paste0("DoubletFinder_results_VU40T_",species,".rds")))

seurat_filtered_list <- readRDS(file.path(cache_dir, paste0("DoubletFinder_results_VU40T_",species,".rds")))


# Combine metadata from all samples
doublets_summary <- purrr::map2_dfr(
  seurat_filtered_list,
  names(seurat_filtered_list),
  ~ {
    meta_df <- .x@meta.data
    
    # If "SampleID" already exists, remove it to avoid name clash
    if ("SampleID" %in% colnames(meta_df)) {
      meta_df <- meta_df[, colnames(meta_df) != "SampleID"]
    }
    
    meta_df %>%
      mutate(SampleID = .y) %>%
      group_by(SampleID, doublet_finder) %>%
      summarise(total_count = n(), .groups = "drop")
  }
) %>%
  group_by(SampleID) %>%
  mutate(countT = sum(total_count)) %>%
  mutate(percent = paste0(round(100 * total_count / countT, 2), "%")) %>%
  dplyr::select(-countT) %>%
  ungroup()
# Save to file
write.table(
  doublets_summary,
  file = file.path(res_dir, paste0("VU40T_doubletfinder_doublets_summary.txt")),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

### DoubletFinder Val


for (i in seq_along(seurat_filtered_list)){
  sample_name <- names(seurat_filtered_list)[i]
  p <- DimPlot(seurat_filtered_list[[i]], reduction = "umap", group.by = "doublet_finder", label = TRUE)
  
  png(file = file.path(plots_dir, paste0("DoubletFinderMarked_UMAP_", sample_name, ".png")), 
  width = 10, height = 5, units = "in", res = 300)
  
  print(p +
        patchwork::plot_annotation(title = paste0("DoubletFinder Classification on UMAP ", sample_name))
  )
  dev.off()
}


## get only singlets
seurat_singlets_list <- list()
for (i in seq_along(seurat_filtered_list)){
  seurat_singlets_list[[i]] <- subset(seurat_filtered_list[[i]], subset = doublet_finder == "Singlet")
  names(seurat_singlets_list)[i] <- names(seurat_filtered_list)[i]
}

saveRDS(seurat_singlets_list, file = file.path(cache_dir, paste0("VU40T_singlets_only_sep_samples_",species,".RDS")))

message("Doublets removed from samples successfully \n
        RDS checkpoint saved to: ", file.path(cache_dir, paste0("VU40T_singlets_only_sep_samples_",species,".RDS")))

