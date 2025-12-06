#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              help = "Input RDS directory to load"),
  make_option(c("-n", "--n_cores"),
              type = "integer",
              default = 1,
              help = "Number of cores to use"),
  make_option(c("-c", "--cachedir"),
              type = "character",
              help = "Directory to save intermediate RDS caches"),
  make_option(c("-o", "--outdir"),
              type = "character",
              help = "Directory to save results")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required args
stopifnot(!is.null(opt$input), !is.null(opt$cachedir), !is.null(opt$outdir))

# -------------------------------------------------------------------
# Helper: temporarily simulate command-line args for sourced scripts
# -------------------------------------------------------------------
set_script_args <- function(args_vec) {
  assign(".commandArgs", c("Rscript", args_vec), envir = .GlobalEnv)
}

get_args <- function(species = NULL) {
  c(
    if (!is.null(species)) c("--species", species) else NULL,
    "--input", opt$input,
    "--n_cores", opt$n_cores,
    "--cachedir", opt$cachedir,
    "--outdir", opt$outdir
  )
}

run_script <- function(script, species = NULL) {
  message("\n---- Running ", script,
          if (!is.null(species)) paste0(" [", species, "]") else "",
          " ----")
  
  set_script_args(get_args(species))
  source(script, local = FALSE)
}

# -------------------------------------------------------------------
# Scripts grouped by whether they take the species argument
# -------------------------------------------------------------------

scripts_with_species <- c(
  "species_mixing.R",
  "pre-process_QC.R",
  "DoubletFinder.R",
  "pre-integration_seurat.R",
  "Integration_seurat.R",
  "epithelial_pseudotime.R",
  "fibroblast_pseudotime.R"
)

scripts_without_species <- c(
  "GSEA_plotting.R",  # no species
  "cellchat.R"        # no species
)

# -------------------------------------------------------------------
# Run pipeline for both human & mouse
# -------------------------------------------------------------------

for (sp in c("human", "mouse")) {
  
  message("\n====================================")
  message("      STARTING SPECIES: ", sp)
  message("====================================")
  
  # scripts requiring species argument
  for (script in scripts_with_species) {
    run_script(script, species = sp)
  }
  
  # scripts that do NOT require species
  run_script("GSEA_plotting.R", species = NULL)
  
  message("========== COMPLETED ", sp, " ==========\n")
}

# cellchat runs ONCE only (still species-independent)
run_script("cellchat.R", species = NULL)

message("All pipeline runs completed successfully.")
