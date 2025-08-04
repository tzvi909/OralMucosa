library(Seurat)
library(ggplot2)
library(readxl)
library(ggridges)
library(dplyr)

git_dir <- "~/OralMucosa/VU40T_analysis"


## to avoid hitting that 20gb quota on my home dir,
proj_dir <- "/rds/projects/g/gendood-3dmucosa/"
cache_dir <- file.path(proj_dir, "rds_cache")

## check for cache dir
if(!(dir.exists(cache_dir))){
  dir.create(cache_dir, recursive = T)
}

# ### human data
# VU40T.combined <- readRDS(file.path(cache_dir, "VU40T_combined_joined_0.3_res_humanOnly.rds"))
# 
# ##mouse only
# VU40T.combined <- readRDS(file.path(cache_dir, "VU40T_combined_joined_0.6_res_mouseOnly.rds"))

### GSEA_data
GSEA_list <- list()

for (i in seq(1,9)){
  GSEA_list[[i]] <- as.data.frame(read_excel(file.path(proj_dir, "scRNAseqAnalysis/Reactome_epi_clusters.xlsx"),sheet = i))
}

### modify lookup df
# Extract current column name
colname <- colnames(GSEA_list[[9]])[1]

# Convert to character vector and prepend column name as first row
new_vector <- c(colname, as.character(GSEA_list[[9]][[1]]))

# Create new dataframe with no column name
GSEA_list[[9]] <- data.frame(Term = new_vector, stringsAsFactors = FALSE)

# Initialize an empty list to store grouped results by term
term_matches <- list()

# Loop over each term in the search list
for (term in GSEA_list[[9]]$Term) {
  # Initialize empty list to collect matches from all clusters
  matches <- list()
  
  # Loop over the 8 clusters
  for (i in seq(1, 8)) {
    df <- GSEA_list[[i]]
    
    match_df <- df[df$Term == term, ]
    
    if (nrow(match_df) > 0) {
      match_df$Cluster <- paste0("Cluster_", i)  # Optional: add cluster info
      matches[[length(matches) + 1]] <- match_df
    }
  }
  
  # Combine matches into one dataframe if any were found
  if (length(matches) > 0) {
    term_matches[[term]] <- do.call(rbind, matches)
  }
}
term_matches[[1]]

# Convert list of data frames to one long data frame
all_terms_df <- dplyr::bind_rows(
  lapply(names(term_matches), function(term) {
    df <- term_matches[[term]]
    df$Term_Full <- term  # Add the term as a column
    df
  }),
  .id = "Term_ID"
)

plot_df <- all_terms_df %>%
  select(Term = Term_Full, Cluster, Fold_Enrichment = `Fold Enrichment`, PValue) %>%
  filter(!is.na(Fold_Enrichment), !is.na(PValue))

plot_df$Cluster <- sapply(plot_df$Cluster, function(cl) {
  cl_num <- as.numeric(sub("Cluster_", "", cl))
  paste0(cl_num - 1)
})


cluster_colors <- c(
  "0" = "#E69F00",  # orange
  "1" = "#56B4E9",  # sky blue
  "2" = "#009E73",  # bluish green
  "3" = "#F0E442",  # yellow
  "4" = "#0072B2",  # blue
  "5" = "#D55E00",  # vermillion
  "6" = "#CC79A7",  # reddish purple
  "7" = "#999999"   # grey
)


png(file.path(git_dir, "Integrated/Plots/HumanOnly_GSEA_cluster_Dotplot_fold_enrich_onXaxis.png")
    , width = 12, height = 5, units = "in", res = 300)
p <- ggplot(plot_df, aes(x = Fold_Enrichment, y = Term, color = Cluster, size = -log10(PValue))) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "Reactome Term",
    size = expression(-log[10](P-value)),
    title = "Cluster-level Pathway Enrichment"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + scale_color_manual(values = cluster_colors)
print(p)
dev.off()

png(file.path(git_dir, "Integrated/Plots/HumanOnly_GSEA_cluster_Dotplot_log_p_onXaxis.png")
    , width = 10, height = 5, units = "in", res = 300)
p <- ggplot(plot_df, aes(x = -log10(PValue), y = Term, color = Cluster, size = Fold_Enrichment)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = expression(-log[10](P-value)),
    y = "Reactome Term",
    size = "Fold Enrichment",
    title = "Cluster-level Pathway Enrichment"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + scale_color_manual(values = cluster_colors)
print(p)
dev.off()
