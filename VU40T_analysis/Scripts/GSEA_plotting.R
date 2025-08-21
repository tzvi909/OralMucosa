library(Seurat)
library(ggplot2)
library(readxl)
library(dplyr)

git_dir <- "~/OralMucosa/VU40T_analysis"


## to avoid hitting that 20gb quota on my home dir,
proj_dir <- "/rds/projects/g/gendood-3dmucosa/"
cache_dir <- file.path(proj_dir, "rds_cache")

## check for cache dir
if(!(dir.exists(cache_dir))){
  dir.create(cache_dir, recursive = T)
}


### GSEA_data_human
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

## filter by top 10 terms and sig p-val

for (i in seq(1, 8)) {
  GSEA_list[[i]] <- GSEA_list[[i]] %>%
    dplyr::filter(FDR <= 0.05) %>%
    dplyr::arrange(PValue) # %>%
  # dplyr::slice_head(n = 10)
}

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

custom_order <- c(
  "R-HSA-6805567~Keratinization",
  "R-HSA-421270~Cell-cell junction organization",
  "R-HSA-168256~Immune System",
  "R-HSA-449147~Signaling by Interleukins",
  "R-HSA-446107~Type I hemidesmosome assembly",
  "R-HSA-3000157~Laminin interactions",
  "R-HSA-1474290~Collagen formation",
  "R-HSA-1474244~Extracellular matrix organization",
  "R-HSA-9006934~Signaling by Receptor Tyrosine Kinases",
  "R-HSA-8953897~Cellular responses to stimuli",
  "R-HSA-9755511~KEAP1-NFE2L2 pathway",
  "R-HSA-156842~Eukaryotic Translation Elongation",
  "R-HSA-9706574~RHOBTB GTPase Cycle",
  "R-HSA-9716542~Signaling by Rho GTPases, Miro GTPases and RHOBTB3",
  "R-HSA-69278~Cell Cycle, Mitotic"
)

plot_df$Term <- factor(plot_df$Term, levels = rev(custom_order))
plot_df$Term <- sub(".*~", "", plot_df$Term)

plot_df$Term

plot_df$Term <- mgsub::mgsub(plot_df$Term, 
                             c("Extracellular matrix", "Receptor Tyrosine Kinases", "hemidesmosome"), 
                             c("ECM", "RTKs", "HD")
)
plot_df$Term
plot_df$Term <- sapply(plot_df$Term, function(x) {
  if (nchar(x) > 25 && grepl(" ", x)) {
    sub("^(.{25,}?)(?:\\s)", "\\1\n", x, perl = TRUE)
  } else {
    x
  }
})
plot_df$Term <- sapply(plot_df$Term, function(x) {
  if (nchar(x) > 35 && grepl(" ", x) && grepl("\n", x)) {
    sub("^(.{35,}?)(?:\\s)", "\\1<br>", x, perl = TRUE)
  } else {
    x
  }
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

## different plotting options

png(file.path(git_dir, "Integrated/Plots/HumanOnly_GSEA_cluster_Dotplot_fold_enrich_onXaxis.png")
    , width = 7, height = 5, units = "in", res = 300)
p <- ggplot(plot_df, aes(x = Fold_Enrichment, y = Term, color = Cluster, size = -log10(PValue))) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "",
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
    , width = 7, height = 5, units = "in", res = 300)
p <- ggplot(plot_df, aes(x = -log10(PValue), y = Term, color = Cluster, size = Fold_Enrichment)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = expression(-log[10](P-value)),
    y = "",
    size = "Fold Enrichment",
    title = "Cluster-level Pathway Enrichment"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + scale_color_manual(values = cluster_colors)
print(p)
dev.off()

png(file.path(git_dir, "Integrated/Plots/HumanOnly_GSEA_cluster_Dotplot_clusts_onXaxis.png")
    , width = 7, height = 5, units = "in", res = 300)
p <- ggplot(plot_df, aes(x = Cluster , y = Term, color = -log10(PValue) , size = Fold_Enrichment)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = "Cluster",
    y = "",
    size = "Fold Enrichment",
    color = expression(-log[10](P-value)),
    title = "Cluster-level Pathway Enrichment"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 0, hjust = 1)
  ) + scale_color_gradient(low = "yellow", high = "red")
print(p)
dev.off()

### WP HN Cancer GSEAs


GSEA_list <- list()

for (i in seq(1,9)){
  GSEA_list[[i]] <- as.data.frame(read_excel(file.path(proj_dir, "scRNAseqAnalysis/Wiki_epi_clusters.xlsx"),sheet = i))
}


# Create new dataframe with no column name
GSEA_list[[9]] <- data.frame(Term = GSEA_list[[9]], stringsAsFactors = TRUE)

## filter by top 10 terms and sig p-val

for (i in seq(1, 8)) {
  GSEA_list[[i]] <- GSEA_list[[i]] %>%
    dplyr::filter(FDR <= 0.05) %>%
    dplyr::arrange(PValue) # %>%
  # dplyr::slice_head(n = 10)
}

# Initialize an empty list to store grouped results by term
term_matches <- list()

# Loop over each term in the search list
for (term in GSEA_list[[9]][,1]) {
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

plot_df$Term <- factor(plot_df$Term, levels = rev(GSEA_list[[9]][,1]))
plot_df$Term <- sub(".*:", "", plot_df$Term)
## shorten some of these terms for the plot
plot_df$Term <- mgsub::mgsub(plot_df$Term, 
                             c("TGF beta receptor", "Alpha 6 beta 4", " beta", "Retinoblastoma gene"), 
                             c("TGFBR", "&alpha;<sub>6</sub>&beta;<sub>4</sub>", "&hyphen;&beta;", "RB")
)
##using html now so wrap with <br> not \n.
## do not recognise html <sub> in char count
plot_df$Term <- sapply(plot_df$Term, function(x) {
  if (nchar(x) > 29 && grepl(" ", x) && !(grepl("<sub>", x))  && !(grepl("&beta;", x))) {
    sub("^(.{15,}?)(?:\\s)", "\\1<br>", x, perl = TRUE)
  } else {
    x
  }
})
plot_df$Term <- sapply(plot_df$Term, function(x) {
  if (nchar(x) > 35 && grepl(" ", x) && grepl("<br>", x) && !(grepl("<sub>", x))) {
    sub("^(.{35,}?)(?:\\s)", "\\1<br>", x, perl = TRUE)
  } else {
    x
  }
})


plot_df$Term <- factor(plot_df$Term, levels = rev(unique(plot_df$Term)))

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

## different plotting options

png(file.path(git_dir, "Integrated/Plots/HumanOnly_HN_Cancer_GSEA_cluster_Dotplot_fold_enrich_onXaxis.png")
    , width = 6, height = 6.3, units = "in", res = 300)
p <- ggplot(plot_df, aes(x = Fold_Enrichment, y = Term, color = Cluster, size = -log10(PValue))) +
  geom_point(alpha = 0.8) +
  # scale_y_discrete(labels = identity) +
  theme_minimal() +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "",
    size = expression(-log[10](P-value)),
    title = "Head and Neck Cancer Cluster-level\nPathway Enrichment"
  ) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  ) + scale_color_manual(values = cluster_colors)
print(p)
dev.off()

png(file.path(git_dir, "Integrated/Plots/HumanOnly_HN_Cancer_GSEA_cluster_Dotplot_log_p_onXaxis.png")
    , width = 6, height = 6.3, units = "in", res = 300)
p <- ggplot(plot_df, aes(x = -log10(PValue), y = Term, color = Cluster, size = Fold_Enrichment)) +
  geom_point(alpha = 0.8) +
  scale_y_discrete(labels = identity) +
  theme_minimal() +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = expression(-log[10](P-value)),
    y = "",
    size = "Fold Enrichment",
    title = "Head and Neck Cancer Cluster-level\nPathway Enrichment"
  ) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  ) + scale_color_manual(values = cluster_colors)
print(p)
dev.off()

png(file.path(git_dir, "Integrated/Plots/HumanOnly_HN_Cancer_GSEA_cluster_Dotplot_clusts_onXaxis.png")
    , width = 6, height = 6.3, units = "in", res = 300)
p <- ggplot(plot_df, aes(x = Cluster , y = Term, color = -log10(PValue) , size = Fold_Enrichment)) +
  geom_point(alpha = 0.8) +
  scale_y_discrete(labels = identity) +
  theme_minimal() +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = "Cluster",
    y = "",
    size = "Fold Enrichment",
    color = expression(-log[10](P-value)),
    title = "Head and Neck Cancer Cluster-level\nPathway Enrichment"
  ) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 8),
    axis.text.x = element_text(angle = 0, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  ) + scale_color_gradient(low = "yellow", high = "red")
print(p)
dev.off()

