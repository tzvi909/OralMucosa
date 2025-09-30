library(Seurat)
library(ggplot2)
library(readxl)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)

### funcs

### for the fibroblast GSEAs
# ---- Parse a single sheet with DAVID "Annotation Cluster" blocks ----
parse_david_sheet <- function(path, sheet) {
  raw <- read_excel(path, sheet = sheet, col_names = FALSE)
  if (nrow(raw) == 0) return(tibble())
  
  names(raw) <- paste0("col", seq_len(ncol(raw)))
  raw <- mutate(raw, row_id = row_number())
  
  # find cluster header rows
  hdr_rows <- which(str_detect(raw$col1 %||% "", 
                               regex("^\\s*Annotation\\s+Cluster\\s*\\d+", ignore_case = TRUE)))
  if (length(hdr_rows) == 0) return(tibble())
  
  hdr_next <- c(hdr_rows[-1], nrow(raw) + 1)
  
  parse_block <- function(h_idx, h_next) {
    cluster_num <- str_match(raw$col1[h_idx], "Annotation\\s+Cluster\\s*(\\d+)")[,2] %>% as.integer()
    
    # Enrichment score: look across that same row in other columns
    score_cell <- raw[h_idx, -1] %>% unlist(use.names = FALSE) %>% as.character()
    score_val <- str_match(score_cell, "Enrichment\\s*Score\\s*:?\\s*([0-9eE.+-]+)")
    enrich_score <- suppressWarnings(as.numeric(na.omit(score_val[,2])[1]))
    
    # find the column-name row (first "Category" below header)
    rng <- seq.int(h_idx + 1, h_next - 1)
    if (!length(rng)) return(tibble())
    name_row <- rng[which(str_to_upper(str_trim(raw$col1[rng] %||% "")) == "CATEGORY")[1]]
    if (is.na(name_row)) return(tibble())
    
    nm <- raw[name_row, , drop = FALSE] %>% unlist(use.names = FALSE) %>% as.character()
    nm <- nm[!is.na(nm)]
    nm[nm == "%"] <- "Percent"
    
    dat_start <- name_row + 1
    dat_end   <- h_next - 1
    if (dat_start > dat_end) return(tibble())
    
    # Only keep known DAVID columns (up to FDR)
    expected_cols <- c("Category","Term","Count","Percent","PValue","Genes",
                       "List Total","Pop Hits","Pop Total","Fold Enrichment",
                       "Bonferroni","Benjamini","FDR")
    
    nm <- nm[seq_len(min(length(nm), length(expected_cols)))]
    block <- raw[dat_start:dat_end, seq_along(nm), drop = FALSE]
    names(block) <- nm
    
    
    block %>%
      filter(!is.na(Category), str_trim(as.character(Category)) != "") %>%
      mutate(
        sheet              = sheet,
        annotation_cluster = cluster_num,
        enrichment_score   = enrich_score,
        .before = 1
      )
  }
  
  map2_dfr(hdr_rows, hdr_next, parse_block)
}


# ---- Read all sheets, optionally exclude by name pattern, coerce numerics ----
read_david_excel <- function(path,
                             exclude_sheets_pattern = "(?i)lookup|key|map|legend|info") {
  all_sheets <- excel_sheets(path)
  if (!length(all_sheets)) return(tibble())
  
  # exclude obvious non-result sheets by name (customizable)
  keep_sheets <- all_sheets[!str_detect(all_sheets, exclude_sheets_pattern)]
  if (!length(keep_sheets)) keep_sheets <- all_sheets  # fallback if filter removes all
  
  out <- map_dfr(keep_sheets, ~{
    message("Parsing sheet: ", .x)
    dat <- parse_david_sheet(path, .x)
    if (!nrow(dat)) message("  -> skipped (no clusters found)")
    dat
  })
  
  if (!nrow(out)) return(tibble())
  
  # coerce numerics where present
  numeric_cols <- c("Count","Percent","PValue","List Total","Pop Hits","Pop Total",
                    "Fold Enrichment","Bonferroni","Benjamini","FDR",
                    "annotation_cluster","enrichment_score")
  for (cc in intersect(numeric_cols, names(out))) {
    out[[cc]] <- suppressWarnings(as.numeric(out[[cc]]))
  }
  if ("Count" %in% names(out)) out$Count <- suppressWarnings(as.integer(out$Count))
  
  desired_order <- c("sheet","annotation_cluster","enrichment_score",
                     "Category","Term","Count","Percent","PValue","Genes",
                     "List Total","Pop Hits","Pop Total","Fold Enrichment",
                     "Bonferroni","Benjamini","FDR")
  
  out %>%
    relocate(any_of(desired_order), .before = 1) %>%
    distinct()
}


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

## replace signalling with rho, miro gtpases etc. with just R-HSA-194315~Signaling by Rho GTPases
## in lookup tbl
GSEA_list[[9]]$Term[GSEA_list[[9]]$Term == "R-HSA-9716542~Signaling by Rho GTPases, Miro GTPases and RHOBTB3"] <- "R-HSA-194315~Signaling by Rho GTPases"

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
  "R-HSA-194315~Signaling by Rho GTPases",
  "R-HSA-69278~Cell Cycle, Mitotic"
)

plot_df$Term <- factor(plot_df$Term, levels = rev(custom_order))
plot_df$Term <- sub(".*~", "", plot_df$Term)


plot_df$Term <- mgsub::mgsub(plot_df$Term, 
                             c("Extracellular matrix", "Receptor Tyrosine Kinases", "hemidesmosome"), 
                             c("ECM", "RTKs", "HD")
)
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
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  ) + scale_color_gradient(low = "blue", high = "red")
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
  ) + scale_color_gradient(low = "blue", high = "red")
print(p)
dev.off()


### Mouse Fibroblast GSEAs

david_excel <- "~/Pathways_fib_clusters.xlsx"

# 2) Parse everything (skips sheets with names matching the exclude pattern,
#    and any sheet that has no Annotation Cluster blocks - i.e. the lookup tbl on last sheet)
david_tidy <- read_david_excel(david_excel)

# 3) (Optional) filter by FDR; rank by enrichment_score for plotting
david_sig <- david_tidy %>%
  filter(!is.na(FDR), FDR < 0.05) %>%
  arrange(desc(enrichment_score), FDR)

# Create new dataframe with no column name
fibro_gsea_lookup <- as.data.frame(read_excel(david_excel, sheet = 14, col_names = T))[,c(1,2,3,11)]
fibro_gsea_lookup$Term <- factor(fibro_gsea_lookup$Term, levels = unique(fibro_gsea_lookup$Term) )
### filter by lookup table

filtered_df <- david_sig[david_sig$Term %in% fibro_gsea_lookup$Term, ]
colnames(filtered_df)[colnames(filtered_df) == "sheet"] <- "Cluster"


plot_df <- filtered_df %>%
  select(Term = Term, Cluster, Fold_Enrichment = `Fold Enrichment`, PValue, enrichment_score) %>%
  filter(!is.na(Fold_Enrichment), !is.na(PValue))

plot_df$Cluster <- sapply(plot_df$Cluster, function(cl) {
  cl_num <- sub("Cluster_", "", cl)
})

# levels_in_data <- fibro_gsea_lookup$Term[fibro_gsea_lookup$Term %in% plot_df$Term]
# 
# plot_df$Term <- factor(
#   plot_df$Term,
#   levels = rev(c(unique(levels_in_data)))  # or just levels_in_data to drop leftovers
# )
# 
# levels_in_data

plt_str_replace_and_wrap <- function(x){
  x <- sub(".*:", "", x)
  x <- sub(".*~", "", x)
  x <- mgsub::mgsub(x, 
                    c("TGF beta receptor", "Alpha 6 beta 4", " beta", "Retinoblastoma gene", 
                      "Insulin-like Growth Factor"," \\(IGF\\)", 
                      " Binding Proteins", " \\(IGFBPs\\)", 
                      "Toll-like receptor", "tyrosine kinase inhibitor", "Extracellular matrix", "NF-kappa B", "Human papillomavirus",
                      "reactive oxygen species", "Interleukin-1", "Receptor Tyrosine Kinases"), 
                    c("TGFBR", "&alpha;<sub>6</sub>&beta;<sub>4</sub>", "&hyphen;&beta;", "RB",
                      "", "IGF" ,"", "IGFBPs", "TLR", "TKI", "ECM", "NF-kB", "HPV", "ROS", "IL-1",
                      "TKRs"))
  
  x <- sapply(x, function(y) {
    if (nchar(y) > 29 && grepl(" ", y) && !(grepl("<sub>", y))  && !(grepl("&beta;", y))) {
      sub("^(.{15,}?)(?:\\s)", "\\1<br>", y, perl = TRUE)
    } else {
      y
    }
  })
  x <- sapply(x, function(y) {
    if (nchar(y) > 35 && grepl(" ", y) && grepl("<br>", y) && !(grepl("<sub>", y))) {
      sub("^(.{35,}?)(?:\\s)", "\\1<br>", y, perl = TRUE)
    } else {
      y
    }
  })
  
}

plot_df$Term <- plt_str_replace_and_wrap(plot_df$Term)
fibro_gsea_lookup$Term <- plt_str_replace_and_wrap(fibro_gsea_lookup$Term)




plot_df$Term <- factor(plot_df$Term, levels = rev(unique(fibro_gsea_lookup$Term)))

cluster_order <- factor(c(0:12))

# if you already set Term levels somewhere, keep them;
# otherwise capture the current order you want for y:
term_order <- levels(plot_df$Term) %||% unique(plot_df$Term)


cluster_colors <- c(
  "0"  = "#E69F00", # orange
  "1"  = "#56B4E9", # sky blue
  "2"  = "#009E73", # bluish green
  "3"  = "#F0E442", # yellow
  "4"  = "#0072B2", # blue
  "5"  = "#D55E00", # vermillion
  "6"  = "#CC79A7", # reddish purple
  "7"  = "#999999", # grey
  "8"  = "#8DD3C7", # teal
  "9"  = "#FB8072", # coral red
  "10" = "#80B1D3", # light blue
  "11" = "#FDB462", # peach
  "12" = "#B3DE69"  # lime green
)

## different plotting options

png(file.path(git_dir, "Integrated/Mouse/Plots/MouseOnly_GSEA_cluster_Dotplot_fold_enrich_onXaxis.png")
    , width = 10, height = 18, units = "in", res = 300)
p <- ggplot(plot_df, aes(x = Fold_Enrichment, y = Term, color = Cluster, size = -log10(PValue))) +
  geom_point(alpha = 0.8) +
  scale_y_discrete(labels = identity) +
  theme_minimal() +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "",
    size = expression(-log[10](P-value)),
    title = "Cluster-level\nPathway Enrichment"
  ) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  ) + scale_color_manual(values = cluster_colors)
print(p)
dev.off()

png(file.path(git_dir, "Integrated/Mouse/Plots/MouseOnly_GSEA_cluster_Dotplot_log_p_onXaxis.png")
    , width = 10, height = 18, units = "in", res = 300)
p <- ggplot(plot_df, aes(x = -log10(PValue), y = Term, color = Cluster, size = Fold_Enrichment)) +
  geom_point(alpha = 0.8) +
  scale_y_discrete(labels = identity) +
  theme_minimal() +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = expression(-log[10](P-value)),
    y = "",
    size = "Fold Enrichment",
    title = "Cluster-level\nPathway Enrichment"
  ) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  ) + scale_color_manual(values = cluster_colors)
print(p)
dev.off()

png(file.path(git_dir, "Integrated/Mouse/Plots/MouseOnly_GSEA_cluster_Dotplot_clusts_onXaxis.png")
    , width = 10, height = 18, units = "in", res = 300)
p <- ggplot(plot_df, aes(x = Cluster , y = Term, color = -log10(PValue) , size = Fold_Enrichment)) +
  geom_point(alpha = 0.8) +
  scale_y_discrete(labels = identity) +
  scale_x_discrete(limits = cluster_order) +
  theme_minimal() +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = "Cluster",
    y = "",
    size = "Fold Enrichment",
    color = expression(-log[10](P-value)),
    title = "Cluster-level\nPathway Enrichment"
  ) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 8),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title  = element_text(hjust = 0.5)
  ) + scale_color_gradient(low = "blue", high = "red")
print(p)
dev.off()


# Collapse each pathway's genes into a sorted string
identical_sets <- filtered_df %>%
  mutate(
    GeneSet = sapply(strsplit(Genes, ",\\s*"), function(x) paste(sort(unique(x)), collapse = ","))
  ) %>%
  group_by(GeneSet) %>%
  summarise(pathways = list(Term), n = n()) %>%
  filter(n > 1)

identical_sets

compute_cluster_overlaps <- function(df) {
  pathway_genes <- df %>%
    distinct(annotation_cluster, Term, Genes) %>%
    mutate(GeneSet = strsplit(Genes, ",\\s*"))
  
  n <- nrow(pathway_genes)
  
  # if fewer than 2 pathways in the cluster → no pairs
  if (n < 2) return(NULL)
  
  pairs <- t(combn(n, 2)) %>%
    as.data.frame() %>%
    setNames(c("i", "j"))
  
  pairs %>%
    rowwise() %>%
    mutate(
      path1 = pathway_genes$Term[i],
      path2 = pathway_genes$Term[j],
      overlap = length(intersect(pathway_genes$GeneSet[[i]], pathway_genes$GeneSet[[j]])),
      union   = length(union(pathway_genes$GeneSet[[i]], pathway_genes$GeneSet[[j]])),
      jaccard = ifelse(union == 0, 0, overlap / union)
    ) %>%
    ungroup() %>%
    mutate(Cluster = unique(df$Cluster))
}

cluster_overlaps <- filtered_df %>%
  group_split(Cluster) %>%
  map_df(compute_cluster_overlaps)

# only keep pairs that actually overlap
cluster_overlaps <- cluster_overlaps %>%
  filter(overlap > 0)

head(cluster_overlaps, 10)
write.csv(cluster_overlaps, file = file.path(git_dir, "/Integrated/Mouse/pathway_gene_overlaps_byClust.csv"))

library(igraph)

get_redundant_groups <- function(df, threshold = 0.7) {
  edges <- df %>%
    filter(jaccard >= threshold) %>%
    select(path1, path2)
  
  if (nrow(edges) == 0) {
    return(tibble(Term = character(), group = integer()))
  }
  
  g <- graph_from_data_frame(edges, directed = FALSE)
  comps <- components(g)$membership
  
  tibble(
    Term = names(comps),
    group = as.integer(comps)
  )
}


redundant_groups <- cluster_overlaps %>%
  group_by(Cluster) %>%
  group_modify(~ get_redundant_groups(.x, threshold = 0.7)) %>%
  ungroup()

redundant_groups

# Step 3: merge redundancy groups back into the original results
filtered_with_groups <- filtered_df %>%
  left_join(redundant_groups, by = c("Cluster", "Term"))

# Step 4: keep only the most significant pathway per redundancy group
# (here by lowest FDR, but you can switch to highest enrichment score if you want)
filtered_unique <- filtered_with_groups %>%
  group_by(Cluster, group) %>%
  slice_min(FDR, n = 1, with_ties = FALSE) %>%
  ungroup()

# View result
filtered_unique



plot_df <- filtered_unique %>%
  select(Term = Term, Cluster, Fold_Enrichment = `Fold Enrichment`, PValue, enrichment_score) %>%
  filter(!is.na(Fold_Enrichment), !is.na(PValue))

plot_df$Cluster <- sapply(plot_df$Cluster, function(cl) {
  cl_num <- sub("Cluster_", "", cl)
})
plot_df$Term <- plt_str_replace_and_wrap(plot_df$Term)

pruned_order<- intersect(fibro_gsea_lookup$Term, plot_df$Term)

plot_df$Term <- factor(plot_df$Term, levels = rev(pruned_order))

cluster_order <- factor(c(0:12))

# if you already set Term levels somewhere, keep them;
# otherwise capture the current order you want for y:
term_order <- levels(plot_df$Term) %||% unique(plot_df$Term)


png(file.path(git_dir, "Integrated/Mouse/Plots/MouseOnly_GSEA_cluster_pruned_Dotplot_clusts_onXaxis.png")
    , width = 10, height = 10, units = "in", res = 300)
p <- ggplot(plot_df, aes(x = Cluster , y = Term, color = -log10(PValue) , size = Fold_Enrichment)) +
  geom_point(alpha = 0.8) +
  scale_y_discrete(labels = identity) +
  scale_x_discrete(limits = cluster_order) +
  theme_minimal() +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = "Cluster",
    y = "",
    size = "Fold Enrichment",
    color = expression(-log[10](P-value)),
    title = "Cluster-level\nPathway Enrichment"
  ) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 8),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title  = element_text(hjust = 0.5)
  ) + scale_color_gradient(low = "blue", high = "red")
print(p)
dev.off()

new_pruned_termsList <- read.csv("~/new_order_fibro_GSEA.csv", header = F, check.names = F)
new_pruned_termsList <- new_pruned_termsList$V1

filtered_new <- filtered_df[filtered_df$Term %in% new_pruned_termsList,]

new_pruned_termsList<-plt_str_replace_and_wrap(new_pruned_termsList)

plot_df <- filtered_new %>%
  select(Term = Term, Cluster, Fold_Enrichment = `Fold Enrichment`, PValue, enrichment_score) %>%
  filter(!is.na(Fold_Enrichment), !is.na(PValue))

plot_df$Cluster <- sapply(plot_df$Cluster, function(cl) {
  cl_num <- sub("Cluster_", "", cl)
})
plot_df$Term <- plt_str_replace_and_wrap(plot_df$Term)



pruned_order<- intersect(fibro_gsea_lookup$Term, plot_df$Term)

plot_df$Term <- factor(plot_df$Term, levels = rev(new_pruned_termsList))

cluster_order <- factor(c(0:12))

# if you already set Term levels somewhere, keep them;
# otherwise capture the current order you want for y:
term_order <- levels(plot_df$Term) %||% unique(plot_df$Term)


png(file.path(git_dir, "Integrated/Mouse/Plots/MouseOnly_GSEA_cluster_pruned_Dotplot_clusts_onXaxis.png")
    , width = 10, height = 12, units = "in", res = 300)
p <- ggplot(plot_df, aes(x = Cluster , y = Term, color = -log10(PValue) , size = Fold_Enrichment)) +
  geom_point(alpha = 0.8) +
  scale_y_discrete(labels = identity) +
  scale_x_discrete(limits = cluster_order) +
  theme_minimal() +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = "Cluster",
    y = "",
    size = "Fold Enrichment",
    color = expression(-log[10](P-value)),
    title = "Cluster-level\nPathway Enrichment"
  ) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 8),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title  = element_text(hjust = 0.5)
  ) + scale_color_gradient(low = "blue", high = "red")
print(p)
dev.off()



