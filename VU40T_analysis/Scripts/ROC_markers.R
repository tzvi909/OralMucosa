###libs
library(Seurat)

## to avoid hitting that 20gb quota on home dir,
proj_dir <- "/rds/projects/g/gendood-3dmucosa/"
analysis_dir <- file.path(proj_dir, "scRNAseqAnalysis/")
git_dir <- file.path(analysis_dir, "OralMucosa/VU40T_analysis")
out_prefix_dir <- file.path(git_dir, "Integrated")

cache_dir <- file.path(proj_dir, "rds_cache")

## check for dirs recursively

chk_dir_list <- list(analysis_dir, git_dir, cache_dir, out_prefix_dir)

for (path in chk_dir_list){
  if(!(dir.exists(path))){
    dir.create(path, recursive = T)
  }
}

species <- c("human", "mouse")

for (i in species){

    res_dir <- file.path(out_prefix_dir, paste0(stringr::str_to_title(i)))

    
    if (i == "human"){
        seu <- readRDS(file.path(cache_dir, "VU40T_combined_joined_0.3_res_humanOnly.rds"))
    }else{
        seu <- readRDS(file.path(cache_dir, "VU40T_combined_joined_0.6_res_mouseOnly.rds")) 
    }
    
    DefaultAssay(seu) <- "RNA"
    
    markers_lfc025 <- FindAllMarkers(seu,
                              only.pos = TRUE,
                              min.pct = 0,
                              logfc.threshold = 0.25,
                              test.use = "roc",
                              densify = T)
    write.csv(markers_lfc025, file.path(res_dir, paste0(i,"_lfc0.25_Markers_roc.csv")))
    markers_all <-  FindAllMarkers(seu,
                              only.pos = TRUE,
                              min.pct = 0,
                              logfc.threshold = 0,
                              test.use = "roc",
                              densify = T)      
    write.csv(markers_all, file.path(res_dir, paste0(i,"_all_Markers_roc.csv")))
}

res_dir <- file.path(out_prefix_dir, paste0(stringr::str_to_title(i)))

# Load canonical human cell-cycle gene lists
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

cc <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes)
VlnPlot(cc, c("S.Score","G2M.Score"))


markers_lfc025 <- read.csv(file.path(res_dir, paste0("human","_lfc0.25_Markers_roc.csv")))
