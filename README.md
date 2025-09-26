# OralMucosa
## VU40T Analysis
An R-based pipeline for analysing scRNAseq data from VU40T human sublingual epithelial cells and 3T3-murine fibroblast cells derived from 3D-raft cultures
---
## 📁 Project Structure
```bash
README.md
VU40T_analysis/
├── Scripts/               # R and shell scripts used in the pipeline
│   ├── createscrnaseqsamplesheet.sh     # automate samplesheet creation for nf-core/scrnaseq input
│   ├── runnf-corescrnaseq.slurm    # slurm shell script to run nf-core/scrnaseq (v. 4.0.0) on raw FastQ data
│   ├── species_mixing.R    # loads in both human and mouse alignments and assigns species labels to cell barcodes
│   ├── pre-integration_seurat.R    # QC for both fibroblasts and epithelial cells
│   ├── DoubletFinder.R     # Pre-processing cont. and single sample analysis
│   ├── Integration_seurat.R    # Integration of preprocessed data from each species seperately, outputs umaps, pseudobulk volcanos and GSEAs and marker dotplots.
│   ├── GSEA_plotting.R     # Script to plot GSEA results from DAVID as dotplots 
│   ├── epithelial_pseudotime.R     # pseudotime analysis for epithelial cells with monocle3
│   ├── fibroblast_pseudotime.R     # pseudotime analysis for fibroblast cells with monocle3
│   ├── runpyscenic.sh   # shellscript to run pyscenic on the curated TF list with gene names mapped onto older aliases in ranking DB
│   └── cellchat.R     # Ligand-receptor cellular communication analysis between epithelial and fibroblast clusters with cellchat
├── human_multiqc_report.html     # from nf-core/scrnaseq pipeline run
├── Mixed_species_plots/     # Output species assignment figures 
├── VU40T_species_assignment_by_reads.csv   # species assignment table
├── Seperate_samples/   # figs and tables from pre-integration_seurat.R -> DoubletFinder.R for both human and mouse alignments.
│   ├── Doublet_detection/  #  doublets summary from doubletfinder on human alignment
│   ├── Plots/  #  Plots from seurat analysis on seperate samples using human alignment only - pre-integration
│   └── Mouse/      # plots and markers from seurat analysis on seperate samples using mouse alignment only - pre-integration 
│   │   ├── Plots/  # final pre integration plots on mouse alignment only
│   │   └── QC_plots/   # diagnostic plots from     pre-integration_seurat.R - mouse alignment only
└── Integrated/   # figs and tables from Integration_seurat.R for human alignments in this dir and Mouse alignments in a subdir.
│   ├── Plots/  #  Plots from seurat analysis on Integrated dataset using human alignment only. includes pseudobulk and bulk plots.
│   ├── Pseudotime_humanOnly/  #  Epithelial pseudotime DEG list which passed monocle3 QC
│   └── Mouse/      #   Integrated Fibroblast (mouse alignment only) tables amd plots
│   │   └── Plots/    #  Plots from seurat analysis on Integrated dataset using mouse alignment only. includes pseudobulk, bulk and pseudotime plots. pseudotime DEGs are in the Pseudotime_MouseOnly/ subdir
```


## 🛠️ Dependencies
This analysis was run using R (v. 4.4.1) and the following R packages:

```bash
sessionInfo()

 [1] clusterProfiler_4.12.0     
 [2] enrichplot_1.24.0          
 [3] msigdbr_24.1.0             
 [4] pagoda2_1.0.12             
 [5] igraph_2.1.4               
 [6] Matrix_1.7-3               
 [7] clustree_0.5.1             
 [8] ggraph_2.2.1               
 [9] patchwork_1.3.1            
[10] DoubletFinder_2.0.6        
[11] parallelly_1.45.0          
[12] biomaRt_2.60.0             
[13] ggridges_0.5.6             
[14] ggplot2_3.5.2              
[15] SeuratWrappers_0.4.0       
[16] monocle3_1.4.26            
[17] SingleCellExperiment_1.26.0
[18] SummarizedExperiment_1.34.0
[19] GenomicRanges_1.56.2       
[20] GenomeInfoDb_1.42.3        
[21] IRanges_2.40.1             
[22] S4Vectors_0.44.0           
[23] MatrixGenerics_1.16.0      
[24] matrixStats_1.5.0          
[25] Biobase_2.66.0             
[26] BiocGenerics_0.52.0        
[27] dplyr_1.1.4                
[28] Seurat_5.3.0               
[29] SeuratObject_5.1.0         
[30] sp_2.2-0  
[31] EnhancedVolcano_1.24.0     
[32] ggrepel_0.9.6
[33] CellChat_2.1.2
```

and Python (v. 4.4.1) libraries/tools:

```bash
pySCENIC v 0.12.1
```
