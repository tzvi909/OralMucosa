#!/bin/bash
#SBATCH --account=gendood-3dmucosa
#SBATCH --time=1-0:0:0
#SBATCH --qos=bbdefault
#SBATCH --mail-type=ALL
#SBATCH --nodes 1
#SBATCH --ntasks 16

ulimit -n 10000

## dependencies
module purge; module load bluebear
module load bear-apps/2022b
module load Nextflow/24.04.2

export NXF_SINGULARITY_CACHEDIR="/rds/projects/g/gendood-3dmucosa/.apptainer"

#main
nextflow run nf-core/scrnaseq -r 4.0.0 \
    --input ../VU40Tsamplesheet.csv \
    --outdir ../LPS_VU40T_QC_and_counts \
    --genome GRCh38 \
    -profile bluebear \
    --aligner cellranger \
    --fasta ../../Refs/hg38.fa \
    --gtf ../../Refs/hg38.ncbiRefSeq.gtf