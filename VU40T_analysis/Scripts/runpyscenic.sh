#!/bin/bash
#SBATCH --account=gendood-3dmucosa
#SBATCH --time=1-0:0:0
#SBATCH --qos=bbdefault
#SBATCH --mail-type=ALL
#SBATCH --nodes 1
#SBATCH --ntasks 48
#SBATCH --mem 120G

set -euo pipefail

## OPENBLAS and DASK settings
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export SLURM_CPU_BIND=cores

# Spill to local scratch and avoid nanny kills

# Node-local temp (fallback to /tmp if SLURM_TMPDIR is undefined)
export TMPDIR=${SLURM_TMPDIR:-/tmp}

export DASK_TEMPORARY_DIRECTORY="$TMPDIR"
export DASK_DISTRIBUTED__WORKER__TEMPORARY_DIRECTORY="$TMPDIR"
export DASK_DISTRIBUTED__WORKER__MEMORY__TARGET="0.80"
export DASK_DISTRIBUTED__WORKER__MEMORY__SPILL="0.90"
export DASK_DISTRIBUTED__WORKER__MEMORY__PAUSE="0.95"
export DASK_DISTRIBUTED__WORKER__MEMORY__TERMINATE="0.99"

# (nice-to-have on 2022.x)
export DASK_DISTRIBUTED__COMM__RETRY__COUNT=8
export DASK_DISTRIBUTED__COMM__TIMEOUTS__CONNECT="60s"
export DASK_DISTRIBUTED__COMM__TIMEOUTS__TCP="120s"

ulimit -n 65536  || true

### modules
module purge
module load bluebear
module load bear-apps/2022a
module load pySCENIC/0.12.1-foss-2022a

### files

cd /rds/projects/g/gendood-3dmucosa/scRNAseqAnalysis/human_epi_mtx

data_dir="../pyscenic_resources"
res_dir="../pyscenic_res"
mkdir -p "$data_dir" "$res_dir"



### download resources as recommended by tool authors
##hg38 TF list
# TF="$data_dir/allTFs_hg38.txt"
# if [[ ! -s "$TF" ]]; then
#   wget -q --show-progress -O "$TF" \
#     "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt"
# fi

## using curated TF list

TF="../OralMucosa/VU40T_analysis/Integrated/curated_TFs_hg38.txt"

##TF rankings and annotation as last index
# urls=(
#   "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
#   "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather"
#   "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
#   "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather"
#   "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
# )
## only need the ranking ones, scoring used for scenic+ which is ChIP-seq based
urls=(
  "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
  "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
  "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
)

for u in "${urls[@]}"; do
  fn="${u##*/}"                 # basename
  dest="$data_dir/$fn"
  echo "Checking $fn ..."
  if [[ ! -s "$dest" ]]; then   # check local file, not URL
    echo "  downloading → $dest"
    wget -q --show-progress -O "$dest" "$u" \
      || { echo "  FAILED: $u"; rm -f "$dest"; exit 1; }
  else
    echo "  already present (non-empty)"
  fi
done
# Download database directly (with wget or curl):

## for ranking and annotation files
fns=( "${urls[@]##*/}" )

### ranking dbs
rdbs=( "${fns[@]:0:${#fns[@]}-1}" )
### motif annotation
ann="${fns[@]: -1}"   # last element (note the space before -1)

## copy loom to scratch dir
cp -f pyscenic_compatible_human_epi_counts.loom "$TMPDIR/pyscenic_compatible_human_epi_counts.loom"

echo "performing XGBoost algo with grn"
pyscenic grn \
    --num_workers 20 \
    --seed 123 \
    -o "$TMPDIR/curatedTFs_pyscenic_comp_human_epi.adjacencies.tsv" \
    --cell_id_attribute "Cell" \
    --gene_attribute "Gene" \
    "$TMPDIR/pyscenic_compatible_human_epi_counts.loom" \
    "$TF"
## copy adjacencies out of scratch to PATH    
cp -f "$TMPDIR/curatedTFs_pyscenic_comp_human_epi.adjacencies.tsv" "$res_dir/curatedTFs_pyscenic_comp_human_epi.adjacencies.tsv"

echo "performing pyscenic ctx step"

pyscenic ctx \
    --num_workers 20 \
    $res_dir/curatedTFs_pyscenic_comp_human_epi.adjacencies.tsv \
    "${rdbs[@]/#/$data_dir/}" \
    --annotations_fname "$data_dir/$ann" \
    --cell_id_attribute "Cell" \
    --gene_attribute "Gene" \
    --expression_mtx_fname pyscenic_compatible_human_epi_counts.loom \
    --mode "dask_multiprocessing" \
    --output $res_dir/curatedTFs_pyscenic_fixed_human_epi_regulons.csv \
    --mask_dropouts \
    --all_modules     

echo "converting regulons to aucell distances"

pyscenic aucell \
        pyscenic_compatible_human_epi_counts.loom \
        $res_dir/curatedTFs_pyscenic_fixed_human_epi_regulons.csv  \
        -o $res_dir/curatedTFs_pyscenic_fixed_human_epi_auc_mtx.csv \
        --seed 123 \
        --cell_id_attribute "Cell" \
        --gene_attribute "Gene" \
        --num_workers 20
