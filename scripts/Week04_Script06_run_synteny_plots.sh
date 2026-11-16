#!/bin/bash
#SBATCH --job-name=GS_synteny
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dario.walker@students.unibe.ch

set -euo pipefail


# ============================================================================
# Settings
# ============================================================================

WORKDIR=/data/users/dwalker/organization_annotation_eukaryote_genomes
COURSEDIR=/data/courses/genome_annotation_resources/2026

GENESPACE_DIR="$WORKDIR/output/comparisons/GENESPACE"
GENESPACE_SIF="$COURSEDIR/containers/genespace_latest.sif"

RSCRIPT="$WORKDIR/scripts/Week04_Script06_synteny_plots.R"

TMPDIR="$WORKDIR/output/comparisons/GENESPACE_plot_tmp"
PLOTDIR="$WORKDIR/plots"

mkdir -p \
    "$TMPDIR" \
    "$PLOTDIR" \
    "$WORKDIR/logs"


# ============================================================================
# Preflight
# ============================================================================

echo "========================================"
echo "Week 5 - GENESPACE synteny visualization"
echo "========================================"


for FILE in \
    "$RSCRIPT" \
    "$GENESPACE_DIR/results/gsParams.rda" \
    "$GENESPACE_SIF"
do

    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Missing or empty required file:"
        echo "$FILE"
        exit 1
    fi

done


echo
echo "GENESPACE parameter file:"
ls -lh "$GENESPACE_DIR/results/gsParams.rda"

echo
echo "Container:"
ls -lh "$GENESPACE_SIF"


# ============================================================================
# Run riparian plots
# ============================================================================

echo
echo "Creating riparian plots..."


apptainer exec \
    --cleanenv \
    --bind /data \
    --bind "$WORKDIR" \
    --bind "$TMPDIR":/temp \
    --env TMPDIR=/temp \
    --pwd "$WORKDIR" \
    "$GENESPACE_SIF" \
    Rscript \
    "$RSCRIPT"


# ============================================================================
# Copy GENESPACE-generated Lu1 vs TAIR10 dotplots
# ============================================================================

echo
echo "Copying Lu1 vs TAIR10 dotplots..."


RAW_DOTPLOT="$GENESPACE_DIR/dotplots/Lu1_vs_TAIR10.rawHits.pdf"

SYNTENIC_DOTPLOT="$GENESPACE_DIR/dotplots/Lu1_vs_TAIR10.syntenicHits.pdf"


for FILE in \
    "$RAW_DOTPLOT" \
    "$SYNTENIC_DOTPLOT"
do

    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Missing expected GENESPACE dotplot:"
        echo "$FILE"
        exit 1
    fi

done


cp \
    "$RAW_DOTPLOT" \
    "$PLOTDIR/11_Lu1_vs_TAIR10_rawHits.pdf"


cp \
    "$SYNTENIC_DOTPLOT" \
    "$PLOTDIR/12_Lu1_vs_TAIR10_syntenicHits.pdf"


# ============================================================================
# Validation
# ============================================================================

echo
echo "========================================"
echo "Final Week 5 plots"
echo "========================================"


for FILE in \
    "$PLOTDIR/09_Lu1_vs_TAIR10_riparian.pdf" \
    "$PLOTDIR/10_all_genomes_TAIR10_riparian.pdf" \
    "$PLOTDIR/11_Lu1_vs_TAIR10_rawHits.pdf" \
    "$PLOTDIR/12_Lu1_vs_TAIR10_syntenicHits.pdf"
do

    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Missing expected final plot:"
        echo "$FILE"
        exit 1
    fi

    ls -lh "$FILE"

done


echo
echo "========================================"
echo "GENESPACE synteny visualization complete"
echo "========================================"

