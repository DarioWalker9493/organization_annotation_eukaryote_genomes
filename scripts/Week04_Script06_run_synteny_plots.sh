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

ls -lh \
    "$GENESPACE_DIR/results/gsParams.rda"


echo
echo "Container:"

ls -lh \
    "$GENESPACE_SIF"


# ============================================================================
# Run R visualization
# ============================================================================

echo
echo "Running GENESPACE riparian plots..."


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
# Convert final PDFs to PNG
# ============================================================================

echo
echo "Converting PDFs to PNG..."


for NUMBER in 09 10
do

    PDF=$(find "$PLOTDIR" \
        -maxdepth 1 \
        -type f \
        -name "${NUMBER}_*.pdf" \
        -print \
        -quit)

    if [[ -z "$PDF" ]]; then

        echo "ERROR: Plot ${NUMBER} PDF not found."
        exit 1

    fi

    PREFIX="${PDF%.pdf}"

    pdftoppm \
        -png \
        -singlefile \
        -r 300 \
        "$PDF" \
        "$PREFIX"

done


# ============================================================================
# Validation
# ============================================================================

echo
echo "========================================"
echo "Final Week 5 plots"
echo "========================================"


ls -lh \
    "$PLOTDIR"/09_* \
    "$PLOTDIR"/10_*


echo
echo "========================================"
echo "GENESPACE synteny visualization complete"
echo "========================================"

