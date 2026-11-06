#!/bin/bash
#SBATCH --job-name=LTR_identity
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dario.walker@students.unibe.ch

set -euo pipefail

COURSEDIR=/data/courses/genome_annotation_resources/2026
WORKDIR=/data/users/dwalker/organization_annotation_eukaryote_genomes
OUTDIR="$WORKDIR/output/LTR_full_length"
PLOT_DIR="$WORKDIR/plots"

mkdir -p "$PLOT_DIR"

cd "$OUTDIR"

apptainer exec -C \
    --bind "$WORKDIR" \
    --bind "$COURSEDIR" \
    --env PLOT_DIR="$PLOT_DIR" \
    --pwd "$OUTDIR" \
    "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
    Rscript --vanilla \
    "$WORKDIR/scripts/Week01_Script03_full_length_LTR_identity.R"