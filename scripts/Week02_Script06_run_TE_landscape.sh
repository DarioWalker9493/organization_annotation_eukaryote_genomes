#!/bin/bash
#SBATCH --job-name=TE_landscape
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dario.walker@students.unibe.ch

set -euo pipefail

COURSEDIR=/data/courses/genome_annotation_resources/2026
WORKDIR=/data/users/dwalker/organization_annotation_eukaryote_genomes

mkdir -p \
    "$WORKDIR/plots" \
    "$WORKDIR/results/TE_annotation"

echo "========================================"
echo "Week 2 - Task 4"
echo "TE divergence landscape and age"
echo "========================================"

apptainer exec -C \
    --bind "$WORKDIR" \
    --bind "$COURSEDIR" \
    --env WORKDIR="$WORKDIR" \
    "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
    Rscript --vanilla \
    "$WORKDIR/scripts/Week02_Script06_TE_landscape.R"

echo
echo "Generated plots:"
ls -lh \
    "$WORKDIR/plots/05_TE_divergence_landscape.png" \
    "$WORKDIR/plots/05_TE_divergence_landscape.pdf"

echo
echo "Generated summaries:"
ls -lh \
    "$WORKDIR/results/TE_annotation/Week02_TE_divergence_landscape.tsv" \
    "$WORKDIR/results/TE_annotation/Week02_TE_landscape_peak_summary.tsv"

echo
echo "========================================"
echo "Task 4 completed"
echo "========================================"