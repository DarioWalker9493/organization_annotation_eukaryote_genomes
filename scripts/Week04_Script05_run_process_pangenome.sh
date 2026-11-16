#!/bin/bash
#SBATCH --job-name=GS_pangenome
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dario.walker@students.unibe.ch

set -euo pipefail

WORKDIR=/data/users/dwalker/organization_annotation_eukaryote_genomes
COURSEDIR=/data/courses/genome_annotation_resources/2026

ANALYSIS_SIF="$COURSEDIR/containers/Annotation_analysis_R4.4.sif"
TMPDIR="$WORKDIR/output/comparisons/Annotation_R_tmp"

mkdir -p "$TMPDIR" "$WORKDIR/logs"

apptainer exec \
    --cleanenv \
    --bind /data \
    --bind "$WORKDIR" \
    --bind "$TMPDIR":/temp \
    --env TMPDIR=/temp \
    --pwd "$WORKDIR" \
    "$ANALYSIS_SIF" \
    Rscript \
    "$WORKDIR/scripts/Week04_Script05_process_pangenome.R"
