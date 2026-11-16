#!/bin/bash
#SBATCH --job-name=GENESPACE
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dario.walker@students.unibe.ch

set -euo pipefail

WORKDIR=/data/users/dwalker/organization_annotation_eukaryote_genomes
COURSEDIR=/data/courses/genome_annotation_resources/2026

GENESPACE_DIR="$WORKDIR/output/comparisons/GENESPACE"
GENESPACE_SIF="$COURSEDIR/containers/genespace_latest.sif"

TMPDIR="$WORKDIR/output/comparisons/GENESPACE_R_tmp"

mkdir -p "$TMPDIR" "$WORKDIR/logs"

cd "$WORKDIR"

apptainer exec \
    --cleanenv \
    --bind /data \
    --bind "$WORKDIR" \
    --bind "$TMPDIR":/temp \
    --env TMPDIR=/temp \
    --pwd "$WORKDIR" \
    "$GENESPACE_SIF" \
    Rscript \
    "$WORKDIR/scripts/Week04_Script04_GENESPACE.R" \
    "$GENESPACE_DIR"
