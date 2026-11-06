#!/bin/bash
#SBATCH --job-name=TE_clades
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
    "$WORKDIR/results/TE_annotation" \
    "$WORKDIR/plots"

echo "========================================"
echo "Week 2 - Task 3"
echo "Copia/Gypsy clade genomic abundance"
echo "========================================"

apptainer exec -C \
    --bind "$WORKDIR" \
    --bind "$COURSEDIR" \
    --env WORKDIR="$WORKDIR" \
    "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
    Rscript --vanilla \
    "$WORKDIR/scripts/Week02_Script03_clade_abundance.R"

echo
echo "========================================"
echo "Task 3 clade abundance analysis completed"
echo "========================================"