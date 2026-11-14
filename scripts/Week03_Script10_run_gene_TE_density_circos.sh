#!/bin/bash
#SBATCH --job-name=Gene_TE_density
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
USERROOT=/data/users/dwalker

GENOME="$WORKDIR/data/assemblies/Lu1_Flye.fasta"

mkdir -p \
    "$WORKDIR/logs" \
    "$WORKDIR/plots" \
    "$WORKDIR/results/gene_annotation"

echo "========================================"
echo "Week 3 - Gene vs TE density"
echo "========================================"

# Ensure assembly index exists
apptainer exec -C \
    --bind "$USERROOT" \
    --bind "$COURSEDIR" \
    "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
    samtools faidx "$GENOME"

echo "Running gene/TE density analysis..."

apptainer exec -C \
    --bind "$USERROOT" \
    --bind "$COURSEDIR" \
    --env WORKDIR="$WORKDIR" \
    "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
    Rscript --vanilla \
    "$WORKDIR/scripts/Week03_Script10_gene_TE_density_circos.R"

echo
echo "Generated outputs:"

ls -lh \
    "$WORKDIR/plots/06_Gene_TE_density_circos.png" \
    "$WORKDIR/plots/06_Gene_TE_density_circos.pdf" \
    "$WORKDIR/results/gene_annotation/Week03_gene_TE_density_top10_scaffolds.tsv" \
    "$WORKDIR/results/gene_annotation/Week03_gene_TE_density_100kb_windows.tsv" \
    "$WORKDIR/results/gene_annotation/Week03_gene_TE_density_correlations.tsv"

echo
echo "========================================"
echo "Gene vs TE density completed"
echo "========================================"
