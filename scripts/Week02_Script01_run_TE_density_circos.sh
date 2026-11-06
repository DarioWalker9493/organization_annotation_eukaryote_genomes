#!/bin/bash
#SBATCH --job-name=TE_density
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
    "$WORKDIR/results/TE_annotation"

echo "========================================"
echo "Week 1 - Task 2.2: TE density"
echo "Genome: $GENOME"
echo "========================================"

# -------------------------------------------------
# Create FASTA index
# -------------------------------------------------

echo "Creating/checking FASTA index..."

apptainer exec -C \
    --bind "$USERROOT" \
    --bind "$COURSEDIR" \
    "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
    samtools faidx "$GENOME"

echo "FASTA index:"
head "$GENOME.fai"

# -------------------------------------------------
# Run circlize analysis
# -------------------------------------------------

echo "Running TE density analysis..."

apptainer exec -C \
    --bind "$USERROOT" \
    --bind "$COURSEDIR" \
    --env WORKDIR="$WORKDIR" \
    "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
    Rscript --vanilla \
    "$WORKDIR/scripts/Week02_Script01_TE_density_circos.R"

echo "========================================"
echo "Task 2.2 completed"
echo "========================================"