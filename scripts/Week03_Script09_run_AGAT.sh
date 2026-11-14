#!/bin/bash
#SBATCH --job-name=AGAT_MAKER
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dario.walker@students.unibe.ch

set -euo pipefail

COURSEDIR=/data/courses/genome_annotation_resources/2026
WORKDIR=/data/users/dwalker/organization_annotation_eukaryote_genomes
FINALDIR="$WORKDIR/output/MAKER/final"

AGAT_SIF="$COURSEDIR/containers/AGAT_v1.4.1.sif"

mkdir -p "$WORKDIR/logs"

cd "$FINALDIR"

INPUT="filtered.genes.renamed.gff3"
OUTPUT="annotation.stat"

if [[ ! -s "$INPUT" ]]; then
    echo "ERROR: Missing or empty input: $INPUT"
    exit 1
fi

echo "========================================"
echo "Week 3 - AGAT annotation statistics"
echo "========================================"

echo
echo "Input:"
ls -lh "$INPUT"

echo
echo "Running AGAT..."

apptainer exec \
    --cleanenv \
    --bind "$COURSEDIR" \
    --bind "$WORKDIR" \
    --pwd "$FINALDIR" \
    "$AGAT_SIF" \
    agat_sp_statistics.pl \
    --gff "$INPUT" \
    --output "$OUTPUT"

echo
echo "========================================"
echo "AGAT completed"
echo "========================================"

ls -lh "$OUTPUT"

echo
echo "Annotation statistics:"
cat "$OUTPUT"
