#!/bin/bash
#SBATCH --job-name=Maker_merge
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dario.walker@students.unibe.ch

set -euo pipefail

COURSEDIR=/data/courses/genome_annotation_resources/2026
WORKDIR=/data/users/dwalker/organization_annotation_eukaryote_genomes
OUTDIR="$WORKDIR/output/MAKER"

MAKER_SIF="$COURSEDIR/containers/MAKER_3.01.03.sif"

mkdir -p \
    "$OUTDIR/tmp" \
    "$WORKDIR/logs"

cd "$OUTDIR"

echo "========================================"
echo "Week 3 - Merge MAKER outputs"
echo "========================================"

shopt -s nullglob

INDEXES=(
    *.maker.output/*_master_datastore_index.log
)

if [[ "${#INDEXES[@]}" -ne 1 ]]; then
    echo "ERROR: Expected exactly one MAKER datastore index."
    echo "Found: ${#INDEXES[@]}"
    printf '%s\n' "${INDEXES[@]}"
    exit 1
fi

INDEX="${INDEXES[0]}"

echo
echo "Using datastore index:"
echo "$INDEX"

# -------------------------------------------------
# Merge complete GFF
# -------------------------------------------------

echo
echo "Merging full MAKER GFF..."

apptainer exec \
    --cleanenv \
    --env TMPDIR=/tmp \
    --bind "$OUTDIR/tmp":/tmp \
    --bind "$COURSEDIR" \
    --bind /data/users/dwalker \
    --pwd "$OUTDIR" \
    "$MAKER_SIF" \
    gff3_merge \
    -d "$INDEX" \
    -o assembly.all.maker.gff

# -------------------------------------------------
# Merge GFF without embedded sequence
# -------------------------------------------------

echo
echo "Merging MAKER GFF without sequence..."

apptainer exec \
    --cleanenv \
    --env TMPDIR=/tmp \
    --bind "$OUTDIR/tmp":/tmp \
    --bind "$COURSEDIR" \
    --bind /data/users/dwalker \
    --pwd "$OUTDIR" \
    "$MAKER_SIF" \
    gff3_merge \
    -n \
    -d "$INDEX" \
    -o assembly.all.maker.noseq.gff

# -------------------------------------------------
# Merge protein and transcript FASTA
# -------------------------------------------------

echo
echo "Merging protein and transcript sequences..."

apptainer exec \
    --cleanenv \
    --env TMPDIR=/tmp \
    --bind "$OUTDIR/tmp":/tmp \
    --bind "$COURSEDIR" \
    --bind /data/users/dwalker \
    --pwd "$OUTDIR" \
    "$MAKER_SIF" \
    fasta_merge \
    -d "$INDEX" \
    -o assembly

# -------------------------------------------------
# Basic checks
# -------------------------------------------------

echo
echo "Generated files:"

ls -lh \
    assembly.all.maker.gff \
    assembly.all.maker.noseq.gff \
    assembly.all.maker.proteins.fasta \
    assembly.all.maker.transcripts.fasta

echo
echo "Gene models:"
awk -F '\t' '$3 == "gene" {n++} END {print n+0}' \
    assembly.all.maker.noseq.gff

echo
echo "mRNA models:"
awk -F '\t' '$3 == "mRNA" {n++} END {print n+0}' \
    assembly.all.maker.noseq.gff

echo
echo "Proteins:"
grep -c '^>' assembly.all.maker.proteins.fasta

echo
echo "Transcripts:"
grep -c '^>' assembly.all.maker.transcripts.fasta

echo
echo "========================================"
echo "MAKER merge completed"
echo "========================================"