#!/bin/bash
#SBATCH --job-name=InterProScan
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dario.walker@students.unibe.ch

set -euo pipefail

COURSEDIR=/data/courses/genome_annotation_resources/2026
WORKDIR=/data/users/dwalker/organization_annotation_eukaryote_genomes
OUTDIR="$WORKDIR/output/MAKER"
FINALDIR="$OUTDIR/final"

protein="assembly.all.maker.proteins.fasta"
INTERPROSCAN_SIF="$COURSEDIR/containers/InterProScan_5.70.sif"

mkdir -p \
    "$FINALDIR/tmp" \
    "$WORKDIR/logs"

cd "$FINALDIR"

echo "========================================"
echo "Week 3 - InterProScan / Pfam"
echo "========================================"

INPUT="${protein}.renamed.fasta"

if [[ ! -s "$INPUT" ]]; then
    echo "ERROR: Protein FASTA not found or empty:"
    echo "$INPUT"
    exit 1
fi

echo "Input:"
ls -lh "$INPUT"

echo
echo -n "Protein sequences: "
grep -c '^>' "$INPUT"

echo
echo "Starting InterProScan..."
echo

apptainer exec \
    --cleanenv \
    --pwd "$FINALDIR" \
    --bind "$COURSEDIR/data/interproscan-5.70-102.0/data":/opt/interproscan/data \
    --bind "$COURSEDIR" \
    --bind "$WORKDIR" \
    --bind "$FINALDIR/tmp":/temp \
    "$INTERPROSCAN_SIF" \
    /opt/interproscan/interproscan.sh \
    -appl pfam \
    --disable-precalc \
    -f TSV \
    --goterms \
    --iprlookup \
    --seqtype p \
    -i "$INPUT" \
    -o output.iprscan

echo
echo "========================================"
echo "InterProScan completed"
echo "========================================"

ls -lh output.iprscan

echo
echo -n "InterProScan result lines: "
wc -l < output.iprscan

echo -n "Proteins with at least one Pfam hit: "
cut -f1 output.iprscan | sort -u | wc -l
