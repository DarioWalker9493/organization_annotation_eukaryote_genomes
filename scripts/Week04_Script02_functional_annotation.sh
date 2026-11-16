#!/bin/bash
#SBATCH --job-name=W4_function
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dario.walker@students.unibe.ch

set -euo pipefail

WORKDIR=/data/users/dwalker/organization_annotation_eukaryote_genomes
COURSEDIR=/data/courses/genome_annotation_resources/2026

FINALDIR="$WORKDIR/output/MAKER/final"
OUTDIR="$WORKDIR/output/comparisons"

PROTEINS="$FINALDIR/assembly.all.maker.proteins.fasta.renamed.filtered.fasta"
GENES="$FINALDIR/filtered.genes.renamed.gff3"
UNIPROT="$COURSEDIR/data/uniprot_viridiplantae_reviewed.fa"

BESTHITS="$OUTDIR/blastp_uniprot.besthits.tsv"
MAKER_SIF="$COURSEDIR/containers/MAKER_3.01.03.sif"

mkdir -p "$OUTDIR/tmp" "$WORKDIR/logs"

for FILE in "$PROTEINS" "$GENES" "$UNIPROT" "$BESTHITS"
do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Missing or empty file: $FILE"
        exit 1
    fi
done

cd "$OUTDIR"

echo "========================================"
echo "Week 4 - UniProt functional annotation"
echo "========================================"

apptainer exec \
    --cleanenv \
    --env TMPDIR=/tmp \
    --bind /data \
    --bind "$WORKDIR" \
    --bind "$OUTDIR/tmp":/tmp \
    --pwd "$OUTDIR" \
    "$MAKER_SIF" \
    maker_functional_fasta \
    "$UNIPROT" \
    "$BESTHITS" \
    "$PROTEINS" \
    > maker_proteins.filtered.Uniprot.fasta

apptainer exec \
    --cleanenv \
    --env TMPDIR=/tmp \
    --bind /data \
    --bind "$WORKDIR" \
    --bind "$OUTDIR/tmp":/tmp \
    --pwd "$OUTDIR" \
    "$MAKER_SIF" \
    maker_functional_gff \
    "$UNIPROT" \
    "$BESTHITS" \
    "$GENES" \
    > filtered.genes.Uniprot.gff3

echo
echo "Outputs:"
ls -lh \
    maker_proteins.filtered.Uniprot.fasta \
    filtered.genes.Uniprot.gff3

echo
echo "========================================"
echo "Functional annotation completed"
echo "========================================"
