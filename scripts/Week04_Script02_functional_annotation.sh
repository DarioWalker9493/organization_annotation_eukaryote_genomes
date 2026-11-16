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

ANNOTATOR="$WORKDIR/scripts/Week04_Script02_annotate_uniprot_gff.py"

OUTPUT_FASTA="$OUTDIR/maker_proteins.filtered.Uniprot.fasta"
OUTPUT_GFF="$OUTDIR/filtered.genes.Uniprot.gff3"
OUTPUT_TABLE="$OUTDIR/Week04_uniprot_function_summary.tsv"

TMPDIR="$OUTDIR/tmp"

mkdir -p \
    "$TMPDIR" \
    "$WORKDIR/logs"

cd "$OUTDIR"


echo "========================================"
echo "Week 4 - Functional annotation"
echo "========================================"


# ============================================================================
# Preflight
# ============================================================================

for FILE in \
    "$PROTEINS" \
    "$GENES" \
    "$UNIPROT" \
    "$BESTHITS" \
    "$ANNOTATOR"
do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Missing or empty file:"
        echo "$FILE"
        exit 1
    fi
done


echo
echo "Inputs:"

ls -lh \
    "$PROTEINS" \
    "$GENES" \
    "$UNIPROT" \
    "$BESTHITS"


# ============================================================================
# 1. Add UniProt functions to filtered protein FASTA
#
# maker_functional_fasta worked correctly in the previous run and is retained.
# ============================================================================

echo
echo "Adding UniProt functions to protein FASTA..."

apptainer exec \
    --cleanenv \
    --env TMPDIR=/tmp \
    --bind /data \
    --bind "$WORKDIR" \
    --bind "$TMPDIR":/tmp \
    --pwd "$OUTDIR" \
    "$MAKER_SIF" \
    maker_functional_fasta \
    "$UNIPROT" \
    "$BESTHITS" \
    "$PROTEINS" \
    > "$OUTPUT_FASTA"


# ============================================================================
# 2. Add UniProt functions to final GFF3
#
# maker_functional_gff fails on this MAKER version with a Perl HASH-reference
# error despite the GFF being structurally valid. We therefore perform the
# same mapping explicitly from the UniProt FASTA and BLAST best-hit table.
# ============================================================================

echo
echo "Adding UniProt functions to final GFF3..."

python3 "$ANNOTATOR" \
    --uniprot "$UNIPROT" \
    --besthits "$BESTHITS" \
    --proteins "$PROTEINS" \
    --gff "$GENES" \
    --output-gff "$OUTPUT_GFF" \
    --output-table "$OUTPUT_TABLE"


# ============================================================================
# 3. Final validation
# ============================================================================

echo
echo "========================================"
echo "Validation"
echo "========================================"

INPUT_PROTEINS=$(grep -c '^>' "$PROTEINS")
OUTPUT_PROTEINS=$(grep -c '^>' "$OUTPUT_FASTA")

INPUT_GENES=$(awk -F '\t' '$3=="gene"{n++} END{print n+0}' "$GENES")
OUTPUT_GENES=$(awk -F '\t' '$3=="gene"{n++} END{print n+0}' "$OUTPUT_GFF")

INPUT_MRNAS=$(awk -F '\t' '$3=="mRNA"{n++} END{print n+0}' "$GENES")
OUTPUT_MRNAS=$(awk -F '\t' '$3=="mRNA"{n++} END{print n+0}' "$OUTPUT_GFF")


echo "Input proteins:    $INPUT_PROTEINS"
echo "Output proteins:   $OUTPUT_PROTEINS"

echo "Input genes:       $INPUT_GENES"
echo "Output genes:      $OUTPUT_GENES"

echo "Input mRNAs:       $INPUT_MRNAS"
echo "Output mRNAs:      $OUTPUT_MRNAS"


if [[ "$INPUT_PROTEINS" -ne "$OUTPUT_PROTEINS" ]]; then
    echo "ERROR: Protein count changed."
    exit 1
fi

if [[ "$INPUT_GENES" -ne "$OUTPUT_GENES" ]]; then
    echo "ERROR: Gene count changed."
    exit 1
fi

if [[ "$INPUT_MRNAS" -ne "$OUTPUT_MRNAS" ]]; then
    echo "ERROR: mRNA count changed."
    exit 1
fi


echo
echo "Features carrying UniProt cross-references:"

grep -c 'UniProtKB:' "$OUTPUT_GFF" || true


echo
echo "Features carrying UniProt functional Notes:"

grep -c 'Note=' "$OUTPUT_GFF" || true


echo
echo "Outputs:"

ls -lh \
    "$OUTPUT_FASTA" \
    "$OUTPUT_GFF" \
    "$OUTPUT_TABLE"


echo
echo "========================================"
echo "Week 4 functional annotation completed"
echo "========================================"

