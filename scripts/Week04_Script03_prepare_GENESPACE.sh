#!/bin/bash
#SBATCH --job-name=GS_prepare
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dario.walker@students.unibe.ch

set -euo pipefail

WORKDIR=/data/users/dwalker/organization_annotation_eukaryote_genomes
COURSEDIR=/data/courses/genome_annotation_resources/2026

ACCESSION="Lu1"

FINALDIR="$WORKDIR/output/MAKER/final"
GENES="$FINALDIR/filtered.genes.renamed.gff3"
LONGEST_PROTEINS="$FINALDIR/assembly.all.maker.proteins.fasta.renamed.longest.fasta"

OUTDIR="$WORKDIR/output/comparisons"
GENESPACE_DIR="$OUTDIR/GENESPACE"

mkdir -p \
    "$GENESPACE_DIR/bed" \
    "$GENESPACE_DIR/peptide" \
    "$GENESPACE_DIR/tmp" \
    "$WORKDIR/logs"

echo "========================================"
echo "Week 4 - Prepare GENESPACE inputs"
echo "========================================"

# -------------------------------------------------
# Lu1 BED
# -------------------------------------------------

grep -P "\tgene\t" "$GENES" > "$GENESPACE_DIR/temp_genes.gff3"

awk 'BEGIN{OFS="\t"}
{
    split($9,a,";");
    split(a[1],b,"=");
    print $1, $4-1, $5, b[2]
}' "$GENESPACE_DIR/temp_genes.gff3" \
    > "$GENESPACE_DIR/bed/${ACCESSION}.bed"

# -------------------------------------------------
# Lu1 longest protein per gene
# Remove MAKER transcript suffix (-RA, -RB, ...)
# -------------------------------------------------

sed '/^>/s/ .*//' "$LONGEST_PROTEINS" \
    | sed '/^>/s/-R[A-Z]*$//' \
    > "$GENESPACE_DIR/peptide/${ACCESSION}.fa"

# -------------------------------------------------
# TAIR10
# -------------------------------------------------

cp "$COURSEDIR/data/TAIR10.bed" \
   "$GENESPACE_DIR/bed/TAIR10.bed"

cp "$COURSEDIR/data/TAIR10.fa" \
   "$GENESPACE_DIR/peptide/TAIR10.fa"

# -------------------------------------------------
# Three Lian et al. accessions
# -------------------------------------------------

for OTHER in Altai_5 Etna_2 Taz_0
do
    cp \
      "$COURSEDIR/data/Lian_et_al/GENESPACE/bed/${OTHER}.bed" \
      "$GENESPACE_DIR/bed/${OTHER}.bed"

    cp \
      "$COURSEDIR/data/Lian_et_al/GENESPACE/peptide/${OTHER}.fa" \
      "$GENESPACE_DIR/peptide/${OTHER}.fa"
done

# -------------------------------------------------
# Validate BED / FASTA IDs
# -------------------------------------------------

echo
echo "Validating GENESPACE inputs..."

for ACC in Lu1 TAIR10 Altai_5 Etna_2 Taz_0
do
    BED="$GENESPACE_DIR/bed/${ACC}.bed"
    FASTA="$GENESPACE_DIR/peptide/${ACC}.fa"

    cut -f4 "$BED" | sort > "$GENESPACE_DIR/${ACC}.bed.ids"
    grep '^>' "$FASTA" \
        | sed 's/^>//' \
        | awk '{print $1}' \
        | sort > "$GENESPACE_DIR/${ACC}.fa.ids"

    echo
    echo "=== $ACC ==="

    echo -n "BED genes: "
    wc -l < "$GENESPACE_DIR/${ACC}.bed.ids"

    echo -n "FASTA proteins: "
    wc -l < "$GENESPACE_DIR/${ACC}.fa.ids"

    echo -n "Duplicate BED IDs: "
    uniq -d "$GENESPACE_DIR/${ACC}.bed.ids" | wc -l

    echo -n "Duplicate FASTA IDs: "
    uniq -d "$GENESPACE_DIR/${ACC}.fa.ids" | wc -l

    echo -n "BED IDs missing from FASTA: "
    comm -23 \
      "$GENESPACE_DIR/${ACC}.bed.ids" \
      "$GENESPACE_DIR/${ACC}.fa.ids" \
      | wc -l

    echo -n "FASTA IDs missing from BED: "
    comm -13 \
      "$GENESPACE_DIR/${ACC}.bed.ids" \
      "$GENESPACE_DIR/${ACC}.fa.ids" \
      | wc -l
done

echo
echo "========================================"
echo "GENESPACE input preparation completed"
echo "========================================"
