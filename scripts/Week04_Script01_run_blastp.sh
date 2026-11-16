#!/bin/bash
#SBATCH --job-name=BLASTP_W4
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dario.walker@students.unibe.ch

set -euo pipefail

WORKDIR=/data/users/dwalker/organization_annotation_eukaryote_genomes
COURSEDIR=/data/courses/genome_annotation_resources/2026

FINALDIR="$WORKDIR/output/MAKER/final"
PROTEINS="$FINALDIR/assembly.all.maker.proteins.fasta.renamed.filtered.fasta"

UNIPROT="$COURSEDIR/data/uniprot_viridiplantae_reviewed.fa"
TAIR10="$COURSEDIR/data/TAIR10_pep_20110103_representative_gene_model"

OUTDIR="$WORKDIR/output/comparisons"

mkdir -p \
    "$OUTDIR/tmp" \
    "$WORKDIR/logs"

cd "$OUTDIR"

echo "========================================"
echo "Week 4 - BLASTP functional comparison"
echo "========================================"

# -------------------------------------------------
# Preflight
# -------------------------------------------------

for FILE in "$PROTEINS" "$UNIPROT" "$TAIR10"
do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Missing or empty input:"
        echo "$FILE"
        exit 1
    fi
done

echo
echo "Input proteins:"
ls -lh "$PROTEINS"

echo -n "Protein sequences: "
grep -c '^>' "$PROTEINS"

echo
echo "Reference databases:"
ls -lh "$UNIPROT" "$TAIR10"

# -------------------------------------------------
# BLAST+
# -------------------------------------------------

module load BLAST+/2.15.0-gompi-2021a

echo
echo "BLAST version:"
blastp -version | head -n 1

# -------------------------------------------------
# UniProt and TAIR10 searches in parallel
# -------------------------------------------------

echo
echo "Starting UniProt BLASTP..."

blastp \
    -query "$PROTEINS" \
    -db "$UNIPROT" \
    -num_threads 10 \
    -outfmt 6 \
    -evalue 1e-5 \
    -max_target_seqs 10 \
    -out "$OUTDIR/blastp_uniprot.tsv" &

uniprot_job=$!

echo "Starting TAIR10 BLASTP..."

blastp \
    -query "$PROTEINS" \
    -db "$TAIR10" \
    -num_threads 10 \
    -outfmt 6 \
    -evalue 1e-5 \
    -max_target_seqs 10 \
    -out "$OUTDIR/blastp_TAIR10.tsv" &

tair_job=$!

# -------------------------------------------------
# Wait for both searches
# -------------------------------------------------

wait "$uniprot_job"
uniprot_status=$?

wait "$tair_job"
tair_status=$?

if [[ "$uniprot_status" -ne 0 ]]; then
    echo "ERROR: UniProt BLASTP failed."
    exit "$uniprot_status"
fi

if [[ "$tair_status" -ne 0 ]]; then
    echo "ERROR: TAIR10 BLASTP failed."
    exit "$tair_status"
fi

# -------------------------------------------------
# Select highest-scoring hit per query
# -------------------------------------------------

echo
echo "Selecting best hits..."

for reference in uniprot TAIR10
do
    sort \
        -k1,1 \
        -k12,12gr \
        -k11,11g \
        "$OUTDIR/blastp_${reference}.tsv" \
        | awk -F '\t' '!seen[$1]++' \
        > "$OUTDIR/blastp_${reference}.besthits.tsv"
done

# -------------------------------------------------
# Summary
# -------------------------------------------------

TOTAL=$(grep -c '^>' "$PROTEINS")

UNIPROT_HITS=$(wc -l < "$OUTDIR/blastp_uniprot.besthits.tsv")
TAIR10_HITS=$(wc -l < "$OUTDIR/blastp_TAIR10.besthits.tsv")

echo
echo "========================================"
echo "BLASTP summary"
echo "========================================"

echo "Total proteins:       $TOTAL"
echo "UniProt best hits:    $UNIPROT_HITS"
echo "TAIR10 best hits:     $TAIR10_HITS"

awk \
    -v hits="$UNIPROT_HITS" \
    -v total="$TOTAL" \
    'BEGIN {
        printf "UniProt hit rate:     %.2f%%\n", 100*hits/total
    }'

awk \
    -v hits="$TAIR10_HITS" \
    -v total="$TOTAL" \
    'BEGIN {
        printf "TAIR10 hit rate:      %.2f%%\n", 100*hits/total
    }'

echo
echo "Output files:"
ls -lh \
    "$OUTDIR/blastp_uniprot.tsv" \
    "$OUTDIR/blastp_uniprot.besthits.tsv" \
    "$OUTDIR/blastp_TAIR10.tsv" \
    "$OUTDIR/blastp_TAIR10.besthits.tsv"

echo
echo "========================================"
echo "Week 4 BLASTP completed"
echo "========================================"
