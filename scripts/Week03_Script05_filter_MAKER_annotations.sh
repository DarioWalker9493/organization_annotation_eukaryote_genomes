#!/bin/bash
#SBATCH --job-name=Maker_filter
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dario.walker@students.unibe.ch

set -euo pipefail

COURSEDIR=/data/courses/genome_annotation_resources/2026
WORKDIR=/data/users/dwalker/organization_annotation_eukaryote_genomes
OUTDIR="$WORKDIR/output/MAKER"
FINALDIR="$OUTDIR/final"

MAKER_SIF="$COURSEDIR/containers/MAKER_3.01.03.sif"
ANALYSIS_SIF="$COURSEDIR/containers/Annotation_analysis_R4.4.sif"

protein="assembly.all.maker.proteins.fasta"
transcript="assembly.all.maker.transcripts.fasta"
gff="assembly.all.maker.noseq.gff"

mkdir -p \
    "$FINALDIR/tmp" \
    "$WORKDIR/logs"

cd "$FINALDIR"

echo "========================================"
echo "Week 3 - MAKER annotation refinement"
echo "========================================"

# -------------------------------------------------
# Preflight
# -------------------------------------------------

for FILE in \
    "${gff}.renamed.gff" \
    "${protein}.renamed.fasta" \
    "${transcript}.renamed.fasta" \
    output.iprscan
do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Required file missing or empty:"
        echo "$FILE"
        exit 1
    fi
done

echo
echo "Inputs:"
ls -lh \
    "${gff}.renamed.gff" \
    "${protein}.renamed.fasta" \
    "${transcript}.renamed.fasta" \
    output.iprscan

# -------------------------------------------------
# 1. Add InterProScan / Pfam annotations to GFF
# -------------------------------------------------

echo
echo "Adding InterProScan annotations to GFF..."

apptainer exec \
    --cleanenv \
    --env TMPDIR=/tmp \
    --bind "$OUTDIR/tmp":/tmp \
    --bind "$COURSEDIR" \
    --bind "$WORKDIR" \
    --pwd "$FINALDIR" \
    "$MAKER_SIF" \
    ipr_update_gff \
    "${gff}.renamed.gff" \
    output.iprscan \
    > "${gff}.renamed.iprscan.gff"

# -------------------------------------------------
# 2. AED cumulative distribution
# -------------------------------------------------

echo
echo "Calculating AED distribution..."

apptainer exec \
    --cleanenv \
    --env TMPDIR=/tmp \
    --bind "$OUTDIR/tmp":/tmp \
    --bind "$COURSEDIR" \
    --bind "$WORKDIR" \
    --pwd "$FINALDIR" \
    "$MAKER_SIF" \
    perl /usr/local/bin/AED_cdf_generator.pl \
    -b 0.025 \
    "${gff}.renamed.gff" \
    > assembly.all.maker.renamed.gff.AED.txt

# -------------------------------------------------
# 3. Quality filter
# AED < 1 and/or Pfam domain
# -------------------------------------------------

echo
echo "Applying MAKER quality filter..."

apptainer exec \
    --cleanenv \
    --env TMPDIR=/tmp \
    --bind "$OUTDIR/tmp":/tmp \
    --bind "$COURSEDIR" \
    --bind "$WORKDIR" \
    --pwd "$FINALDIR" \
    "$MAKER_SIF" \
    perl /usr/local/bin/quality_filter.pl \
    -s \
    "${gff}.renamed.iprscan.gff" \
    > "${gff}_iprscan_quality_filtered.gff"

# -------------------------------------------------
# 4. Keep only gene-related features
# -------------------------------------------------

echo
echo "Extracting gene features..."

grep -P \
    "\tgene\t|\tCDS\t|\texon\t|\tfive_prime_UTR\t|\tthree_prime_UTR\t|\tmRNA\t" \
    "${gff}_iprscan_quality_filtered.gff" \
    > filtered.genes.renamed.gff3

sed -i '1i##gff-version 3' filtered.genes.renamed.gff3

# -------------------------------------------------
# 5. Extract retained mRNA IDs
# -------------------------------------------------

echo
echo "Extracting retained mRNA IDs..."

grep -P "\tmRNA\t" filtered.genes.renamed.gff3 \
    | awk '{print $9}' \
    | cut -d ';' -f1 \
    | sed 's/ID=//g' \
    > list.txt

# -------------------------------------------------
# 6. Filter transcript FASTA
# -------------------------------------------------

echo
echo "Filtering transcript FASTA..."

apptainer exec \
    --cleanenv \
    --bind "$COURSEDIR" \
    --bind "$WORKDIR" \
    --pwd "$FINALDIR" \
    "$ANALYSIS_SIF" \
    seqkit grep \
    -f list.txt \
    "${transcript}.renamed.fasta" \
    -o "${transcript}.renamed.filtered.fasta"

# -------------------------------------------------
# 7. Filter protein FASTA
# -------------------------------------------------

echo
echo "Filtering protein FASTA..."

apptainer exec \
    --cleanenv \
    --bind "$COURSEDIR" \
    --bind "$WORKDIR" \
    --pwd "$FINALDIR" \
    "$ANALYSIS_SIF" \
    seqkit grep \
    -f list.txt \
    "${protein}.renamed.fasta" \
    -o "${protein}.renamed.filtered.fasta"

# -------------------------------------------------
# Summary
# -------------------------------------------------

echo
echo "========================================"
echo "Annotation filtering summary"
echo "========================================"

echo
echo "Feature types retained:"
cut -f3 filtered.genes.renamed.gff3 \
    | grep -v '^#' \
    | sort \
    | uniq -c

echo
echo "Gene models:"
awk -F '\t' '$3=="gene"{n++} END{print n+0}' \
    filtered.genes.renamed.gff3

echo "mRNA models:"
awk -F '\t' '$3=="mRNA"{n++} END{print n+0}' \
    filtered.genes.renamed.gff3

echo "Filtered proteins:"
grep -c '^>' "${protein}.renamed.filtered.fasta"

echo "Filtered transcripts:"
grep -c '^>' "${transcript}.renamed.filtered.fasta"

echo
echo "Files:"
ls -lh \
    "${gff}.renamed.iprscan.gff" \
    assembly.all.maker.renamed.gff.AED.txt \
    "${gff}_iprscan_quality_filtered.gff" \
    filtered.genes.renamed.gff3 \
    "${protein}.renamed.filtered.fasta" \
    "${transcript}.renamed.filtered.fasta"

echo
echo "========================================"
echo "MAKER annotation refinement completed"
echo "========================================"
