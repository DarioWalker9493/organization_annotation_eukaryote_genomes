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

MAKER_SIF="$COURSEDIR/containers/MAKER_3.01.03.sif"

UNIPROT="$COURSEDIR/data/uniprot_viridiplantae_reviewed.fa"
BESTHITS="$OUTDIR/blastp_uniprot.besthits.tsv"

PROTEINS="$FINALDIR/assembly.all.maker.proteins.fasta.renamed.filtered.fasta"

# Use the complete MAKER-produced GFF here.
BASE_GFF="$FINALDIR/assembly.all.maker.noseq.gff.renamed.gff"

# InterProScan output generated in Week 3.
IPRSCAN="$FINALDIR/output.iprscan"

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
    "$UNIPROT" \
    "$BESTHITS" \
    "$PROTEINS" \
    "$BASE_GFF" \
    "$IPRSCAN"
do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Missing or empty required file:"
        echo "$FILE"
        exit 1
    fi
done

echo
echo "Inputs:"
ls -lh \
    "$UNIPROT" \
    "$BESTHITS" \
    "$PROTEINS" \
    "$BASE_GFF" \
    "$IPRSCAN"


# ============================================================================
# 1. Add UniProt functions to filtered protein FASTA
# ============================================================================

echo
echo "Adding UniProt functions to filtered protein FASTA..."

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
    > "$OUTDIR/maker_proteins.filtered.Uniprot.fasta"


# ============================================================================
# 2. Add UniProt functions to complete renamed MAKER GFF
# ============================================================================

echo
echo "Adding UniProt functions to complete MAKER GFF..."

apptainer exec \
    --cleanenv \
    --env TMPDIR=/tmp \
    --bind /data \
    --bind "$WORKDIR" \
    --bind "$TMPDIR":/tmp \
    --pwd "$OUTDIR" \
    "$MAKER_SIF" \
    maker_functional_gff \
    "$UNIPROT" \
    "$BESTHITS" \
    "$BASE_GFF" \
    > "$OUTDIR/assembly.all.maker.renamed.Uniprot.gff"


# ============================================================================
# 3. Add InterProScan / Pfam / GO information
# ============================================================================

echo
echo "Adding InterProScan / Pfam / GO annotations..."

apptainer exec \
    --cleanenv \
    --env TMPDIR=/tmp \
    --bind /data \
    --bind "$WORKDIR" \
    --bind "$TMPDIR":/tmp \
    --pwd "$OUTDIR" \
    "$MAKER_SIF" \
    ipr_update_gff \
    "$OUTDIR/assembly.all.maker.renamed.Uniprot.gff" \
    "$IPRSCAN" \
    > "$OUTDIR/assembly.all.maker.renamed.Uniprot.iprscan.gff"


# ============================================================================
# 4. Apply same Week-3 quality filter
# ============================================================================

echo
echo "Applying MAKER quality filter..."

apptainer exec \
    --cleanenv \
    --env TMPDIR=/tmp \
    --bind /data \
    --bind "$WORKDIR" \
    --bind "$TMPDIR":/tmp \
    --pwd "$OUTDIR" \
    "$MAKER_SIF" \
    perl /usr/local/bin/quality_filter.pl \
    -s \
    "$OUTDIR/assembly.all.maker.renamed.Uniprot.iprscan.gff" \
    > "$OUTDIR/assembly.all.maker.renamed.Uniprot.iprscan.quality_filtered.gff"


# ============================================================================
# 5. Keep final gene-model features
# ============================================================================

echo
echo "Extracting final gene features..."

grep -P \
    "\tgene\t|\tmRNA\t|\tCDS\t|\texon\t|\tfive_prime_UTR\t|\tthree_prime_UTR\t" \
    "$OUTDIR/assembly.all.maker.renamed.Uniprot.iprscan.quality_filtered.gff" \
    > "$OUTDIR/filtered.genes.Uniprot.gff3"

sed -i '1i##gff-version 3' \
    "$OUTDIR/filtered.genes.Uniprot.gff3"


# ============================================================================
# 6. Summary
# ============================================================================

echo
echo "========================================"
echo "Functional annotation summary"
echo "========================================"

echo -n "Input filtered proteins: "
grep -c '^>' "$PROTEINS"

echo -n "Functional FASTA proteins: "
grep -c '^>' "$OUTDIR/maker_proteins.filtered.Uniprot.fasta"

echo -n "Final genes: "
awk -F '\t' '$3=="gene"{n++} END{print n+0}' \
    "$OUTDIR/filtered.genes.Uniprot.gff3"

echo -n "Final mRNAs: "
awk -F '\t' '$3=="mRNA"{n++} END{print n+0}' \
    "$OUTDIR/filtered.genes.Uniprot.gff3"

echo -n "GFF features with UniProt Note: "
grep -c 'Note=' \
    "$OUTDIR/filtered.genes.Uniprot.gff3" || true

echo -n "GFF features with Pfam Dbxref: "
grep -c 'Pfam:' \
    "$OUTDIR/filtered.genes.Uniprot.gff3" || true

echo -n "GFF features with GO terms: "
grep -c 'Ontology_term=' \
    "$OUTDIR/filtered.genes.Uniprot.gff3" || true

echo
echo "Outputs:"

ls -lh \
    "$OUTDIR/maker_proteins.filtered.Uniprot.fasta" \
    "$OUTDIR/assembly.all.maker.renamed.Uniprot.gff" \
    "$OUTDIR/assembly.all.maker.renamed.Uniprot.iprscan.gff" \
    "$OUTDIR/assembly.all.maker.renamed.Uniprot.iprscan.quality_filtered.gff" \
    "$OUTDIR/filtered.genes.Uniprot.gff3"

echo
echo "========================================"
echo "Week 4 functional annotation completed"
echo "========================================"

