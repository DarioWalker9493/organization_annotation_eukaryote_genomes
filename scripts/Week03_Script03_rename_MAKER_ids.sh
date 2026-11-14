#!/bin/bash
#SBATCH --job-name=Maker_rename
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
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

MAKER_SIF="$COURSEDIR/containers/MAKER_3.01.03.sif"

protein="assembly.all.maker.proteins.fasta"
transcript="assembly.all.maker.transcripts.fasta"
gff="assembly.all.maker.noseq.gff"

PREFIX="LU1"

mkdir -p \
    "$FINALDIR/tmp" \
    "$WORKDIR/logs"

echo "========================================"
echo "Week 3 - Rename MAKER gene models"
echo "========================================"
echo "Prefix: $PREFIX"
echo

# -------------------------------------------------
# Copy original merged annotation files
# -------------------------------------------------

cp "$OUTDIR/$gff" \
   "$FINALDIR/${gff}.renamed.gff"

cp "$OUTDIR/$protein" \
   "$FINALDIR/${protein}.renamed.fasta"

cp "$OUTDIR/$transcript" \
   "$FINALDIR/${transcript}.renamed.fasta"

cd "$FINALDIR"

echo "Original MAKER IDs:"
cut -f9 "${gff}.renamed.gff" \
    | cut -d ';' -f1 \
    | sed 's/ID=//g' \
    | head -n 10 || true

# -------------------------------------------------
# Generate ID mapping
# -------------------------------------------------

echo
echo "Generating ID map..."

apptainer exec \
    --cleanenv \
    --env TMPDIR=/tmp \
    --bind "$OUTDIR/tmp":/tmp \
    --bind "$COURSEDIR" \
    --bind "$WORKDIR" \
    --pwd "$FINALDIR" \
    "$MAKER_SIF" \
    maker_map_ids \
    --prefix "$PREFIX" \
    --justify 7 \
    "${gff}.renamed.gff" \
    > id.map

echo
echo "ID map preview:"
head -n 10 id.map

# -------------------------------------------------
# Rename IDs in GFF
# -------------------------------------------------

echo
echo "Renaming GFF IDs..."

apptainer exec \
    --cleanenv \
    --env TMPDIR=/tmp \
    --bind "$OUTDIR/tmp":/tmp \
    --bind "$COURSEDIR" \
    --bind "$WORKDIR" \
    --pwd "$FINALDIR" \
    "$MAKER_SIF" \
    map_gff_ids \
    id.map \
    "${gff}.renamed.gff"

# -------------------------------------------------
# Rename protein FASTA IDs
# -------------------------------------------------

echo
echo "Renaming protein IDs..."

apptainer exec \
    --cleanenv \
    --env TMPDIR=/tmp \
    --bind "$OUTDIR/tmp":/tmp \
    --bind "$COURSEDIR" \
    --bind "$WORKDIR" \
    --pwd "$FINALDIR" \
    "$MAKER_SIF" \
    map_fasta_ids \
    id.map \
    "${protein}.renamed.fasta"

# -------------------------------------------------
# Rename transcript FASTA IDs
# -------------------------------------------------

echo
echo "Renaming transcript IDs..."

apptainer exec \
    --cleanenv \
    --env TMPDIR=/tmp \
    --bind "$OUTDIR/tmp":/tmp \
    --bind "$COURSEDIR" \
    --bind "$WORKDIR" \
    --pwd "$FINALDIR" \
    "$MAKER_SIF" \
    map_fasta_ids \
    id.map \
    "${transcript}.renamed.fasta"

# -------------------------------------------------
# Checks
# -------------------------------------------------

echo
echo "Renamed GFF IDs:"
cut -f9 "${gff}.renamed.gff" \
    | cut -d ';' -f1 \
    | sed 's/ID=//g' \
    | head -n 10 || true

echo
echo "Renamed protein IDs:"
grep '^>' "${protein}.renamed.fasta" | head -n 10 || true

echo
echo "Renamed transcript IDs:"
grep '^>' "${transcript}.renamed.fasta" | head -n 10 || true

echo
echo "Counts after renaming:"

awk -F '\t' '$3=="gene"{n++} END{print "Genes:",n+0}' \
    "${gff}.renamed.gff"

awk -F '\t' '$3=="mRNA"{n++} END{print "mRNAs:",n+0}' \
    "${gff}.renamed.gff"

echo -n "Proteins: "
grep -c '^>' "${protein}.renamed.fasta"

echo -n "Transcripts: "
grep -c '^>' "${transcript}.renamed.fasta"

echo
echo "========================================"
echo "MAKER ID renaming completed"
echo "========================================"
