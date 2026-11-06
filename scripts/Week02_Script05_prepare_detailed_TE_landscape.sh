#!/bin/bash
#SBATCH --job-name=TE_landprep
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dario.walker@students.unibe.ch

set -euo pipefail

COURSEDIR=/data/courses/genome_annotation_resources/2026
WORKDIR=/data/users/dwalker/organization_annotation_eukaryote_genomes

ANNO_DIR="$WORKDIR/output/EDTA/assembly.fasta.mod.EDTA.anno"

RM_OUT="$ANNO_DIR/assembly.fasta.mod.out"
RM_NONZERO="$ANNO_DIR/assembly.fasta.mod.nonzero.out"

echo "========================================"
echo "Week 2 - Task 4"
echo "Detailed TE landscape preparation"
echo "========================================"

# -------------------------------------------------
# Preflight
# -------------------------------------------------

if [[ ! -s "$RM_OUT" ]]; then
    echo "ERROR: RepeatMasker output not found:"
    echo "$RM_OUT"
    exit 1
fi

echo
echo "Detailed RepeatMasker input:"
ls -lh "$RM_OUT"

# -------------------------------------------------
# Remove zero-divergence matches
# -------------------------------------------------

echo
echo "Removing records with 0% divergence..."

awk 'NR == 1 || $2 != 0' \
    "$RM_OUT" \
    > "$RM_NONZERO"

echo
echo "Original lines:"
wc -l "$RM_OUT"

echo
echo "Non-zero divergence lines:"
wc -l "$RM_NONZERO"

# -------------------------------------------------
# Parse divergence landscape
# -------------------------------------------------

echo
echo "Running parseRM.pl..."

cd "$ANNO_DIR"

apptainer exec -C \
    --bind "$WORKDIR" \
    --bind "$COURSEDIR" \
    --pwd "$ANNO_DIR" \
    "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
    perl "$COURSEDIR/scripts/05-parseRM.pl" \
    -i "$(basename "$RM_NONZERO")" \
    -l 50,1 \
    -v

# -------------------------------------------------
# Check detailed outputs
# -------------------------------------------------

RNAME="$RM_NONZERO.landscape.Div.Rname.tab"
RFAM="$RM_NONZERO.landscape.Div.Rfam.tab"
RCLASS="$RM_NONZERO.landscape.Div.Rclass.tab"

echo
echo "Generated landscape tables:"

ls -lh \
    "$RNAME" \
    "$RFAM" \
    "$RCLASS"

echo
echo "----------------------------------------"
echo "Rname preview"
echo "----------------------------------------"

if [[ -s "$RNAME" ]]; then
    head -n 12 "$RNAME"
else
    echo "WARNING: Rname table is empty."
fi

echo
echo "----------------------------------------"
echo "Rfam preview"
echo "----------------------------------------"

if [[ -s "$RFAM" ]]; then
    head -n 12 "$RFAM"
else
    echo "WARNING: Rfam table is empty."
fi

echo
echo "========================================"
echo "Detailed TE landscape preparation done"
echo "========================================"