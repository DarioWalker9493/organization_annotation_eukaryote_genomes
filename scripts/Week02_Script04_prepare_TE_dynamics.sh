#!/bin/bash
#SBATCH --job-name=TE_dynamics
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

RM_OUT="$ANNO_DIR/assembly.fasta.mod.EDTA.RM.out"
RM_NONZERO="$ANNO_DIR/assembly.fasta.mod.EDTA.RM.nonzero.out"

echo "========================================"
echo "Week 2 - Task 4"
echo "TE divergence landscape preparation"
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
echo "RepeatMasker input:"
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
echo "Filtered RepeatMasker file:"
ls -lh "$RM_NONZERO"

echo
echo "Original lines:"
wc -l "$RM_OUT"

echo "Non-zero divergence lines:"
wc -l "$RM_NONZERO"

# -------------------------------------------------
# Parse RepeatMasker divergence
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
# Show generated files
# -------------------------------------------------

echo
echo "Files generated/modified in annotation directory:"
find "$ANNO_DIR" \
    -maxdepth 1 \
    -type f \
    -newer "$RM_NONZERO" \
    -printf "%f\t%k KB\n" \
    | sort

echo
echo "========================================"
echo "TE divergence parsing completed"
echo "========================================"