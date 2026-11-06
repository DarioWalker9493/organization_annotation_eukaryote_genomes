#!/bin/bash
#SBATCH --job-name=TEsorter_lib
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dario.walker@students.unibe.ch

set -euo pipefail

COURSEDIR=/data/courses/genome_annotation_resources/2026
WORKDIR=/data/users/dwalker/organization_annotation_eukaryote_genomes

TELIB="$WORKDIR/output/EDTA/assembly.fasta.mod.EDTA.TElib.fa"
OUTDIR="$WORKDIR/output/TEsorter_library"

mkdir -p "$OUTDIR"

echo "========================================"
echo "Week 2 - Task 3"
echo "Library-level Copia/Gypsy classification"
echo "========================================"

# -------------------------------------------------
# Preflight
# -------------------------------------------------

if [[ ! -s "$TELIB" ]]; then
    echo "ERROR: EDTA TE library not found:"
    echo "$TELIB"
    exit 1
fi

echo "Input TE library:"
ls -lh "$TELIB"

# -------------------------------------------------
# Extract Copia sequences
# -------------------------------------------------

echo
echo "Extracting Copia sequences..."

apptainer exec -C \
    --bind "$WORKDIR" \
    --bind "$COURSEDIR" \
    "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
    seqkit grep -r -p "Copia" \
    "$TELIB" \
    > "$OUTDIR/Copia_sequences.fa"

# -------------------------------------------------
# Extract Gypsy sequences
# -------------------------------------------------

echo
echo "Extracting Gypsy sequences..."

apptainer exec -C \
    --bind "$WORKDIR" \
    --bind "$COURSEDIR" \
    "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
    seqkit grep -r -p "Gypsy" \
    "$TELIB" \
    > "$OUTDIR/Gypsy_sequences.fa"

# -------------------------------------------------
# Check extracted datasets
# -------------------------------------------------

echo
echo "Extracted sequence counts:"

echo -n "Copia: "
grep -c '^>' "$OUTDIR/Copia_sequences.fa"

echo -n "Gypsy: "
grep -c '^>' "$OUTDIR/Gypsy_sequences.fa"

cd "$OUTDIR"

# -------------------------------------------------
# Run TEsorter on Copia
# -------------------------------------------------

echo
echo "Running TEsorter on Copia..."

apptainer exec -C \
    --bind "$WORKDIR" \
    --bind "$COURSEDIR" \
    "$COURSEDIR/containers/TEsorter_1.3.0.sif" \
    TEsorter \
    Copia_sequences.fa \
    -db rexdb-plant \
    -p "$SLURM_CPUS_PER_TASK"

# -------------------------------------------------
# Run TEsorter on Gypsy
# -------------------------------------------------

echo
echo "Running TEsorter on Gypsy..."

apptainer exec -C \
    --bind "$WORKDIR" \
    --bind "$COURSEDIR" \
    "$COURSEDIR/containers/TEsorter_1.3.0.sif" \
    TEsorter \
    Gypsy_sequences.fa \
    -db rexdb-plant \
    -p "$SLURM_CPUS_PER_TASK"

# -------------------------------------------------
# Final check
# -------------------------------------------------

echo
echo "Generated classification files:"

ls -lh \
    Copia_sequences.fa.rexdb-plant.cls.tsv \
    Gypsy_sequences.fa.rexdb-plant.cls.tsv

echo
echo "========================================"
echo "Task 3 TEsorter classification completed"
echo "========================================"