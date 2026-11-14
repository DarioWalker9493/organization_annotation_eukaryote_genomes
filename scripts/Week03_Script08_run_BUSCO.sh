#!/bin/bash
#SBATCH --job-name=BUSCO_MAKER
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dario.walker@students.unibe.ch

set -euo pipefail

WORKDIR=/data/users/dwalker/organization_annotation_eukaryote_genomes
FINALDIR="$WORKDIR/output/MAKER/final"

protein="assembly.all.maker.proteins.fasta"
transcript="assembly.all.maker.transcripts.fasta"

PROTEIN_INPUT="${protein}.renamed.longest.fasta"
TRANSCRIPT_INPUT="${transcript}.renamed.longest.fasta"

mkdir -p "$WORKDIR/logs"

cd "$FINALDIR"

echo "========================================"
echo "Week 3 - BUSCO annotation assessment"
echo "========================================"

for FILE in "$PROTEIN_INPUT" "$TRANSCRIPT_INPUT"
do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Missing or empty input:"
        echo "$FILE"
        exit 1
    fi
done

echo
echo "Inputs:"
ls -lh "$PROTEIN_INPUT" "$TRANSCRIPT_INPUT"

echo
echo -n "Protein sequences: "
grep -c '^>' "$PROTEIN_INPUT"

echo -n "Transcript sequences: "
grep -c '^>' "$TRANSCRIPT_INPUT"

echo
echo "Loading BUSCO..."

module load BUSCO/5.4.2-foss-2021a
export PYTHONNOUSERSITE=1

echo
echo "BUSCO version:"
busco --version

# -------------------------------------------------
# Protein BUSCO
# -------------------------------------------------

echo
echo "========================================"
echo "Running BUSCO on longest proteins"
echo "========================================"

rm -rf busco_proteins

busco \
    -i "$PROTEIN_INPUT" \
    -l brassicales_odb10 \
    -o busco_proteins \
    -m proteins \
    -c "$SLURM_CPUS_PER_TASK"

# -------------------------------------------------
# Transcript BUSCO
# -------------------------------------------------

echo
echo "========================================"
echo "Running BUSCO on longest transcripts"
echo "========================================"

rm -rf busco_transcripts

busco \
    -i "$TRANSCRIPT_INPUT" \
    -l brassicales_odb10 \
    -o busco_transcripts \
    -m transcriptome \
    -c "$SLURM_CPUS_PER_TASK"

# -------------------------------------------------
# Summary
# -------------------------------------------------

echo
echo "========================================"
echo "BUSCO summaries"
echo "========================================"

echo
echo "--- Proteins ---"
find busco_proteins \
    -maxdepth 2 \
    -name 'short_summary*.txt' \
    -exec cat {} \;

echo
echo "--- Transcripts ---"
find busco_transcripts \
    -maxdepth 2 \
    -name 'short_summary*.txt' \
    -exec cat {} \;

echo
echo "========================================"
echo "BUSCO assessment completed"
echo "========================================"
