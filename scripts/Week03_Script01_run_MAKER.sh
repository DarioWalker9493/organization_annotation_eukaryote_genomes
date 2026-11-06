#!/bin/bash
#SBATCH --job-name=Maker_gene_annotation
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dario.walker@students.unibe.ch

set -euo pipefail

COURSEDIR=/data/courses/genome_annotation_resources/2026
WORKDIR=/data/users/dwalker/organization_annotation_eukaryote_genomes
OUTDIR="$WORKDIR/output/MAKER"

MAKER_SIF="$COURSEDIR/containers/MAKER_3.01.03.sif"
AUGUSTUS_CONFIG_PATH="$COURSEDIR/tools/MAKER/augustus_config"

echo "========================================"
echo "Week 3 - Gene Annotation with MAKER"
echo "========================================"
echo "Workdir: $WORKDIR"
echo "Output:  $OUTDIR"
echo "CPUs:    $SLURM_CPUS_PER_TASK"
echo "========================================"

# -------------------------------------------------
# Preflight
# -------------------------------------------------

GENOME="$WORKDIR/data/assemblies/Lu1_Flye.fasta"
TRINITY="$WORKDIR/data/transcripts/Lu1_Trinity.fasta"
TELIB="$WORKDIR/output/EDTA/assembly.fasta.mod.EDTA.TElib.fa"

for FILE in \
    "$GENOME" \
    "$TRINITY" \
    "$TELIB" \
    "$WORKDIR/maker_opts.ctl" \
    "$WORKDIR/maker_bopts.ctl" \
    "$WORKDIR/maker_exe.ctl" \
    "$WORKDIR/maker_evm.ctl"
do
    if [[ ! -e "$FILE" ]]; then
        echo "ERROR: Required file not found:"
        echo "$FILE"
        exit 1
    fi
done

echo
echo "Input files:"
ls -lh \
    "$GENOME" \
    "$TRINITY" \
    "$TELIB"

# -------------------------------------------------
# Prepare MAKER output directory
# -------------------------------------------------

unset SINGULARITY_BIND
unset APPTAINER_BIND
unset SINGULARITY_TMPDIR

mkdir -p \
    "$OUTDIR/tmp" \
    "$WORKDIR/logs"

# Copy the edited control files into the MAKER run directory
cp \
    "$WORKDIR/maker_opts.ctl" \
    "$WORKDIR/maker_bopts.ctl" \
    "$WORKDIR/maker_exe.ctl" \
    "$WORKDIR/maker_evm.ctl" \
    "$OUTDIR/"

cd "$OUTDIR"

echo
echo "MAKER control files:"
ls -lh maker_*.ctl

echo
echo "Starting MAKER..."
echo

# -------------------------------------------------
# Run MAKER with MPI
# -------------------------------------------------

apptainer exec \
    --cleanenv \
    --env AUGUSTUS_CONFIG_PATH="$AUGUSTUS_CONFIG_PATH" \
    --bind /data \
    --bind "$WORKDIR" \
    --bind "$OUTDIR/tmp":/TMP \
    --pwd "$OUTDIR" \
    "$MAKER_SIF" \
    mpiexec \
    -launcher fork \
    -n "$SLURM_CPUS_PER_TASK" \
    maker \
    -mpi \
    --ignore_nfs_tmp \
    -TMP /TMP \
    maker_opts.ctl \
    maker_bopts.ctl \
    maker_exe.ctl \
    maker_evm.ctl

echo
echo "========================================"
echo "MAKER run finished"
echo "========================================"

# -------------------------------------------------
# Basic output check
# -------------------------------------------------

echo
echo "MAKER output directories:"
find "$OUTDIR" \
    -maxdepth 1 \
    -type d \
    -name "*.maker.output" \
    -print

echo
echo "Datastore indexes:"
find "$OUTDIR" \
    -maxdepth 2 \
    -name "*_master_datastore_index.log" \
    -print

echo
echo "========================================"
echo "Week 3 MAKER job completed"
echo "========================================"