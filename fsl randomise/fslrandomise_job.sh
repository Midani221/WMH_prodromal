#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N randomise
#$ -l h_rt=72:00:00
#$ -l h_vmem=24G
#$ -pe smp 1
#$ -j y

# Load FSL
export FSLDIR=/usr/local/fsl-5.0.7
source "${FSLDIR}/etc/fslconf/fsl.sh"
export PATH="${FSLDIR}/bin:${PATH}"

set -e

# Work in the directory from which qsub was submitted
WORKDIR="${SGE_O_WORKDIR:-$PWD}"
cd "$WORKDIR"

# Input files
INPUT_4D="input_4D.nii.gz"
DESIGN_MAT="design.mat"
DESIGN_CON="design.con"
MASK="WMH_mask_5percent.nii.gz"

# Output prefix
OUTPUT="results/model_mask5"

mkdir -p "$(dirname "$OUTPUT")"

# Confirm that all inputs exist
for FILE in "$INPUT_4D" "$DESIGN_MAT" "$DESIGN_CON" "$MASK"; do
    if [ ! -f "$FILE" ]; then
        echo "ERROR: Cannot find $FILE"
        exit 1
    fi
done

echo "Starting randomise..."
echo "Working directory: $WORKDIR"
echo "Output prefix: $OUTPUT"

randomise \
    -i "$INPUT_4D" \
    -o "$OUTPUT" \
    -d "$DESIGN_MAT" \
    -t "$DESIGN_CON" \
    -m "$MASK" \
    -n 10000 \
    -T \
    -x

echo "Randomise finished successfully."
echo "Results saved with prefix: ${WORKDIR}/${OUTPUT}"