#!/usr/bin/env bash
#$ -S /bin/bash
#$ -N extract_wmh_roi
#$ -cwd
#$ -j y
#$ -l h_vmem=24G
#$ -l h_rt=04:00:00

set -uo pipefail

# Directories containing sub-*/ses-*/output folders, # ADD ALL THE ROOTS OF INTEREST
SEARCH_ROOTS=(
  

)

OUTPUT_DIR="$PWD"

# Load FSL if needed
if ! command -v fslstats >/dev/null 2>&1; then
    export FSLDIR=/usr/local/fsl-5.0.7
    source "${FSLDIR}/etc/fslconf/fsl.sh"
    export PATH="${FSLDIR}/bin:${PATH}"
fi

if ! command -v fslstats >/dev/null 2>&1 ||
   ! command -v fslmaths >/dev/null 2>&1; then
    echo "ERROR: fslstats or fslmaths was not found."
    exit 1
fi

DATESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTCSV="${OUTPUT_DIR}/results2t1roi_combined_volume_${DATESTAMP}.csv"

path_component() {
    local path="$1"
    local prefix="$2"
    local value

    value=$(printf '%s\n' "$path" |
        tr '/' '\n' |
        grep "^${prefix}" |
        tail -n 1)

    echo "${value:-NA}"
}

csv_quote() {
    local value="$1"
    value=${value//\"/\"\"}
    printf '"%s"' "$value"
}

# Validate search roots
VALID_ROOTS=()

for root in "${SEARCH_ROOTS[@]}"; do
    if [[ -d "$root" ]]; then
        VALID_ROOTS+=("$(cd "$root" && pwd -P)")
    else
        echo "WARNING: Directory not found: $root"
    fi
done

if [[ ${#VALID_ROOTS[@]} -eq 0 ]]; then
    echo "ERROR: None of the search directories exist."
    exit 1
fi

# Find results2t1roi_deep files
TMPFILE=$(mktemp "${OUTPUT_DIR}/t1roi_deep_files.XXXXXX")
trap 'rm -f "$TMPFILE"' EXIT

for root in "${VALID_ROOTS[@]}"; do
    find "$root" -type f \
        \( -name "results2t1roi_deep.nii.gz" \
        -o -name "results2t1roi_deep.nii" \) \
        2>/dev/null >> "$TMPFILE"
done

mapfile -t DEEP_FILES < <(sort -u "$TMPFILE")

if [[ ${#DEEP_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No results2t1roi_deep files were found."
    exit 1
fi

echo "subject,session,combined_file,volume_ml" > "$OUTCSV"

echo "Found ${#DEEP_FILES[@]} deep WMH files."

for deep_image in "${DEEP_FILES[@]}"; do
    outdir=$(dirname "$deep_image")

    # Find the matching periventricular file
    if [[ -f "$outdir/results2t1roi_perivent.nii.gz" ]]; then
        perivent_image="$outdir/results2t1roi_perivent.nii.gz"
    elif [[ -f "$outdir/results2t1roi_perivent.nii" ]]; then
        perivent_image="$outdir/results2t1roi_perivent.nii"
    else
        echo "WARNING: Missing periventricular mask in $outdir"
        continue
    fi

    subject=$(path_component "$outdir" "sub-")
    session=$(path_component "$outdir" "ses-")

    combined_image="$outdir/results2t1roi_combined.nii.gz"

    # Add deep and periventricular masks, then ensure binary output
    if ! fslmaths "$deep_image" \
        -add "$perivent_image" \
        -bin "$combined_image"; then
        echo "WARNING: Could not create combined mask for $subject $session"
        continue
    fi

    # fslstats -V returns: voxel_count volume_mm3
    read -r voxel_count volume_mm3 <<< \
        "$(fslstats "$combined_image" -V)"

    volume_ml=$(awk -v volume="$volume_mm3" '
        BEGIN {
            if (volume ~ /^[0-9.eE+-]+$/) {
                printf "%.6f", volume / 1000
            } else {
                print "NA"
            }
        }
    ')

    row="$(csv_quote "$subject"),$(csv_quote "$session")"
    row+=",$(csv_quote "$combined_image"),$volume_ml"

    echo "$row" >> "$OUTCSV"
    echo "Processed: $subject $session - ${volume_ml} mL"
done

echo "CSV saved to: $OUTCSV"