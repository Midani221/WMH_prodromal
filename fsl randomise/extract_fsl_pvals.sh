#!/usr/bin/env bash

# ------------------------------------------------------------
# Extract summary statistics from FSL randomise outputs
# ------------------------------------------------------------

DESIGNS_DIR="designs"
OUT_CSV="${DESIGNS_DIR}/randomise_corrected_stats.csv"

# If needed on your HPC:
export FSLDIR=/usr/local/fsl-5.0.7
source ${FSLDIR}/etc/fslconf/fsl.sh
export PATH=${FSLDIR}/bin:${PATH}

# ------------------------------------------------------------
# CSV header
# ------------------------------------------------------------

echo "Model,\
TFCE_Pos_Max1mP,TFCE_Pos_MinP,TFCE_Pos_SigVoxels,TFCE_Pos_PeakT,TFCE_Pos_MNI_X,TFCE_Pos_MNI_Y,TFCE_Pos_MNI_Z,\
TFCE_Neg_Max1mP,TFCE_Neg_MinP,TFCE_Neg_SigVoxels,TFCE_Neg_PeakT,TFCE_Neg_MNI_X,TFCE_Neg_MNI_Y,TFCE_Neg_MNI_Z,\
Vox_Pos_Max1mP,Vox_Pos_MinP,Vox_Pos_SigVoxels,Vox_Pos_PeakT,Vox_Pos_MNI_X,Vox_Pos_MNI_Y,Vox_Pos_MNI_Z,\
Vox_Neg_Max1mP,Vox_Neg_MinP,Vox_Neg_SigVoxels,Vox_Neg_PeakT,Vox_Neg_MNI_X,Vox_Neg_MNI_Y,Vox_Neg_MNI_Z" \
> "$OUT_CSV"


# ------------------------------------------------------------
# Function to extract statistics
#
# Arguments:
#   $1 = corrected p image (stored as 1-p)
#   $2 = corresponding tstat image
#
# Output:
# max1mp,minp,nvox,peakT,x,y,z
# ------------------------------------------------------------

extract_stats () {

    corrp_img="$1"
    tstat_img="$2"

    if [[ ! -f "$corrp_img" ]]; then
        echo "NA,NA,NA,NA,NA,NA,NA"
        return
    fi

    # --------------------------------------------------------
    # Maximum corrected 1-p and minimum corrected p
    # --------------------------------------------------------

    max1mp=$(fslstats "$corrp_img" -R | awk '{print $2}')

    minp=$(awk -v x="$max1mp" \
        'BEGIN {printf "%.6f", 1-x}')


    # --------------------------------------------------------
    # Check whether anything survives corrected p < 0.05
    #
    # randomise corrp images contain 1-p:
    # 1-p >= 0.95 corresponds to corrected p <= 0.05
    # --------------------------------------------------------

    is_sig=$(awk -v x="$max1mp" \
        'BEGIN {if (x >= 0.95) print 1; else print 0}')


    if [[ "$is_sig" -eq 0 ]]; then

        echo "${max1mp},${minp},0,NA,NA,NA,NA"
        return

    fi


    # --------------------------------------------------------
    # If corresponding tstat exists:
    #
    # cluster:
    #   - defines significant voxels from corrp image
    #   - threshold = 0.95
    #   - uses tstat image as associated statistic
    #   - reports coordinates in mm
    # --------------------------------------------------------

    if [[ -f "$tstat_img" ]]; then

        cluster_output=$(
            cluster \
                -i "$corrp_img" \
                -t 0.95 \
                -c "$tstat_img" \
                --mm \
                --scalarname="1-p"
        )

        # ----------------------------------------------------
        # Total number of significant voxels
        #
        # Column 2 = number of voxels in each cluster
        # Sum across all significant clusters
        # ----------------------------------------------------

        nvox=$(
            echo "$cluster_output" |
            awk 'NR > 1 {sum += $2} END {print sum+0}'
        )


        # ----------------------------------------------------
        # With -c tstat:
        #
        # columns are:
        # 10 = maximum associated t-statistic
        # 11 = peak X coordinate (mm)
        # 12 = peak Y coordinate (mm)
        # 13 = peak Z coordinate (mm)
        #
        # Find cluster containing largest t-statistic
        # ----------------------------------------------------

        peak_info=$(
            echo "$cluster_output" |
            awk '
                NR > 1 {
                    if (!found || $10 > maxT) {
                        maxT = $10
                        x = $11
                        y = $12
                        z = $13
                        found = 1
                    }
                }
                END {
                    if (found)
                        print maxT, x, y, z
                }
            '
        )

        if [[ -n "$peak_info" ]]; then

            peakT=$(echo "$peak_info" | awk '{print $1}')
            peakX=$(echo "$peak_info" | awk '{print $2}')
            peakY=$(echo "$peak_info" | awk '{print $3}')
            peakZ=$(echo "$peak_info" | awk '{print $4}')

        else

            peakT="NA"
            peakX="NA"
            peakY="NA"
            peakZ="NA"

        fi

    else

        # corrp image exists but corresponding tstat missing
        nvox="NA"
        peakT="NA"
        peakX="NA"
        peakY="NA"
        peakZ="NA"

    fi


    echo "${max1mp},${minp},${nvox},${peakT},${peakX},${peakY},${peakZ}"
}


# ------------------------------------------------------------
# Loop over every model directory
# ------------------------------------------------------------

for directory in "${DESIGNS_DIR}"/*; do

    [[ -d "$directory" ]] || continue

    model=$(basename "$directory")

    echo "Processing: $model"


    # --------------------------------------------------------
    # Corrected probability images
    # --------------------------------------------------------

    tfce_pos="${directory}/model_mask5_tfce_corrp_tstat1.nii"
    tfce_neg="${directory}/model_mask5_tfce_corrp_tstat2.nii"

    vox_pos="${directory}/model_mask5_vox_corrp_tstat1.nii"
    vox_neg="${directory}/model_mask5_vox_corrp_tstat2.nii"


    # --------------------------------------------------------
    # Corresponding raw t-statistic images
    #
    # Both TFCE and voxelwise inference originate from the
    # same randomise tstat1 / tstat2 images.
    # --------------------------------------------------------

    tstat_pos="${directory}/model_mask5_tstat1.nii"
    tstat_neg="${directory}/model_mask5_tstat2.nii"


    # --------------------------------------------------------
    # Extract
    # --------------------------------------------------------

    tfce_pos_stats=$(extract_stats "$tfce_pos" "$tstat_pos")
    tfce_neg_stats=$(extract_stats "$tfce_neg" "$tstat_neg")

    vox_pos_stats=$(extract_stats "$vox_pos" "$tstat_pos")
    vox_neg_stats=$(extract_stats "$vox_neg" "$tstat_neg")


    # --------------------------------------------------------
    # Write one row per model
    # --------------------------------------------------------

    echo "${model},\
${tfce_pos_stats},\
${tfce_neg_stats},\
${vox_pos_stats},\
${vox_neg_stats}" \
    >> "$OUT_CSV"

done


echo ""
echo "=============================================="
echo "Done."
echo "Results saved to:"
echo "$OUT_CSV"
echo "=============================================="