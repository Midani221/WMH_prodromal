This is most of the code that I have used during my thesis for my MSc in Advanced Neuroimaging. I have created this to be part of my submission.
The analysis was done in 4 main steps, first of all, the primary imaging patient sheets were prepared on excel and raw imaging data (DICOM format) was downloaded
from the repository.
This was followed by these steps:
1. Bidsification: When downloading from PPMI, the data is split unevenlly across two or more folders, with patients having images in either folders. You need to:
   - Run Mover_final.py to fix image locations
   - v_fixer to fix nomenclature for easy saving
   - bidsify_multiple to bidsify each session for each participant and save in bids format across all.
   - If bidsify_multiple does not work for some case, you can use manual conversion to NIFTI using MRIcrogl.exe.
2. Pipeline run:
   - Pipeline run using bash on the HPC.
   - Collate the images if you need to download them for QC, then run the shiny app for QC to make it faster. There is QC integrated within the pipeline but it can
      be slow.
   - Extract volume using fslstats. The attached script works for total native volume from native space.
3. Clinical Analysis:
   - This is a bit disorganised due to multiple runs and diagnostics with AI, overall, sheets should be downloaded from the repository and the code will clean it
     using predetermined variable names.
   - FDA_version2 will run automatically to carry out functional principal component analysis using PACE. FDA original does it manually.
   - global_wmh_models gives outputs for all linear models.
   - Graphing can be done separately depending on question of interest.
4. FSL randomise
   - Directory files is generated from file of interest that has subject and session names
   - design matrix is generated from a corresponding file with the outcomes and demographics of interest, it also generates the unique text files for each variable
     to account for missingness. You then need to turn it into mat file, show in extra steps.
   - You need to retrieve the files and make the 4D image using fsl merge, shown in extra steps
   - You need to create a mask, 5% mask example is shown.
   - Run fslrandomise_job.sh job. Decide on the computing requirement depending on number of files in 4D image. ~900 require at least 32 gbs.

Creator: To be added later
