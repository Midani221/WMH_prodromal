#!/bin/bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Section of the script that sets common variables we need for all cluster jobs         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Job name
MYJOBNAME="prodromal3.1_v1.1.1" ;

# Job base directory from where all else happens
export BASE= "BASE PATH of files" ;

# Log directory for standard output and errors
LOGDIR=${BASE}/logs/${MYJOBNAME} ;

# Set this LOGDIR as a directory and create directory if it doesnt exist
[ -d ${LOGDIR} ] || mkdir -p ${LOGDIR} ; chmod 755 ${LOGDIR} ;

# Set a input directory
export SUBJECTS_DIR=${BASE}/PATH to PRODROMAL/data ; # Adjust to put correct path of prodromal folder

# Create a variable with our number of subjects to run
totaltasks=$(ls ${SUBJECTS_DIR} | wc -l) ;

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                Section of the script that sends individual jobs to the cluster              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Set out QSUB SGE block
qsub -N ${MYJOBNAME} -e ${LOGDIR} -o ${LOGDIR} -t 1-${totaltasks} -wd ${BASE} -j y <<'EOT'
#$ -S /bin/bash
#$ -V
#$ -l h_vmem=24G,hostslots=1


# Load in our subjects to a variable
SUBJECTS=($(ls -1d ${SUBJECTS_DIR}/sub-*)) ;

# Pass these subjects to individual variables for use on cluster
SUBJECT=${SUBJECTS[${SGE_TASK_ID}-1]} ;

# Get a base/subject name from our subject image 
SUBID=$(basename ${SUBJECT}| awk -F. '{print $1}') ;

# Print our subject ID
echo -e "Subject ID: ${SUBID}\n" ;

# Print our input and output folders
echo -e "Input: ${SUBJECT}\n" ;

# Write out our WML apptainer command
# CMD="/usr/local/apptainer-1.3.4/bin/apptainer run --bind "${SUBJECTS_DIR}":/data ${JOBBASEDIR}/enigma-pd-wml-1.1.1.sif -s ${SUBID}"

CMD="/usr/local/apptainer-1.3.4/bin/apptainer run --bind "${SUBJECTS_DIR}":/data 
${BASE}/enigma-pd-wml-1.1.1.sif -s ${SUBID}"

# Print our command we want and execute it
echo "Running: ${CMD}" && ${CMD};

# Change permissions
chmod -R 755 ${SUBJECTS_DIR}/derivatives ;

# Print that job is complete
echo "${MYJOBNAME} complete"

# End of SGE block
EOT
