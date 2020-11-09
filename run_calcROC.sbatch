#!/bin/bash
#SBATCH --job-name=calcROC
#SBATCH --output=calc-ROC_%A_%a.out
#SBATCH --error=calc-ROC_%A_%a.err
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16000
#SBATCH --time=1000:00:00
#SBATCH --partition=2018allq

# Define and create a unique scratch directory for this job

pwd; hostname; date
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
echo "SLURM_JOBID="$SLURM_JOBID

cp ${SLURM_SUBMIT_DIR}/calcROC.R ${SCRATCH_DIRECTORY}
cp ~/data/ENA_SRP044740/tidydata/*rds ${SCRATCH_DIRECTORY}

# Call script
Rscript --vanilla --quiet --slave calcROC.R ${feature}

# After the job is done we copy our output back to our directory
mkdir -p ${SLURM_SUBMIT_DIR}/feature-imp/${feature}_all-samples_${SLURM_JOBID}/ && cp ${SCRATCH_DIRECTORY}/*.RDS ${SLURM_SUBMIT_DIR}/feature-imp/${feature}_all-samples_${SLURM_JOBID}/
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}

exit 0
