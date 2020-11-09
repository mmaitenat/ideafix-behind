#!/bin/bash
#SBATCH --job-name=so
#SBATCH --output=so-models_%A_%a.out
#SBATCH --error=so-models_%A_%a.err
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30000
#SBATCH --time=1000:00:00
#SBATCH --partition=2018allq


# Define and create a unique scratch directory for this job

pwd; hostname; date
echo ${model}
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}

cp ${SLURM_SUBMIT_DIR}/trainTestModels.R ${SCRATCH_DIRECTORY}
cp ~/data/ENA_SRP044740/tidydata/*rds ${SCRATCH_DIRECTORY}
cp ${SLURM_SUBMIT_DIR}/models/test_set/* ${SCRATCH_DIRECTORY}

# In this script the seed is not necessary but I keep it because of consistency issues
export SEED=1

test_sample=${sample}
train_sample=ALL
echo "Train_sample = All but ${test_sample}"
echo "Test_sample=${test_sample}"

# Call script
Rscript --vanilla --quiet --slave trainTestModels.R ${model} ${SEED} ${train_sample} ${test_sample}

# After the job is done we copy our output back to our directory
mkdir -p ${SLURM_SUBMIT_DIR}/results/train_${train_sample}_test_${test_sample}_${SLURM_JOBID}/ && cp ${SCRATCH_DIRECTORY}/*.RDS ${SLURM_SUBMIT_DIR}/results/train_${train_sample}_test_${test_sample}_${SLURM_JOBID}/
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}

exit 0