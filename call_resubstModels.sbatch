#!/bin/bash
#SBATCH --job-name=resubst
#SBATCH --output=resubst-models_%A.out
#SBATCH --error=resubst-models_%A.err
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
#SBATCH --time=1000:00:00
#SBATCH --partition=2018allq


# Define and create a unique scratch directory for this job

pwd; hostname; date
echo ${model}
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
seed=1

cp ${SLURM_SUBMIT_DIR}/resubstModels.R ${SCRATCH_DIRECTORY}
cp ~/data/ENA_SRP044740/tidydata/*rds ${SCRATCH_DIRECTORY}
cp ${SLURM_SUBMIT_DIR}/models/* ${SCRATCH_DIRECTORY}


Rscript --vanilla --quiet --slave resubstModels.R ${model} ${seed} ${sample}

# After the job is done we copy our output back to our directory
mkdir -p ${SLURM_SUBMIT_DIR}/results/${sample}_${SLURM_JOBID}/ && cp ${SCRATCH_DIRECTORY}/*.RDS ${SLURM_SUBMIT_DIR}/results/${sample}_${SLURM_JOBID}/
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}

exit 0