#!/bin/bash

## This script launches jobs for calculating model performance
## You set the possible values from the command line
## Use the script as:
## 		bash launch_calcError.sh --samples='sample1 sample2' --model='NB NN RF xgbTree logReg' --error='resubstitution kcv testset'

### PARSING OF THE ARGUMENTS ###

for i in "$@"; do
    case $i in
    -m=* | --model=*)
        _MODELS="${i#*=}"
        shift
        ;;
    -s=* | --samples=*)
        _SAMPLES="${i#*=}"
        shift
        ;;
    -e=* | --error=*)
        _ERROR="${i#*=}"
        shift
        ;;
    -h | --help)
        echo 'NAME'
        echo '  launch_calcError.sh - launch multiple jobs in the queue'
        echo
        echo 'SYNOPSIS'
        echo '  launch_calcError.sh [OPTIONS...]'
        echo
        echo 'DESCRIPTION'
        echo '  The options collect the list of possible values and then all the combinations are put in the queue'
        echo
        echo 'OPTIONS:'
        echo '  -m, --model'
        echo '          List of learning algorithms. Valid values are NB, NN, RF, xgbTree and logReg. Example (remember single quotes around the values) --model=NN NB'
        echo '  -s, --samples'
        echo '          List of samplenames of data to be used as a test-set. Example (remember single quotes around the values) -s=sample1 sample3 sample9'
        echo '  -e, --error'
        echo '          Error-estimation to be computed. Valid values as kcv, resubstitution and testset. Example (remember single quotes around the values) -e=resubstitution'
        exit 0
        shift
        ;;
    *)
        echo 'The option you passed is not known. Use -h or --help to get the valid arguments'
        exit 1
        ;;
    esac
done

######## END PARSING ARGUMENTS ########

#################################################
## Determine script to be run and its options ###
#################################################

if [ $_ERROR == "resubstitution" ]; then
    scriptname=call_resubstModels.sbatch
elif [ $_ERROR == "kcv" ]; then
    scriptname=call_trainModels.sbatch
elif [ $_ERROR == "testset" ]; then
    scriptname=call_trainTestModels.sbatch
else
    echo 'The error option you passed is not known. Use one from kcv, resubstitution or testset.'
    exit 1
fi

#################################
## Launch all the combinations ##
#################################

# We can pass variables to the slurm script using the --export option
for model in $_MODELS; do
    for sample in $_SAMPLES; do
        if [ $_ERROR == "testset" ]; then
            # note that this call can be set into one with resubst-error, as both scripts share parameters
            sbatch --export=model=$model,sample=${sample} ${scriptname}
        elif [[ $_ERROR == "kcv" ]]; then
            for task in {1..10}; do
                sbatch --export=model=$model,sample=${sample},task=$task ${scriptname}
            done
        else # resubstitution error
            sbatch --export=model=$model,sample=${sample} ${scriptname}
        fi
    done
done
