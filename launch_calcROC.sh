#!/bin/bash

## This script launches all the combinations of parameters for a job
## You set the possible values from the command line
## Use the script as:
## 		bash launch_calcROC.sh --features='frag_length allele_freq'

### PARSING OF THE ARGUMENTS ###

for i in "$@"
do
    case $i in
        -f=*|--features=*)
            _FEATURES="${i#*=}"
            shift
        ;;
        -h|--help)
            echo 'NAME'
            echo '  launch_calcROC.sh - launch multiple jobs in the queue'
            echo
            echo 'SYNOPSIS'
            echo '  launch_calcROC.sh [OPTIONS...]'
            echo
            echo 'DESCRIPTION'
            echo '  The options collect the list of possible values and then all the combinations are put in the queue'
            echo
            echo 'OPTIONS:'
            echo '  -f, --features'
            echo '          List of features to calculate power of. Example (remember single quotes around the values) -f=frag_length allele_freq'
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

#################################
## Launch all the combinations ##
#################################


for feature in $_FEATURES
do
  sbatch --export=feature=$feature run_calcROC.sbatch
done
