#!/bin/bash

#### Breast cancer data

# Preprocess data
./preprocess_data.sh

# Obtain non-deaminations by intersecting matched FF and FFPE samples
./compress_data.sh
./intersect_data.sh

# Prepare variant data matrices and labels with our method
./create_data_matrices.sh

# Univariate feature analysis
./launch_calcROC.sh --features='frag_length allele_freq alt_bases norm_alt_bases ref_bases norm_ref_bases base_qual base_qual_frac pos_from_end map_qual FdeamC SOB SBGuo SBGATK norm_pos_from_end hp_length'
/opt/R-3.5.1/bin/Rscript --vanilla --quiet --slave calcJMI.R

# Label variants by LearnReadOrientationModel
./run_LROM.sh

# Estimate damage
./estimate_damage.sh

# Get damage numbers and make plots
/opt/R-3.5.1/bin/Rscript --vanilla --quiet --slave make_damage_figs.R

# Train models and estimate their performance
# Code for running in SLURM environment
./launch_calcError --samples='sample1 sample2 sample8 sample10 sample7 sample12 sample13 sample11 sample3 sample5 sample4 sample6 sample9 sample6_2 sample1_2 sample11_2 sample12_2 sample10_2 sample5_2 sample8_2 sample1_3 sample2_2 sample10_3 sample7_2 sample13_2 sample4_2 sample3_2' --model='NB NN RF xgbTree logReg' --error='resubstitution kcv testset'

# Make data plots and performance plots
/opt/R-3.5.1/bin/Rscript --vanilla --quiet --slave make-figs.R

# Train final models
/opt/R-3.5.1/bin/Rscript --vanilla --quiet --slave train_final_models.R


#### Colon and liver cancer data

cd colon_and_liver

# Preprocess data
./preprocess_data.sh

# Obtain non-deaminations by intersecting matched FF and FFPE samples
./compress_data.sh
./intersect_data.sh

# Prepare variant data matrices and labels with our method
./create_data_matrices.sh

# Estimate damage
./estimate_damage.sh

# Get damage numbers and make damage and data plots
/opt/R-3.5.1/bin/Rscript --vanilla --quiet --slave make_damage_figs.R

# Test the final models, run other state-of-the-art methods and plot ROC curves
cd ..
/opt/R-3.5.1/bin/Rscript --vanilla --quiet --slave compare_with_other_tools.R