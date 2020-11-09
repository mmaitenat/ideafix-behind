#!/bin/bash

# Change the following variables according to your system
DIR="~/data/colon_liver/"
OUTDIR="~/data/colon_liver/tidydata"
pair_filename="maxwell_FFPE_FF_tumor_pairs.csv"

FFPE_samples=($(awk -F "," '{print $1}' $pair_filename))
FF_samples=($(awk -F "," '{print $2}' $pair_filename))
len=${#FFPE_samples[@]}

for ((i = 0; i < $len; i++)); do
    echo "iter ${i}"
    Rscript prepare_data_matrices.R --samplename ${FFPE_samples[$i]} --outfolder ${OUTDIR} --deam.filename "${DIR}${FFPE_samples[${i}]}/aligned/${FFPE_samples[${i}]}_vs_${FF_samples[${i}]}_filtermarks_annotated.vcf" --deam.HRun.filename "${DIR}${FFPE_samples[${i}]}/aligned/${FFPE_samples[${i}]}_vs_${FF_samples[${i}]}_HRun.vcf" --deam.Polyx.filename "${DIR}${FFPE_samples[${i}]}/aligned/${FFPE_samples[${i}]}_vs_${FF_samples[${i}]}_VCFPolyx.vcf"  --somatic.mut.filename "${DIR}${FFPE_samples[$i]}/aligned/${FFPE_samples[${i}]}_${FF_samples[${i}]}_isec.vcf" --somatic.mut.source FFPE --somatic.mut.HRun.filename "${DIR}${FFPE_samples[$i]}/aligned/${FFPE_samples[$i]}_HRun.vcf" --somatic.mut.Polyx.filename "${DIR}${FFPE_samples[$i]}/aligned/${FFPE_samples[$i]}_VCFPolyx.vcf"
done
