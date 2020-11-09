#!/bin/bash

# Change the following variables according to your system
DIR="~/data/ENA_SRP044740/"
OUTDIR="~/data/ENA_SRP044740/tidydata"
pair_filename="sample_pairs_with_names.csv"

FFPE_samples=($(awk -F "," '{print $2}' $pair_filename))
samplenames=($(awk -F "," '{print $3}' $pair_filename))
len=${#FFPE_samples[@]}


for ((i = 0; i < $len; i++)); do
    echo "iter ${i}"
    Rscript prepare_data_matrices.R --samplename ${samplenames[$i]} --outfolder ${OUTDIR} --deam.filename "${DIR}${FFPE_samples[${i}]}/aligned/${FFPE_samples[${i}]}_FFPE_vs_FF_filtermarks_annotated.vcf" --deam.HRun.filename "${DIR}${FFPE_samples[${i}]}/aligned/${FFPE_samples[${i}]}_FFPE_vs_FF_HRun.vcf" --deam.Polyx.filename "${DIR}${FFPE_samples[${i}]}/aligned/${FFPE_samples[${i}]}_FFPE_vs_FF_VCFPolyx.vcf"  --somatic.mut.filename "${DIR}${FFPE_samples[$i]}/aligned/0002.vcf" --somatic.mut.source FFPE --somatic.mut.HRun.filename "${DIR}${FFPE_samples[$i]}/aligned/${FFPE_samples[$i]}_HRun.vcf" --somatic.mut.Polyx.filename "${DIR}${FFPE_samples[$i]}/aligned/${FFPE_samples[$i]}_VCFPolyx.vcf"
done

