#!/bin/bash
DIR="~/data/ENA_SRP044740/"
pair_filename="FFPE-FF_pairs.csv"
FFPE_samples=($(awk -F "," '{print $1}' $pair_filename))
FF_samples=($(awk -F "," '{print $2}' $pair_filename))
len=${#FFPE_samples[@]}

for ((i = 0; i < $len; i++)); do
  bcftools isec -f PASS -p "${DIR}${FFPE_samples[$i]}/aligned/" -o "${FFPE_samples[$i]}_FFPE_FF_intersection.vcf.gz" "${DIR}${FFPE_samples[$i]}/aligned/${FFPE_samples[$i]}_filtermarks_annotated.vcf.gz" "${DIR}${FF_samples[$i]}/aligned/${FF_samples[$i]}_filtermarks_annotated.vcf.gz"
done

# Note that the output is not the "${FFPE_samples[$i]}_FFPE_FF_intersection.vcf.gz" file I was expecting but 4 vcf files: 0000.vcf, 0001.vcf, 0002.vcf and 0003.vcf. Intersection is 0002.vcf