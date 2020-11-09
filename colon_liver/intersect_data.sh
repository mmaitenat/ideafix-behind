#!/bin/bash

# Change the following variables according to your system
pair_filename="maxwell_FFPE_FF_tumor_pairs.csv"
DIR="~/data/colon_liver/"

while read pair; do
  FFPE=`echo $pair | awk -F, '{print $1}'`
  FF=`echo $pair | awk -F, '{print $2}'`
  echo "FFPE filename is $FFPE and FF filename is $FF"
  bcftools isec -f PASS -p "${DIR}${FFPE}/aligned/" -w1 "${DIR}${FFPE}/aligned/${FFPE}_filtermarks_annotated.vcf.gz" "${DIR}${FF}/aligned/${FF}_filtermarks_annotated.vcf.gz" && mv "${DIR}${FFPE}/aligned/0002.vcf" "${DIR}${FFPE}/aligned/${FFPE}_${FF}_isec.vcf"
done <${pair_filename}

exit 0
