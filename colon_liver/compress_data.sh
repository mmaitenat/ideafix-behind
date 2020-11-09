#!/bin/bash

# This is for bcftools isec to works

# Change the following variables according to your system
DIR="~/data/colon_liver/"
run_accessions=run_accessions.txt

cat $run_accessions | while read sample; do
  bcftools view -Oz -o "${DIR}${sample}/aligned/${sample}_filtermarks_annotated.vcf.gz" "${DIR}${sample}/aligned/${sample}_filtermarks_annotated.vcf"
  htsfile "${DIR}${sample}/aligned/${sample}_filtermarks_annotated.vcf.gz"
  bcftools index "${DIR}${sample}/aligned/${sample}_filtermarks_annotated.vcf.gz"
done

exit 0
