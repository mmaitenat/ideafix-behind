#!/bin/bash

# Change the following variables according to your system
DIR="~/data/ENA_SRP044740/"
REF="~/genomes/hg19/ucsc.hg19.fasta"
gatk_cmd="/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -jar ~/gatk-4.0.8.1/gatk-package-4.0.8.1-local.jar"
run_accessions=run_accessions.txt

cat $run_accessions | while read sample; do
    cd ${DIR}/${sample}/aligned/
    BAMFILE="${sample}.bam"

    ${gatk_cmd} CollectF1R2Counts -R $REF -I $BAMFILE \
        -alt-table tumor-alt.tsv \
        -ref-hist tumor-ref.metrics \
        -alt-hist tumor-alt.metrics

    ${gatk_cmd} LearnReadOrientationModel \
        -alt-table tumor-alt.tsv \
        -ref-hist tumor-ref.metrics \
        -alt-hist tumor-alt.metrics \
        -O artifact-prior.tsv

    ${gatk_cmd} Mutect2 -R $REF -I $BAMFILE \
        -tumor ${sample} \
        -O "${sample}_gatk_FFPE_annot.vcf" \
        --annotation StrandBiasBySample \
        --orientation-bias-artifact-priors artifact-prior.tsv

    ${gatk_cmd} FilterMutectCalls \
        -V "${sample}_gatk_FFPE_annot.vcf" \
        -O "${sample}_gatk_FFPE_annot_filtermarks.vcf"
done
