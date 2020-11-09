#!/bin/bash

# Change the following variables according to your system
DIR="~/data/colon_liver/"
run_accessions=run_accessions.txt

## Depth of coverage

# Change the following variables according to your system
BEDFILE="~/exomes/truseq-exome-targeted-regions-manifest-v1-2.bed"
samples_str=$(cat $run_accessions | tr '\n' ' ')
TEMPDIR=temp
currentdir=$(pwd)

mkdir -p ${TEMPDIR}
cd ${TEMPDIR}

# Create symlinks for easyness in giving the bam files to samtools depth
cat $run_accessions | while read sample; do
    ln -s $DIR/${sample}/aligned/${sample}.bam $sample
done

samtools depth -b $BEDFILE -m 0 -Q 10 -q 10 -a ${samples_str} >colon_liver.depth
mv colon_liver.depth ${DIR}
cd $currentdir

## Fragment length

cat $run_accessions | while read sample; do
    samtools view ${DIR}/${sample}/aligned/${sample}.bam | awk '{print $9}' >${DIR}/${sample}/aligned/${sample}_fragment_length.txt
    echo "Finished sample ${sample}"
done

## Duplication levels

# Change the following variables according to your system
fastqc="~/FastQC/fastqc"

cat $run_accessions | while read sample; do
    echo "Analyzing quality of sample ${sample}"
    cd ${DIR}${sample}
    fastq_files=($(find . -type f -name "EGAR*.fastq.gz" | awk -F'/' '{print $2}'))
    R1=${fastq_files[0]}
    R2=${fastq_files[1]}
    ${fastqc} ${R1} ${R2}
done

# Duplication levels are manually extracted from the generated html files

## Damage levels (GIV score) by Damage-estimator

# Change the following variables according to your system

TOOLDIR="~/Damage-estimator-master/"
TEMPDIR="~/damage_temp/"
REF="~/genomes/hg19/ucsc.hg19.fasta"

cat $run_accessions | while read sample; do
    echo "Analizing sample ${sample}"
    mkdir -p ${TEMPDIR}${sample}
    perl ${TOOLDIR}/split_mapped_reads.pl --bam ${DIR}/${sample}/aligned/${sample}.bam -genome $REF -mpileup1 ${TEMPDIR}${sample}/file1.mpileup -mpileup2 ${TEMPDIR}${sample}/file2.mpileup -Q 10 -q 20
    perl ${TOOLDIR}/estimate_damage.pl --mpileup1 ${TEMPDIR}${sample}/file1.mpileup --mpileup2 ${TEMPDIR}${sample}/file2.mpileup --id ${sample} --qualityscore 20 > ${DIR}${sample}/damage.out
    rm ${TEMPDIR}${sample}/file1.mpileup ${TEMPDIR}${sample}/file2.mpileup
    echo "Finished sample ${sample}"
done

exit 0