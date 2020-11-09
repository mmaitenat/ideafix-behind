#!/bin/bash
DIR="/media/NGS3/Bonnet_data/data/"
REF=/var/tmp/mtellaetxe/Data/RefGenomes/hg19/ucsc.hg19.fasta
picard_cmd="java -jar /home/mtellaetxe/picard.jar"
gatk3_cmd="/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -jar /opt/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"
gatk4="/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -jar /home/mtellaetxe/gatk-4.0.8.1/gatk-package-4.0.8.1-local.jar"
vcfpolyx_cmd="java -jar /opt/jvarkit/vcfpolyx.jar"
snpsift="/opt/snpEff/SnpSift.jar"
dbsnp="/var/tmp/mtellaetxe/Data/variants/common_all_20180423.vcf.gz"
pairs="/home/mtellaetxe/FFPE/data/sample_datasets/Bonnet/maxwell_FFPE_FF_tumor_pairs.csv"

while read pair; do
  tumor=`echo $pair | awk -F, '{print $1}'`
  normal=`echo $pair | awk -F, '{print $2}'`
  tumor_bamfile="${DIR}${tumor}/aligned/${tumor}.bam"
  normal_bamfile="${DIR}${normal}/aligned/${normal}.bam"
  echo "Tumor filename is ${tumor_bamfile} and normal filename is ${normal_bamfile}"

  ${gatk4} Mutect2 \
  -R $REF \
  -I ${tumor_bamfile} \
  -tumor ${tumor} \
  -I ${normal_bamfile} \
  -normal ${normal} \
  -O "${DIR}${tumor}/aligned/${tumor}_vs_${normal}.vcf" \
  --annotation StrandBiasBySample

  # Mark TP and FPs with FilterMutectCalls. Calls that are likely true positives get the PASS label in the FILTER field, and calls that are likely false positives are labeled with the reason(s) for filtering in the FILTER field of the VCF
  echo "Filtering deaminations on sample ${tumor}"
  $gatk4 FilterMutectCalls \
  -V "${DIR}${tumor}/aligned/${tumor}_vs_${normal}.vcf" \
  -O "${DIR}${tumor}/aligned/${tumor}_vs_${normal}_filtermarks.vcf"

  # Annotate SNPs
  echo "Annotating SNPs on sample ${sample}"
  java -jar $snpsift annotate $dbsnp "${DIR}${tumor}/aligned/${tumor}_vs_${normal}_filtermarks.vcf" > "${DIR}${tumor}/aligned/${tumor}_vs_${normal}_filtermarks_annotated.vcf"

  # Annotate homopolymers
  ## VCFPolyX
  ### el skip-filtered hace que las variantes que no son PASS aparezcan en el vcf pero no calcula el POLYX para ellas
  $vcfpolyx_cmd -R ${REF} \
  --skip-filtered \
  -o "${DIR}${tumor}/aligned/${tumor}_vs_${normal}_VCFPolyx.vcf" "${DIR}${tumor}/aligned/${tumor}_vs_${normal}_filtermarks_annotated.vcf"

  ## VariantAnnotator - HomopolymerRun
  # The homopolymer length can only me computed for biallelic sites. Hence, the output vcf will be shorter than the input one
  $gatk3_cmd -T VariantAnnotator \
  -R ${REF} \
  -V "${DIR}${tumor}/aligned/${tumor}_vs_${normal}_filtermarks_annotated.vcf" \
  -o "${DIR}${tumor}/aligned/${tumor}_vs_${normal}_HRun.vcf" \
  -A HomopolymerRun \
  --reference_window_stop 200

  echo "Finished with sample ${tumor}"
done <$pairs

exit 0
