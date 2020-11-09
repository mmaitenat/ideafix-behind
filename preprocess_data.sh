#!/bin/bash
# Change the following variables according to your system
DIR="~/data/ENA_SRP044740/"
REF="~/genomes/hg19/ucsc.hg19.fasta"
picard_cmd="java -jar /opt/picard/build/libs/picard.jar"
gatk3_cmd="java -jar /opt/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"
gatk4="/opt/gatk-4.0.8.1/gatk"
vcfpolyx_cmd="java -jar /opt/jvarkit/dist/vcfpolyx.jar"
snpsift="/opt/snpEff/SnpSift.jar"
dbsnp="~/variants/common_all_20180423.vcf.gz"
pair_filename="FFPE-FF_pairs.csv"
run_accessions="run_accessions.txt"

cat $run_accessions | while read sample; do
    echo "Starting to work on sample ${sample}"
    cd ${DIR}${sample}
    R1=${sample}"_1.fastq.gz"
    R2=$(echo "${R1/_1/_2}")

    # Trim adapter sequences
    echo "Trimming adapters on sample ${sample}"
    trimmomatic PE -threads 6 $R1 $R2 $sample"_trimmed_R1.fastq" $sample"_unpaired_R1.fastq" $sample"_trimmed_R2.fastq" $sample"_unpaired_R2.fastq" ILLUMINACLIP:/opt/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10:10:TRUE SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:36
    mkdir aligned

    # Align data
    echo "Aligning data on sample ${sample}"
    bwa mem -t 6 $REF ${sample}"_trimmed_R1.fastq" ${sample}"_trimmed_R2.fastq" | samtools sort -O BAM -o aligned/${sample}.bam -

    # Add read groups
    echo "Adding read-groups on sample ${sample}"
    $picard_cmd AddOrReplaceReadGroups \
        I=aligned/${sample}".bam" \
        O=aligned/${sample}"_rg.bam" \
        RGLB=lib1 \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM=${sample}
    mv aligned/${sample}_rg.bam aligned/${sample}.bam

    # Index
    samtools index aligned/${sample}.bam

    # Single-mode variant calling
    echo "Calling variants on sample ${sample}"
    $gatk4 Mutect2 \
        -R $REF \
        -I aligned/${sample}.bam \
        -tumor ${sample} \
        -O aligned/${sample}.vcf \
        --annotation StrandBiasBySample

    # Mark TP and FPs with FilterMutectCalls. Calls that are likely true positives get the PASS label in the FILTER field, and calls that are likely false positives are labeled with the reason(s) for filtering in the FILTER field of the VCF
    echo "Filtering variants on sample ${sample}"
    $gatk4 FilterMutectCalls \
        -V aligned/${sample}.vcf \
        -O aligned/${sample}_filtermarks.vcf

    # Annotate SNPs
    echo "Annotating SNPs on sample ${sample}"
    java -jar $snpsift annotate $dbsnp aligned/${sample}_filtermarks.vcf >aligned/${sample}_filtermarks_annotated.vcf

    # Annotate homopolymers
    ## VCFPolyX
    $vcfpolyx_cmd -R ${REF} \
        --skip-filtered \
        -o aligned/${sample}_VCFPolyx.vcf aligned/${sample}_filtermarks_annotated.vcf

    ## VariantAnnotator - HomopolymerRun
    # The homopolymer length can only me computed for biallelic sites. Hence, the output vcf will be shorter than the input one
    $gatk3_cmd -T VariantAnnotator \
        -R ${REF} \
        -V aligned/${sample}_filtermarks_annotated.vcf \
        -o aligned/${sample}_HRun.vcf \
        -A HomopolymerRun \
        --reference_window_stop 200

    # Paired variant calling only for FFPE samples
    if cut -d"," -f1 $pair_filename | grep -q ${sample}; then
        echo "Calling deaminations on sample ${sample}"
        FF_sample=$(grep ${sample} ${pair_filename} | cut -d"," -f2)
        FF_bam_filename="${DIR}${FF_sample}/aligned/${FF_sample}.bam" #FF samples are preprocessed before FFPE
        $gatk4 Mutect2 \
            -R $REF \
            -I aligned/${sample}.bam \
            -tumor ${sample} \
            -I ${FF_bam_filename} \
            -normal ${FF_sample} \
            -O aligned/${sample}_FFPE_vs_FF.vcf \
            --annotation StrandBiasBySample

        # Mark TP and FPs with FilterMutectCalls. Calls that are likely true positives get the PASS label in the FILTER field, and calls that are likely false positives are labeled with the reason(s) for filtering in the FILTER field of the VCF
        echo "Filtering deaminations on sample ${sample}"
        $gatk4 FilterMutectCalls \
            -V aligned/${sample}_FFPE_vs_FF.vcf \
            -O aligned/${sample}_FFPE_vs_FF_filtermarks.vcf

        # Annotate SNPs
        echo "Annotating deaminations on sample ${sample}"
        java -jar $snpsift annotate $dbsnp aligned/${sample}_FFPE_vs_FF_filtermarks.vcf >aligned/${sample}_FFPE_vs_FF_filtermarks_annotated.vcf

        # Annotate homopolymers
        ## VCFPolyX
        $vcfpolyx_cmd -R ${REF} \
            --skip-filtered \
            -o aligned/${sample}_FFPE_vs_FF_VCFPolyx.vcf aligned/${sample}_FFPE_vs_FF_filtermarks_annotated.vcf

        ## VariantAnnotator - HomopolymerRun
        # The homopolymer length can only me computed for biallelic sites. Hence, the output vcf will be shorter than the input one
        $gatk3_cmd -T VariantAnnotator \
            -R ${REF} \
            -V aligned/${sample}_FFPE_vs_FF_filtermarks_annotated.vcf \
            -o aligned/${sample}_FFPE_vs_FF_HRun.vcf \
            -A HomopolymerRun \
            --reference_window_stop 200
    fi
    echo "Finished with sample ${sample}"
done
