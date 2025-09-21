#!/bin/bash
# -------------------------------
# Bulk RNA-seq Pipeline (Paired-End, Multiple Samples, gzipped FASTQ)
# Author: Maroof
# Date: $(date)
# -------------------------------

set -euo pipefail
SECONDS=0

# -------------------------------
# User-defined variables
# -------------------------------
FASTQ_DIR="data"

# Output directories
OUTDIR="results"
QC_DIR="$OUTDIR/qc"
TRIM_DIR="$OUTDIR/trimmed"
ALIGN_DIR="$OUTDIR/alignment"
COUNT_DIR="$OUTDIR/counts"

# Reference files
GENOME_INDEX="HISAT2/grch38/genome"
GTF="../hg38/Homo_sapiens.GRCh38.106.gtf"

# Tools
TRIMMOMATIC="~/Desktop/demo/tools/Trimmomatic-0.39/trimmomatic-0.39.jar"

# Threads
THREADS=4

# -------------------------------
# Create output directories
# -------------------------------
mkdir -p $QC_DIR $TRIM_DIR $ALIGN_DIR $COUNT_DIR

# -------------------------------
# Loop through paired FASTQ files
# Assumes naming: sample_R1.fastq.gz and sample_R2.fastq.gz
# -------------------------------
for R1 in $FASTQ_DIR/*_R1.fastq.gz; do
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    R2="$FASTQ_DIR/${SAMPLE}_R2.fastq.gz"

    echo "==============================="
    echo "Processing sample: $SAMPLE"
    echo "==============================="

    # STEP 1: FastQC (before trimming)
    fastqc $R1 $R2 -o $QC_DIR

    # STEP 2: Trimming with Trimmomatic (paired-end mode)
    TRIM_R1_PAIRED="$TRIM_DIR/${SAMPLE}_R1_trimmed_paired.fastq.gz"
    TRIM_R1_UNPAIRED="$TRIM_DIR/${SAMPLE}_R1_trimmed_unpaired.fastq.gz"
    TRIM_R2_PAIRED="$TRIM_DIR/${SAMPLE}_R2_trimmed_paired.fastq.gz"
    TRIM_R2_UNPAIRED="$TRIM_DIR/${SAMPLE}_R2_trimmed_unpaired.fastq.gz"

    java -jar $TRIMMOMATIC PE -threads $THREADS \
        $R1 $R2 \
        $TRIM_R1_PAIRED $TRIM_R1_UNPAIRED \
        $TRIM_R2_PAIRED $TRIM_R2_UNPAIRED \
        TRAILING:10 -phred33
    echo "Trimmomatic finished for $SAMPLE"

    # QC after trimming (paired reads only)
    fastqc $TRIM_R1_PAIRED $TRIM_R2_PAIRED -o $QC_DIR

    # STEP 3: Alignment with HISAT2
    BAM="$ALIGN_DIR/${SAMPLE}.bam"
    hisat2 -q --rna-strandness RF -x $GENOME_INDEX \
        -1 $TRIM_R1_PAIRED -2 $TRIM_R2_PAIRED -p $THREADS \
        | samtools sort -@ $THREADS -o $BAM
    samtools index $BAM
    echo "HISAT2 finished for $SAMPLE"

    # STEP 4: Quantification with featureCounts
    featureCounts -T $THREADS -p -B -C -a $GTF \
        -o $COUNT_DIR/${SAMPLE}_featurecounts.txt \
        $BAM
    echo "featureCounts finished for $SAMPLE"

done

# -------------------------------
# Completion message
# -------------------------------
duration=$SECONDS
echo "All paired-end samples processed in $(($duration / 60)) minutes and $(($duration % 60)) seconds."
