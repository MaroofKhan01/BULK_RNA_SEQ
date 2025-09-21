#!/bin/bash
# -------------------------------
# Bulk RNA-seq Pipeline (Single-End, Multiple Samples, gzipped FASTQ)
# Author: Maroof Khan
# Date: $(date)
# -------------------------------

set -euo pipefail
SECONDS=0

# -------------------------------
# User-defined variables
# -------------------------------
# Input directory containing FASTQ.GZ files
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
THREADS=30

# -------------------------------
# Create output directories
# -------------------------------
mkdir -p $QC_DIR $TRIM_DIR $ALIGN_DIR $COUNT_DIR

# -------------------------------
# Loop through each FASTQ.GZ file
# -------------------------------
for FASTQ in $FASTQ_DIR/*.fastq.gz; do
    SAMPLE=$(basename "$FASTQ" .fastq.gz)
    echo "==============================="
    echo "Processing sample: $SAMPLE"
    echo "==============================="

    # STEP 1: FastQC (before trimming)
    fastqc $FASTQ -o $QC_DIR

    # STEP 2: Trimming with Trimmomatic
    TRIMMED_FASTQ="$TRIM_DIR/${SAMPLE}_trimmed.fastq.gz"
    java -jar $TRIMMOMATIC SE -threads $THREADS \
        $FASTQ $TRIMMED_FASTQ TRAILING:10 -phred33
    echo "Trimmomatic finished for $SAMPLE"

    # QC after trimming
    fastqc $TRIMMED_FASTQ -o $QC_DIR

    # STEP 3: Alignment with HISAT2
    BAM="$ALIGN_DIR/${SAMPLE}.bam"
    hisat2 -q --rna-strandness R -x $GENOME_INDEX \
        -U $TRIMMED_FASTQ -p $THREADS \
        | samtools sort -@ $THREADS -o $BAM
    samtools index $BAM
    echo "HISAT2 finished for $SAMPLE"

    # STEP 4: Quantification with featureCounts
    featureCounts -T $THREADS -S 2 \
        -a $GTF \
        -o $COUNT_DIR/${SAMPLE}_featurecounts.txt \
        $BAM
    echo "featureCounts finished for $SAMPLE"

done

# -------------------------------
# Completion message
# -------------------------------
duration=$SECONDS
echo "All samples processed in $(($duration / 60)) minutes and $(($duration % 60)) seconds."
