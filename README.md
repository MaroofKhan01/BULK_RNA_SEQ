# 🧬 Bulk RNA-seq Pipeline (Paired-End)

[![Linux](https://img.shields.io/badge/OS-Linux-blue?logo=linux&logoColor=white)](https://www.linux.org/)  
A simple, reproducible **bulk RNA-seq analysis pipeline** for **paired-end `.fastq.gz` data**, written in **Bash**.  

This pipeline automates **QC → Trimming → Alignment → Quantification** using widely used bioinformatics tools.  

---

## 📌 Workflow Overview
1. **Quality Control** → [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
2. **Read Trimming** → [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)  
3. **Alignment** → [HISAT2](https://daehwankimlab.github.io/hisat2/) + [SAMtools](http://www.htslib.org/)  
4. **Quantification** → [featureCounts](http://subread.sourceforge.net/)  

---

## 📂 Directory Structure
RNASeq_pipeline/
│── data/ # Input FASTQ.gz files (R1/R2)
│── HISAT2/ # HISAT2 genome index
│── hg38/ # GTF annotation file
│── results/ # Output directory
│ ├── qc/ # FastQC reports
│ ├── trimmed/ # Trimmed FASTQ.gz files
│ ├── alignment/ # BAM files
│ └── counts/ # Gene count tables
│── rnaseq_pipeline_pe.sh # Main pipeline script (paired-end)
│── README.md # Documentation


---

## ⚙️ Requirements
- **Operating System**: Linux / macOS  
- **Installed tools** in `$PATH`:  
  - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
  - [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)  
  - [HISAT2](https://daehwankimlab.github.io/hisat2/)  
  - [SAMtools](http://www.htslib.org/)  
  - [featureCounts](http://subread.sourceforge.net/)  

---

## 📥 Input Files
- **Paired-end FASTQ files** in `.fastq.gz` format inside the `data/` directory.  
- Naming convention:  
sample1_R1.fastq.gz
sample1_R2.fastq.gz
sample2_R1.fastq.gz
sample2_R2.fastq.gz



- **Reference files**:  
- HISAT2 genome index → `HISAT2/grch38/genome.*`  
- GTF annotation file → `hg38/Homo_sapiens.GRCh38.106.gtf`  

---

## 🚀 Usage
1. Make the script executable:
 ```
 chmod +x rnaseq_pipeline_pe.sh
Run the pipeline:

./rnaseq_pipeline_pe.sh
The pipeline will automatically:

Detect R1/R2 pairs in data/

Run FastQC → Trimming → Alignment → Quantification for each sample

📊 Output
FastQC reports → results/qc/

Trimmed paired FASTQ.gz files → results/trimmed/

Aligned BAM files → results/alignment/

FeatureCounts tables → results/counts/

📝 Notes
The script assumes stranded RNA-seq with RF orientation.

Change --rna-strandness RF in HISAT2 to FR if your library prep is forward-stranded.

featureCounts options:

-p → count fragments (paired-end)

-B → require proper pairs

-C → exclude chimeric reads

Increase THREADS in the script for faster performance.

Currently produces separate count files per sample.

A merge step can be added to generate a combined count matrix for downstream analysis in DESeq2 or edgeR.

