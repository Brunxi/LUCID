# LUCID: Analysis Pipeline for Fungal Pathogens

LUCID is a comprehensive bioinformatics analysis pipeline capable of detecting optimal targets for RNA interference-based control strategies in any phytopathogenic fungi, regardless of their lifestyle. By integrating transcriptomic data and comparative genomics, LUCID identifies genes that are:
1. Highly expressed (obligate biotrophs) or differentially expressed during infection (non-obligate biotrophs)
2. Conserved across a selected group of phytopathogenic fungi
3. Categorized into two key groups:
   - **Conserved Essential Proteins (CEPs)**: Conserved, highly expressed proteins with known functions essential for pathogenicity
   - **Conserved Non-Annotated Proteins (CNAPs)**: Novel conserved proteins with high expression but lacking functional annotation

This dual-approach pipeline provides valuable targets for developing effective and specific RNA interference-based control methods against fungal plant pathogens.

---

## Table of Contents
1. [Introduction](#introduction)
2. [Pipeline Overview](#pipeline-overview)
3. [Installation and Dependencies](#installation-and-dependencies)
4. [Data Preparation](#data-preparation)
5. [Orthogroup Analysis](#orthogroup-analysis)
6. [RNA-seq Processing](#rna-seq-processing)
7. [Target Identification](#target-identification)
8. [Results and Output Files](#results-and-output-files)
9. [Troubleshooting](#troubleshooting)

---

## Introduction

LUCID (Lifestyle-Unrestricted Conserved Interference-Directed targets) is a bioinformatics pipeline designed to identify the most promising targets for RNA interference-based control strategies in phytopathogenic fungi. 

The pipeline's key strength is its ability to analyze any phytopathogenic fungus regardless of its lifestyle (obligate biotrophs, hemibiotrophs, or necrotrophs). LUCID integrates transcriptomics data with comparative genomics to identify genes that are:

1. **Highly expressed or infection-specific**: The pipeline adaptively uses either TPM-based expression analysis (for obligate biotrophs) or differential expression analysis (for non-obligate biotrophs)

2. **Evolutionary conserved**: Using OrthoFinder to identify orthologs conserved across user-selected groups of phytopathogenic fungi

3. **Functionally important**: Through homology searches against a database of known pathogenicity factors

The final output categorizes potential targets into two groups:
- **Conserved Essential Proteins (CEPs)**: Proteins that are conserved, highly expressed, and have homology to known essential pathogenicity factors
- **Conserved Non-Annotated Proteins (CNAPs)**: Conserved, highly expressed proteins without current functional annotation that represent novel potential targets

This comprehensive approach makes LUCID a powerful tool for developing targeted, effective RNA interference strategies against diverse fungal plant pathogens.

---

## Pipeline Overview

```
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│  Installation   │────►│ Data Preparation │────►│   Orthogroup    │
│  & Dependencies │     │   (Proteomes,    │     │    Analysis     │
└─────────────────┘     │  Genomes, Reads) │     │  (OrthoFinder)  │
                        └─────────────────┘     └────────┬────────┘
                                                          │
                                                          ▼
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│  Results and    │◄────│     Target      │◄────│    RNA-seq      │
│  Output Files   │     │  Identification  │     │   Processing    │
└─────────────────┘     └─────────────────┘     └─────────────────┘
```

---

## Installation and Dependencies

### Conda Environment Setup

Create a conda environment with the required dependencies:

```bash
conda create -n lucid-env
conda activate lucid-env

# Install bioinformatics tools
conda install -c bioconda diamond
conda install -c bioconda orthofinder

# Install R and required packages
conda install -c conda-forge r-base
```

### R Package Requirements

Once R is installed, the following libraries need to be available:

```R
# Install these packages in R
install.packages(c("dplyr", "data.table", "lattice", "ggplot2", "tidyr", "DT", "httr", "jsonlite"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SummarizedExperiment", "Rsubread", "GenomicRanges", "Biostrings"))
```

Required R libraries:
- dplyr - For data manipulation
- data.table - For fast data import and manipulation
- SummarizedExperiment - For handling genomic data
- lattice - For visualization
- Rsubread - For RNA-seq alignment and counting
- ggplot2 - For visualization
- GenomicRanges - For genomic interval manipulation
- tidyr - For data reshaping
- DT - For interactive tables
- httr - For HTTP requests
- Biostrings - For biological sequence manipulation
- jsonlite - For JSON parsing

---

## Data Preparation

### Directory Structure

Create the following directory structure for your analysis:

```
LUCID/
├── data/
│   ├── genome/       # Genome files go here
│   ├── proteomes/    # Proteome files go here
│   ├── reads/        # RNA-seq reads go here
│   └── references/   # Contains OrthoFinder output and essential protein database
└── script/           # Analysis scripts
```

**Important Note:** The `references/` directory contains a database of essential proteins for pathogenesis that is used as a reference for detecting Conserved Essential Proteins (CEPs) through comparative genomics using DIAMOND.

### Genome Files

**Location**: `LUCID/data/genome/`

**Required files**:
- `genome.fa` - Genome sequence in FASTA format
- `genome.gtf` - Genome annotation in GTF format

**Important**: Files must be named exactly as specified above.

### Proteome Files

**Location**: `LUCID/data/proteomes/`

**Required format**: 
- FASTA files containing protein sequences
- **Must be downloaded from ENSEMBL FUNGI**

### RNA-seq Data

**Location**: `LUCID/data/reads/`

**Required naming convention**:
- Files must end with `*_1.fastq` and `*_2.fastq` for paired-end reads
- The prefix must match the identifiers used in `coldata.csv` (if applicable)

**Options for obtaining data**:
1. Use your own FASTQ files (place in the reads directory)
2. Download from NCBI SRA using the provided script

---

## Orthogroup Analysis

OrthoFinder is used to identify orthologous groups across multiple fungal proteomes.

### Running OrthoFinder

Execute the OrthoFinder analysis script:

```bash
bash LUCID/script/run_orthofinder.sh
```

**Input**: Proteome files in `LUCID/data/proteomes/`
**Output**: Orthogroup analysis in `LUCID/data/references/`

You can modify input/output paths within the script if needed.

---

## RNA-seq Processing

### Downloading SRA Data (Optional)

To download RNA-seq data from NCBI's Sequence Read Archive:

```bash
bash LUCID/script/downloadSRA.sh
```

Edit this script to include specific SRA accession numbers of interest.

### Sample Metadata for Differential Expression

For differential expression analysis (non-obligate biotrophs), create a `coldata.csv` file in `LUCID/data/genome/` with this format:

```csv
illumina_code,type,treatment
SRR6924534,paired-end,inf
SRR6924535,paired-end,inf
SRR6924536,paired-end,inf
SRR6924547,paired-end,ctrl
SRR6924548,paired-end,ctrl
SRR6924549,paired-end,ctrl
```

Where:
- `illumina_code`: SRA accession number
- `type`: Sequencing type (paired-end or single-end)
- `treatment`: Condition (inf = infection, ctrl = control)

---

## Results and Output Files

The analysis pipeline generates several important output files:

### Common Output Files

- `CEPs.fasta`: Sequences of identified Core Effector Proteins
- `CNAPs.fasta`: Sequences of identified Conserved Non-Annotated Proteins
- `CEPs_ids.txt`: List of CEP protein identifiers
- `CNAPs_ids.txt`: List of CNAP protein identifiers
- `CEPs_uniprot_annotations.tsv`: UniProt annotation information for CEPs
- `CNAPs_uniprot_annotations.tsv`: UniProt annotation information for CNAPs
- `analysis_summary.tsv`: Summary statistics of the analysis
- `diamond_results.txt`: Results from DIAMOND homology search

### For Obligate Biotrophs

- `tpm_matrix.tsv`: Gene expression data in TPM format
- `highly_expressed_genes.tsv`: List of genes above the TPM threshold

### For Non-Obligate Biotrophs

- `count_matrix.tsv`: Raw read count matrix
- `significant_genes_inf_vs_ctrl.tsv`: Genes significantly upregulated during infection
- `pca_plot.png`: Principal Component Analysis visualization of samples

These files provide valuable information for designing RNA interference experiments targeting conserved, highly expressed genes in phytopathogenic fungi.

---

## How LUCID Works

LUCID integrates several bioinformatics approaches to identify key pathogenicity-related proteins that could serve as targets for RNA interference-based control strategies:

1. **Ortholog Identification**: Using OrthoFinder to identify conserved genes across multiple fungal species.

2. **Expression Analysis**: Two alternative approaches based on fungal lifestyle:
   - **For obligate biotrophs**: TPM-based expression quantification (since these fungi cannot be cultured in vitro)
   - **For non-obligate biotrophs**: Differential expression analysis between infection and control conditions using DESeq2

3. **Target Selection**: Integrates orthology, expression, and homology data to identify:
   - **Core Effector Proteins (CEPs)**: Conserved proteins essential for pathogenicity
   - **Conserved Non-Annotated Proteins (CNAPs)**: Novel conserved proteins that lack functional annotation

4. **Homology Search**: Uses DIAMOND to compare proteins against a database of known essential pathogenicity proteins.

5. **Annotation**: Maps identified proteins to UniProt for functional annotation and provides comprehensive output files.

This integrated approach helps identify promising RNAi targets that are:
- Conserved across multiple pathogenic fungi (orthology analysis)
- Highly expressed or upregulated during infection (expression analysis)
- Functionally related to pathogenicity (homology to essential proteins)
- Potentially novel virulence factors (identification of CNAPs)

---

## Troubleshooting

### Common Issues

- **Error with OrthoFinder**: Ensure all proteome files are valid FASTA format and properly downloaded from ENSEMBL FUNGI
- **RNA-seq analysis error**: Verify that genome files are correctly named as `genome.fa` and `genome.gtf`
- **SRA download failure**: Check internet connection and SRA accession numbers

### Getting Help

For additional assistance, please refer to the manuscript or contact the authors.
