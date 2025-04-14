# LUCID: Locating Unique Candidate of Infection Determinants, a novel computational approach for target location and dsRNA design
 
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
6. [Target Identification](#target-identification)
7. [dsRNA Design](#dsrna-design)
8. [Results and Output Files](#results-and-output-files)
9. [Troubleshooting](#troubleshooting)

---

## Introduction

LUCID (Locating Unique Candidate of Infection Determinants, a novel computational approach for target location and dsRNA design) is a bioinformatics pipeline designed to identify the most promising targets for RNA interference-based control strategies in phytopathogenic fungi. 

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
└─────────┬───────┘     └─────────────────┘     └─────────────────┘
          │
          ▼
┌─────────────────┐
│    dsRNA        │
│    Design       │
└─────────────────┘
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
conda install -c bioconda bowtie

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

### Python Package Requirements

For the dsRNA design component, install the following Python packages:

```bash
pip install biopython pandas primer3-py matplotlib numpy sifi21-plus
```

---

## Data Preparation

### Directory Structure

The following directory structure is necessary for your analysis:

```
LUCID/
├── data/
│   ├── genome/       # Genome files go here
│   ├── proteomes/    # Proteome files go here
│   ├── reads/        # RNA-seq reads go here
│   └── references/   # Contains OrthoFinder output and essential protein database
├── output/           # Directory for output files
└── script/           # Analysis scripts
```

**Important Note:** The `references/` directory contains a database of essential proteins for pathogenesis that is used as a reference for detecting Conserved Essential Proteins (CEPs) through comparative genomics using DIAMOND.

### Genome Files
Genome assembly and annotation file of the fungal species under investigation for transcriptomic analysis.

**Location**: `LUCID/data/genome/`

**Required files**:
- `genome.fa` - Genome sequence in FASTA format
- `genome.gtf` - Genome annotation in GTF format

**Important**: Files must be named exactly as specified above.

### Proteome Files
Proteomes of the fungal species under investigation along with those of a group of phytopathogenic fungi, which will be used for conservation analysis with OrthoFinder.

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

### Downloading SRA Data (Optional)

To download RNA-seq data from NCBI's Sequence Read Archive:

```bash
bash LUCID/script/downloadSRA.sh
```

Edit this script to include specific SRA accession numbers of interest:

```bash
#!/bin/bash
# Output directory
output_dir="../data/reads/"
# Check if directory exists, create it if it doesn't
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
    echo "Directory $output_dir created successfully."
else
    echo "Directory $output_dir already exists."
fi
# List of accessions to download (modify this variable with your accessions)
accessions=("SRR6924534" "SRR6924535" "SRR6924536" "SRR6924547" "SRR6924548" "SRR6924549")

```

### Sample Metadata for Differential Expression

For differential expression analysis (non-obligate biotrophs), create a `coldata.csv` file in `LUCID/data/genome/` with this format:

```csv
illumina_code,type,treatment
SRRXXXXXXX,paired-end,inf
SRRXXXXXXX,paired-end,inf
SRRXXXXXXX,paired-end,inf
SRRXXXXXXX,paired-end,ctrl
SRRXXXXXXX,paired-end,ctrl
SRRXXXXXXX,paired-end,ctrl
```

Where:
- `illumina_code`: SRA accession number
- `type`: Sequencing type (paired-end or single-end)
- `treatment`: Condition (inf = infection, ctrl = control)

---

## Orthogroup Analysis

OrthoFinder is used to identify orthologous groups across multiple fungal proteomes, those previously located in `LUCID/data/proteomes/`

### Running OrthoFinder

Execute the OrthoFinder analysis script:

```bash
bash LUCID/script/run_orthofinder.sh
```

**Input**: Proteome files in `LUCID/data/proteomes/`  
**Output**: Orthogroup analysis in `LUCID/data/references/Orthogroups.tsv`

You can modify input/output paths within the script if needed.

---

## Target Identification

LUCID identifies RNAi targets using two different approaches depending on the fungal lifestyle:

### For Obligate Biotrophic Fungi

For obligate biotrophs, the pipeline uses a TPM (Transcripts Per Million) threshold approach to identify highly expressed genes, since these fungi cannot be cultured outside their host, making differential expression analysis impractical.

**Note:** The `coldata.csv` file is NOT required for this analysis method.

Execute:

```bash
Rscript /LUCID/script/lucid_biotrophs.R \
  <genome_gtf_dir> \
  <reads_dir> \
  <tpm_threshold> \
  <essential_fasta_path> \
  <orthofinder_tsv_path> \
  <output_dir> \
  <proteome_path>
```

**Parameters**:
- `genome_gtf_dir`: Directory containing genome files (must contain genome.fa and genome.gtf) (/LUCID/data/genome/)
- `reads_dir`: Directory containing RNA-seq files (*_1.fastq, *_2.fastq) (/LUCID/data/reads/)
- `tpm_threshold`: TPM threshold for identifying highly expressed genes
- `essential_fasta_path`: Path to essential proteins FASTA (/LUCID/data/references/Essential_protein_database.fasta)
- `orthofinder_tsv_path`: Path to OrthoFinder results TSV (/LUCID/data/references/Orthogroups.tsv)
- `output_dir`: Directory for output files (/LUCID/output/)
- `proteome_path`: Path to proteome FASTA of the fungus under study  (/LUCID/data/proteomes/your_proteome.fasta)

**Analysis approach:**
1. Calculates TPM values for all genes
2. Selects genes above the specified TPM threshold
3. Compares selected genes with conserved orthologs
4. Identifies those with homology to essential pathogenesis proteins

### For Non-Obligate Biotrophic Fungi

For hemibiotrophic and necrotrophic fungi, the pipeline uses differential expression analysis to compare gene expression during infection versus in vitro growth.

**Note:** This approach requires the `/LUCID/data/genome/coldata.csv` file described in the RNA-seq Processing section.

Execute:

```bash
Rscript /LUCID/script/lucid_differential.R \
  <genome_gtf_dir> \
  <reads_dir> \
  <log2fc_threshold> \
  <essential_fasta_path> \
  <orthofinder_tsv_path> \
  <output_dir> \
  <proteome_path>
```

**Parameters:**
- Same as above, except `log2fc_threshold` replaces `tpm_threshold`
- `log2fc_threshold`: Minimum log2 fold change value to consider a gene significantly upregulated during infection

**Analysis approach:**
1. Performs differential expression analysis between infection and control conditions using DESeq2
2. Selects genes with log2FC above threshold (upregulated during infection)
3. Compares selected genes with conserved orthologs
4. Identifies those with homology to essential pathogenesis proteins

---

## dsRNA Design

After identifying target genes (CEPs and CNAPs), the next step is to design dsRNA molecules for RNAi-based control. LUCID includes a specialized script for this purpose that identifies the optimal region for dsRNA amplification.

### Running the dsRNA Design Script

After target identification, use the `phase2-dsrna.sh` script to design optimal primers for dsRNA amplification:

```bash
bash LUCID/script/phase2-dsrna.sh <input.fasta> <output_dir> <transcriptome.fasta>
```

**Parameters**:
- `input.fasta`: FASTA file containing target gene sequences (typically CEPs.fasta or CNAPs.fasta)
- `output_dir`: Directory to store outputs
- `transcriptome.fasta`: Reference transcriptome for siRNA specificity analysis

**Example**:
```bash
bash LUCID/script/phase2-dsrna.sh ./output/CEPs.fasta ./output/dsrna_results/ ./data/genome/transcriptome.fasta
```

### What the Script Does

The `phase2-dsrna.sh` script performs several key functions:

1. **Builds a Bowtie database** from the transcriptome for siRNA specificity analysis
2. **Analyzes target sequences** to identify accessible regions for siRNA binding
3. **Identifies the optimal 200-300 bp window** with the highest density of effective siRNAs
4. **Designs PCR primers** to amplify this optimal window using Primer3

### Output Files

The script generates three main output files for each input sequence:

1. **siRNA Analysis Plot** (`./output_dir/plots/{sequence_id}_analysis.png`):
   - Graphical visualization of siRNA density and accessibility along the sequence
   - Highlighted selected window for dsRNA amplification

2. **Primer Design File** (`./output_dir/primers/{sequence_id}_primers.txt`):
   - Multiple primer pairs for amplifying the selected region
   - Details including sequence, TM, GC content, and product size

3. **JSON Results** (`./output_dir/{sequence_id}.json`):
   - Detailed information on siRNA positions, efficiency scores, and accessibility

### Analysis Algorithm

The script employs a modified version of sifi21 to:

1. Evaluate each potential siRNA for:
   - Target site accessibility
   - Strand-specific binding
   - End stability
   - Terminal nucleotide composition
   - Overall efficacy score

2. Identify a 200-300 bp window with the highest density of effective siRNAs

3. Design primers to amplify this region using Primer3 with parameters optimized for dsRNA production:
   - Optimal primer size: 20 bp (range: 18-25 bp)
   - Optimal TM: 60°C (range: 57-63°C)
   - GC content: 20-80%
   - Product size: 200-300 bp

---

## Results and Output Files

The analysis pipeline generates several important output files:

### Common Output Files

- `CEPs.fasta`: Sequences of identified Conserved Essential Proteins
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

### For dsRNA Design

- Primer design files (in `./output_dir/primers/`)
- siRNA analysis plots (in `./output_dir/plots/`)
- JSON results with detailed siRNA information

These files provide valuable information for designing RNA interference experiments targeting conserved, highly expressed genes in phytopathogenic fungi.

---

## Troubleshooting

### Common Issues

- **Error with OrthoFinder**: Ensure all proteome files are valid FASTA format and properly downloaded from ENSEMBL FUNGI
- **RNA-seq analysis error**: Verify that genome files are correctly named as `genome.fa` and `genome.gtf`
- **SRA download failure**: Check internet connection and SRA accession numbers
- **dsRNA design errors**: Ensure the input FASTA contains valid sequences and Python dependencies are installed

### Getting Help

For additional assistance, please refer to the manuscript or contact the authors.
