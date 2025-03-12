#!/bin/bash

# Script to run OrthoFinder analysis on proteomes in the data/proteomes directory
# Requirements: OrthoFinder installed and in PATH

# Default parameters
PROTEOMES_DIR="../data/proteomes"
OUTPUT_DIR="../data/references"
ORTHOFINDER_THREADS=4

# Function to display help message
show_help() {
    echo "Usage: run_orthofinder.sh [OPTIONS]"
    echo
    echo "This script runs OrthoFinder analysis on protein FASTA files in the proteomes directory."
    echo
    echo "Options:"
    echo "  -i, --input DIR          Input directory with proteome files (default: data/proteomes)"
    echo "  -o, --output DIR         Output directory for results (default: data/references)"
    echo "  -p, --threads NUM        Number of threads for OrthoFinder (default: 4)"
    echo "  -h, --help               Display this help message and exit"
    echo
    echo "Example usage:"
    echo "  ./run_orthofinder.sh -p 8"
    echo
    echo "Notes:"
    echo "  - Place all proteome FASTA files (.fa, .faa, or .fasta) in the input directory"
    echo "  - The script will copy the Orthogroups.tsv file to the output directory for LUCID"
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            PROTEOMES_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -p|--threads)
            ORTHOFINDER_THREADS="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Error: Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Check if OrthoFinder is installed
if ! command -v orthofinder &> /dev/null; then
    echo "Error: OrthoFinder is not installed or not in PATH."
    echo "Please install OrthoFinder: https://github.com/davidemms/OrthoFinder"
    exit 1
fi

# Check if input directory exists and contains files
if [ ! -d "$PROTEOMES_DIR" ]; then
    echo "Error: Proteomes directory does not exist: $PROTEOMES_DIR"
    echo "Please create it and add protein FASTA files (.fa, .faa, or .fasta)"
    exit 1
fi

# Count proteome files
num_proteomes=$(find "$PROTEOMES_DIR" -type f \( -name "*.fa" -o -name "*.faa" -o -name "*.fasta" \) | wc -l)

if [ "$num_proteomes" -eq 0 ]; then
    echo "Error: No proteome files found in $PROTEOMES_DIR"
    echo "Please add protein FASTA files (.fa, .faa, or .fasta) to this directory"
    exit 1
fi

echo "Found $num_proteomes proteome files in $PROTEOMES_DIR"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run OrthoFinder
echo "Running OrthoFinder with $ORTHOFINDER_THREADS threads..."
echo "This may take some time depending on the number and size of proteomes."

orthofinder -f "$PROTEOMES_DIR" -t "$ORTHOFINDER_THREADS" -S diamond

if [ $? -ne 0 ]; then
    echo "Error: OrthoFinder failed. Check the error messages above."
    exit 1
fi

echo "OrthoFinder analysis complete."

# Find the most recent results directory
latest_results=$(find "$PROTEOMES_DIR/OrthoFinder" -name "Results_*" -type d | sort | tail -n1)

if [ -d "$latest_results" ]; then
    # Copy Orthogroups.tsv to output directory
    cp "$latest_results/Orthogroups/Orthogroups.tsv" "$OUTPUT_DIR/"
    echo "Copied Orthogroups.tsv to $OUTPUT_DIR for use with LUCID"
    
    # Optional: Copy other useful files
    cp "$latest_results/Orthogroups/Orthogroups_SingleCopyOrthologues.txt" "$OUTPUT_DIR/" 2>/dev/null || true
    cp "$latest_results/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv" "$OUTPUT_DIR/" 2>/dev/null || true
    
    echo "Analysis complete. OrthoFinder results are in: $latest_results"
    echo "Key files have been copied to: $OUTPUT_DIR"
else
    echo "Warning: Could not find OrthoFinder results directory."
    echo "Please check: $PROTEOMES_DIR/OrthoFinder"
fi