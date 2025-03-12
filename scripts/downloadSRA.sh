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
accessions=("SRR6924534" "SRR6924535" "SRR6924536")
# Function to get the correct ENA URL
get_ena_url() {
    local acc=$1
    local prefix=${acc:0:6}
    
    # Take the last digit of the accession and format it as "00X"
    local last_digit=${acc: -1}
    local subdir="00${last_digit}"
    echo "https://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/${subdir}/${acc}"
}
# Function to verify if a file exists on the server
check_file_exists() {
    local url=$1
    wget --spider --quiet "$url"
    return $?  # Returns 0 if exists, 1 if not
}
# Function to download a single accession
download_accession() {
    local acc=$1
    local outdir=$2
    echo "Processing $acc..."
    # Get the base URL
    base_url=$(get_ena_url "$acc")
    # Build file URLs
    url_1="${base_url}/${acc}_1.fastq.gz"
    url_2="${base_url}/${acc}_2.fastq.gz"
    # Download Read 1
    if check_file_exists "$url_1"; then
        echo "Downloading ${acc}_1.fastq.gz..."
        wget --no-check-certificate "$url_1" -O "${outdir}/${acc}_1.fastq.gz"
    else
        echo "Error: ${acc}_1.fastq.gz not found."
    fi
    # Download Read 2
    if check_file_exists "$url_2"; then
        echo "Downloading ${acc}_2.fastq.gz..."
        wget --no-check-certificate "$url_2" -O "${outdir}/${acc}_2.fastq.gz"
    else
        echo "Error: ${acc}_2.fastq.gz not found."
    fi
}
# Create output directory if it doesn't exist
mkdir -p "$output_dir"
# Process accessions defined in the variable
for acc in "${accessions[@]}"; do
    download_accession "$acc" "$output_dir"
    echo "------------------------------------"
done
echo "All downloads completed."
