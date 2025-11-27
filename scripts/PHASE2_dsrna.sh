#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input.fasta> <output_dir> <transcriptome.fasta>"
    echo "Example: $0 sequence.fa ./results/ ./transcriptome.fasta"
    exit 1
fi
# Convertir rutas relativas a absolutas
input_fasta=$(readlink -f "$1")
output_dir=$(readlink -f "$2")
transcriptome=$(readlink -f "$3")

# Check if input files exist
if [ ! -f "$input_fasta" ]; then
    echo "Error: Input fasta file does not exist: $input_fasta"
    exit 1
fi

if [ ! -f "$transcriptome" ]; then
    echo "Error: Transcriptome file does not exist: $transcriptome"
    exit 1
fi

# Create output directory and bowtie directory
mkdir -p "$output_dir"
bowtie_dir="${output_dir}/bowtie_db"
mkdir -p "$bowtie_dir"

# Build Bowtie database
echo "Building Bowtie database from transcriptome..."
bowtie-build "$transcriptome" "${bowtie_dir}/transcriptome_db"

# Create Python script for primer design
cat > "${output_dir}/get_primers.py" << 'EOL'
import argparse
from Bio import SeqIO
import pandas as pd
import primer3
from primer3 import bindings  # Esta línea es redundante, con import primer3 es suficiente
import matplotlib.pyplot as plt
import os
import tempfile
import json
import sifi21_plus
import numpy as np  # Faltaba esta importación

def process_fasta_and_run_sifi(query_sequences, output_dir, bowtie_db):
    temp_dir = tempfile.mkdtemp()
    json_contents = []
    
    try:
        for seq_record in SeqIO.parse(query_sequences, "fasta"):
            seq_name = seq_record.id
            seq_sequence = str(seq_record.seq)
            
            temp_file_path = os.path.join(temp_dir, f"{seq_name}.fasta")
            with open(temp_file_path, "w") as temp_file:
                temp_file.write(f">{seq_name}\n{seq_sequence}\n")
            
            json_output = os.path.join(output_dir, f"{seq_name}.json")
            try:
                sifiDesign = sifi21_plus.SifiPipeline(
                    bowtie_db=bowtie_db,
                    query_sequences=temp_file_path,
                    output_directory=output_dir,
                    mode=0,
                    sirna_size=21,
                    mismatches=0,
                    accessibility_check=True,
                    accessibility_window=8,
                    strand_check=True,
                    end_check=True,
                    end_stability_treshold=1.0,
                    target_site_accessibility_treshold=0.1,
                    terminal_check=True,
                    no_efficience=False
                )
                
                sifiDesign.run_pipeline
                
                if os.path.exists(json_output) and os.path.getsize(json_output) > 0:
                    with open(json_output) as f:
                        content = json.load(f)
                        if isinstance(content, list):
                            json_contents.extend(content)
                
            except Exception as e:
                print(f"Error processing {seq_name}: {str(e)}")
                continue

        if not json_contents:
            raise ValueError("No valid results obtained from sifi21INTA")
            
        return pd.DataFrame(json_contents)
        
    finally:
        for f in os.listdir(temp_dir):
            os.remove(os.path.join(temp_dir, f))
        os.rmdir(temp_dir)

def get_primers(query_sequences, output_dir, bowtie_db):
    """Process sequences and design primers for siRNA-accessible regions"""
    
    # Create output directories
    plots_dir = os.path.join(output_dir, "plots")
    primers_dir = os.path.join(output_dir, "primers")
    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(primers_dir, exist_ok=True)

    # Get siRNA data
    try:
        df = process_fasta_and_run_sifi(query_sequences, output_dir, bowtie_db)
    except Exception as e:
        print(f"Error processing sequences: {str(e)}")
        return {}

    results = {}
    
    # Process each sequence
    for seq_record in SeqIO.parse(query_sequences, "fasta"):
        query_name = seq_record.id
        sequence = str(seq_record.seq)
        seq_len = len(sequence)
        print(f"Processing sequence: {query_name}")

        try:
            # Filter effective siRNAs with accessibility
            filtered_df = df[
                (df['query_name'] == query_name) & 
                (df['is_efficient'] == True) & 
                (df['target_site_accessibility'] == True)
            ]
            
            if filtered_df.empty:
                print(f"No effective siRNAs found for {query_name}")
                continue

            # Calculate density vector
            density = [0] * seq_len
            accessibility_vector = np.zeros(seq_len)
            positions_used = set()

            # Fill vectors
            for _, row in filtered_df.iterrows():
                pos = int(row['sirna_position'])
                if pos < seq_len - 20:
                    # Update density
                    for i in range(21):
                        density[pos + i] += 1
                    
                    # Update accessibility (only once per position)
                    if pos not in positions_used:
                        accessibility_vector[pos] = float(row['accessibility_value'])
                        positions_used.add(pos)

            # Find optimal window
            max_density = 0
            start_pos = 0
            end_pos = 0
            
            for window_len in range(200, 401):
                for start in range(seq_len - window_len + 1):
                    end = start + window_len
                    window_density = sum(density[start:end])
                    if window_density > max_density:
                        max_density = window_density
                        start_pos = start
                        end_pos = end

            # Design primers
            primers = primer3.bindings.design_primers(
                seq_args={
                    'SEQUENCE_ID': query_name,
                    'SEQUENCE_TEMPLATE': sequence,
                    'SEQUENCE_INCLUDED_REGION': [start_pos, end_pos - start_pos]
                },
                global_args={
                    'PRIMER_OPT_SIZE': 20,
                    'PRIMER_MIN_SIZE': 18,
                    'PRIMER_MAX_SIZE': 25,
                    'PRIMER_OPT_TM': 60.0,
                    'PRIMER_MIN_TM': 57.0,
                    'PRIMER_MAX_TM': 63.0,
                    'PRIMER_MIN_GC': 20.0,
                    'PRIMER_MAX_GC': 80.0,
                    'PRIMER_PRODUCT_SIZE_RANGE': [[200, 300]]
                }
            )

            # Save primers
            primer_file = os.path.join(primers_dir, f"{query_name}_primers.txt")
            with open(primer_file, "w") as f:
                for idx, pair in enumerate(primers['PRIMER_PAIR']):
                    left = primers['PRIMER_LEFT'][idx]
                    right = primers['PRIMER_RIGHT'][idx]
                    f.write(f"PAIR {idx + 1}:\n")
                    f.write(f"  LEFT_PRIMER: {left['SEQUENCE']}, TM: {left['TM']}, GC%: {left['GC_PERCENT']}\n")
                    f.write(f"  RIGHT_PRIMER: {right['SEQUENCE']}, TM: {right['TM']}, GC%: {right['GC_PERCENT']}\n")
                    f.write(f"  PRODUCT_SIZE: {pair['PRODUCT_SIZE']}\n\n")

            # Create plot
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10), height_ratios=[2, 1])
            
            # Upper panel: siRNA density
            ax1.plot(range(seq_len), density, 'b-', linewidth=2, label='Effective siRNAs')
            ax1.axvspan(start_pos, end_pos, alpha=0.2, color='orange', 
                        label=f'Selected window ({start_pos}-{end_pos})')
            ax1.set_title(f'siRNA Density and Accessibility - {query_name}')
            ax1.set_ylabel('Number of effective siRNAs')
            ax1.legend()
            ax1.grid(True, alpha=0.3)

            # Lower panel: accessibility
            non_zero_mask = accessibility_vector > 0
            if np.any(non_zero_mask):
                positions = np.where(non_zero_mask)[0]
                values = accessibility_vector[non_zero_mask]
                ax2.plot(positions, values, 'r-', linewidth=2, label='Accessibility')
                ax2.axhline(y=0.1, color='r', linestyle='--', alpha=0.5, 
                            label='Accessibility threshold')
            
            ax2.set_xlabel('Position')
            ax2.set_ylabel('Accessibility value')
            ax2.legend()
            ax2.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            # Save plot
            plot_file = os.path.join(plots_dir, f"{query_name}_analysis.png")
            plt.savefig(plot_file)
            plt.close()

            # Store results
            results[query_name] = {
                'density': density,
                'window': (start_pos, end_pos),
                'primers': primers,
                'plot_file': plot_file,
                'primer_file': primer_file
            }

        except Exception as e:
            print(f"Error processing {query_name}: {str(e)}")
            continue

    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Design primers for siRNA-accessible regions")
    parser.add_argument("--query_sequences", required=True, help="Input FASTA file")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--bowtie_db", required=True, help="Bowtie database path")
    args = parser.parse_args()
    get_primers(**vars(args))
EOL

# Run the analysis
echo "Processing $input_fasta..."
python3 "${output_dir}/get_primers.py" \
    --query_sequences "$input_fasta" \
    --output_dir "$output_dir" \
    --bowtie_db "${bowtie_dir}/transcriptome_db"

# Check if the execution was successful
if [ $? -eq 0 ]; then
    echo "Analysis completed successfully. Results are in: $output_dir"
    echo "Bowtie database is in: $bowtie_dir"
else
    echo "Error during analysis"
    exit 1
fi

# Clean up Bowtie files if needed
rm -rf "$bowtie_dir"  # Uncomment this line if you want to remove the Bowtie database after analysis
