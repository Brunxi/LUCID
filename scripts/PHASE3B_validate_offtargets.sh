#!/bin/bash

# PHASE 3B: dsRNA Off-Target Kmer Validation using dsRNAmax

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DSRNAMAX="${SCRIPT_DIR}/../dsRNAmax/dsRNAmax"

if [ ! -f "$DSRNAMAX" ]; then
    echo "Error: dsRNAmax not found at: $DSRNAMAX"
    exit 1
fi

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <dsrna.fasta> <offtarget.fasta> [--kmer-len N] [--output-dir DIR]"
    exit 1
fi

DSRNA_FASTA=$(readlink -f "$1")
OFFTARGET_FASTA=$(readlink -f "$2")
shift 2

KMER_LEN=21
OUTPUT_DIR="./output/phase3b"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --kmer-len) KMER_LEN="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

mkdir -p "$OUTPUT_DIR"

OFFTARGET_NAME=$(basename "$OFFTARGET_FASTA" | sed 's/\.[^.]*$//')

# Split FASTA and process each sequence
grep "^>" "$DSRNA_FASTA" | sed 's/^>//' | while read -r TRANSCRIPT_NAME; do
    
    echo "========================================"
    echo "Analyzing: ${TRANSCRIPT_NAME}"
    echo "========================================"
    
    # Extract this sequence
    awk -v name="$TRANSCRIPT_NAME" '
        /^>/ { if ($0 ~ name) {print; getline; print; exit} }
    ' "$DSRNA_FASTA" > "${OUTPUT_DIR}/temp_${TRANSCRIPT_NAME}.fasta"
    
    DSRNA_SEQ=$(grep -v "^>" "${OUTPUT_DIR}/temp_${TRANSCRIPT_NAME}.fasta" | tr -d '\n' | tr -d ' ')
    DSRNA_LEN=${#DSRNA_SEQ}
    
    if [ "$DSRNA_LEN" -lt "$KMER_LEN" ]; then
        echo "Skipping: sequence too short (${DSRNA_LEN} bp)"
        echo ""
        rm -f "${OUTPUT_DIR}/temp_${TRANSCRIPT_NAME}.fasta"
        continue
    fi
    
    NUM_KMERS=$((DSRNA_LEN - KMER_LEN + 1))
    
    echo "Length: ${DSRNA_LEN} bp, ${NUM_KMERS} kmers"
    echo ""
    
    # Run dsRNAmax to check off-targets
    echo "Running dsRNAmax validation..."
    "$DSRNAMAX" \
        -targets "${OUTPUT_DIR}/temp_${TRANSCRIPT_NAME}.fasta" \
        -offTargets "$OFFTARGET_FASTA" \
        -kmerLen $KMER_LEN \
        -constructLen 100 \
        > "${OUTPUT_DIR}/${TRANSCRIPT_NAME}_dsrnamax_${OFFTARGET_NAME}.log" 2>&1
    
    # Extract results from dsRNAmax log
    TOTAL_KMERS=$(grep -oP '\d+(?= target kmers loaded)' "${OUTPUT_DIR}/${TRANSCRIPT_NAME}_dsrnamax_${OFFTARGET_NAME}.log" | head -n1)
    KMERS_REMOVED=$(grep "Total off-target-matching kmers removed:" "${OUTPUT_DIR}/${TRANSCRIPT_NAME}_dsrnamax_${OFFTARGET_NAME}.log" | grep -oP '\d+$')
    
    if [ -z "$KMERS_REMOVED" ]; then
        KMERS_REMOVED=0
    fi
    
    CLEAN_KMERS=$((TOTAL_KMERS - KMERS_REMOVED))
    
    echo "Results:"
    echo "  Total kmers: ${TOTAL_KMERS}"
    echo "  Off-target matches: ${KMERS_REMOVED}"
    echo "  Clean kmers: ${CLEAN_KMERS}"
    echo ""
    
    if [ "$KMERS_REMOVED" -eq 0 ]; then
        echo "Status: SAFE - No off-target matches in ${OFFTARGET_NAME}"
        cp "${OUTPUT_DIR}/temp_${TRANSCRIPT_NAME}.fasta" \
           "${OUTPUT_DIR}/${TRANSCRIPT_NAME}_validated_no_${OFFTARGET_NAME}.fasta"
    else
        echo "Status: WARNING - ${KMERS_REMOVED} off-target matches in ${OFFTARGET_NAME}"
        echo ""
        echo "Percentage compromised: $(awk "BEGIN {printf \"%.1f\", 100*$KMERS_REMOVED/$TOTAL_KMERS}")%"
    fi
    
    echo ""
    echo "Output files:"
    echo "  - ${TRANSCRIPT_NAME}_dsrnamax_${OFFTARGET_NAME}.log"
    if [ "$KMERS_REMOVED" -eq 0 ]; then
        echo "  - ${TRANSCRIPT_NAME}_validated_no_${OFFTARGET_NAME}.fasta"
    fi
    echo ""
    
    rm -f "${OUTPUT_DIR}/temp_${TRANSCRIPT_NAME}.fasta"
done

echo "========================================"
echo "Analysis complete for all sequences"
echo "Results in: ${OUTPUT_DIR}/"
echo "========================================"
