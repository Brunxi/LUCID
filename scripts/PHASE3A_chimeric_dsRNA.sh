#!/bin/bash

# PHASE 3A: Chimeric Multi-Target dsRNA Design

cat << "EOF"
╔════════════════════════════════════════════════════════════════╗
║                  LUCID - PHASE 3A                              ║
║        Chimeric dsRNA Multi-Target Design                      ║
╚════════════════════════════════════════════════════════════════╝
EOF

# Detect dsRNAmax in LUCID/dsRNAmax/
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DSRNAMAX="${SCRIPT_DIR}/../dsRNAmax/dsRNAmax"

if [ ! -f "$DSRNAMAX" ] || [ ! -x "$DSRNAMAX" ]; then
    echo "❌ Error: dsRNAmax not found at: $DSRNAMAX"
    echo ""
    echo "Install dsRNAmax:"
    echo "  cd $(dirname $SCRIPT_DIR)"
    echo "  git clone https://github.com/sfletc/dsRNAmax.git"
    echo "  cd dsRNAmax && go build"
    echo ""
    exit 1
fi

echo "✓ Using dsRNAmax: $DSRNAMAX"

if [ "$#" -lt 2 ]; then
    echo ""
    echo "Usage: $0 <multiple_targets.fasta> <off_target_transcriptome.fasta> [OPTIONS]"
    echo ""
    echo "Design a SINGLE optimized dsRNA targeting MULTIPLE genes (CEPs/CNAPs)"
    echo "while avoiding off-targets in beneficial organisms."
    echo ""
    echo "Required arguments:"
    echo "  multiple_targets.fasta         : FASTA with 2+ target sequences"
    echo "  off_target_transcriptome.fasta : Beneficial organism transcriptome"
    echo ""
    echo "Optional arguments:"
    echo "  --construct-len N   : dsRNA length (default: 300)"
    echo "  --kmer-len N        : Kmer length (default: 21)"
    echo "  --ot-kmer-len N     : Off-target kmer length (default: 21)"
    echo "  --iterations N      : Optimization iterations (default: 100)"
    echo "  --output-dir DIR    : Output directory (default: ./output/phase3a/)"
    echo ""
    echo "Example:"
    echo "  $0 selected_CEPs.fasta Cucumis_melo.cds.fa"
    echo ""
    exit 1
fi

# Argumentos obligatorios
TARGET_FASTA=$(readlink -f "$1")
OFFTARGET_FASTA=$(readlink -f "$2")
shift 2

# Parámetros por defecto
CONSTRUCT_LEN=300
KMER_LEN=21
OT_KMER_LEN=21
ITERATIONS=100
OUTPUT_DIR="./output/phase3a"

# Parsear opcionales
while [[ $# -gt 0 ]]; do
    case "$1" in
        --construct-len) CONSTRUCT_LEN="$2"; shift 2 ;;
        --kmer-len) KMER_LEN="$2"; shift 2 ;;
        --ot-kmer-len) OT_KMER_LEN="$2"; shift 2 ;;
        --iterations) ITERATIONS="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        *) echo "Error: Unknown option: $1"; exit 1 ;;
    esac
done

# Verificar archivos
if [ ! -f "$TARGET_FASTA" ]; then
    echo "❌ Error: Target file not found: $TARGET_FASTA"
    exit 1
fi

if [ ! -f "$OFFTARGET_FASTA" ]; then
    echo "❌ Error: Off-target file not found: $OFFTARGET_FASTA"
    exit 1
fi

# Crear directorio
mkdir -p "$OUTPUT_DIR"

# Analizar input
NUM_TARGETS=$(grep -c "^>" "$TARGET_FASTA")

echo ""
echo "═══════════════════════════════════════════════════════════"
echo "Configuration:"
echo "───────────────────────────────────────────────────────────"
echo "  Number of targets:       $NUM_TARGETS genes"
echo "  Target file:             $(basename $TARGET_FASTA)"
echo "  Off-target file:         $(basename $OFFTARGET_FASTA)"
echo "  dsRNA construct length:  ${CONSTRUCT_LEN} bp"
echo "  Target kmer length:      ${KMER_LEN} nt"
echo "  Off-target kmer length:  ${OT_KMER_LEN} nt"
echo "  Optimization iterations: ${ITERATIONS}"
echo "═══════════════════════════════════════════════════════════"
echo ""

echo "🎯 Target genes:"
grep "^>" "$TARGET_FASTA" | sed 's/^>/  • /'
echo ""

echo "🔬 Running dsRNAmax chimeric design..."
echo ""

# Ejecutar dsRNAmax
"$DSRNAMAX" \
    -targets "$TARGET_FASTA" \
    -offTargets "$OFFTARGET_FASTA" \
    -kmerLen $KMER_LEN \
    -otKmerLen $OT_KMER_LEN \
    -constructLen $CONSTRUCT_LEN \
    -iterations $ITERATIONS \
    -csv "${OUTPUT_DIR}/chimeric_dsrna_stats.csv" \
    2>&1 | tee "${OUTPUT_DIR}/chimeric_design.log"

EXIT_CODE=${PIPESTATUS[0]}

if [ $EXIT_CODE -eq 0 ] && grep -q "dsRNA sense-arm sequence" "${OUTPUT_DIR}/chimeric_design.log"; then
    echo ""
    echo "═══════════════════════════════════════════════════════════"
    echo "✅ SUCCESS: Chimeric dsRNA designed!"
    echo "═══════════════════════════════════════════════════════════"
    echo ""
    
    # Extraer secuencia
    grep "dsRNA sense-arm sequence" -A 1 "${OUTPUT_DIR}/chimeric_design.log" | tail -1 > "${OUTPUT_DIR}/chimeric_sequence.txt"
    
    CHIMERIC_SEQ=$(cat "${OUTPUT_DIR}/chimeric_sequence.txt")
    CHIMERIC_LEN=${#CHIMERIC_SEQ}
    GC_CONTENT=$(grep "dsRNA sense-arm sequence" "${OUTPUT_DIR}/chimeric_design.log" | grep -oP '\d+\.\d+(?=% GC)')
    KMERS_REMOVED=$(grep "Total off-target-matching kmers removed:" "${OUTPUT_DIR}/chimeric_design.log" | grep -oP '\d+$')
    MEDIAN_HITS=$(grep "Median of kmer hits" "${OUTPUT_DIR}/chimeric_design.log" | grep -oP '\d+')
    
    # Crear FASTA
    echo ">Chimeric_dsRNA_${NUM_TARGETS}targets_${CHIMERIC_LEN}bp_OTvalidated" > "${OUTPUT_DIR}/chimeric_dsrna.fasta"
    echo "$CHIMERIC_SEQ" >> "${OUTPUT_DIR}/chimeric_dsrna.fasta"
    
    echo "📊 Chimeric Design Summary:"
    echo "───────────────────────────────────────────────────────────"
    echo "  Chimeric dsRNA length:    ${CHIMERIC_LEN} bp"
    echo "  GC content:               ${GC_CONTENT}%"
    echo "  Targets combined:         ${NUM_TARGETS} genes"
    echo "  Median kmer hits/target:  ${MEDIAN_HITS}"
    echo "  Off-target kmers removed: ${KMERS_REMOVED}"
    echo ""
    echo "✓ This SINGLE dsRNA targets ALL ${NUM_TARGETS} genes simultaneously"
    echo ""
    
    echo "📈 Coverage per target:"
    echo "───────────────────────────────────────────────────────────"
    grep -A $((NUM_TARGETS + 2)) "TARGET SEQUENCE HEADER" "${OUTPUT_DIR}/chimeric_design.log" | \
        grep -v "^+" | grep -v "TARGET SEQUENCE" | grep -v "^$" | \
        awk '{printf "  • %-30s %s matches\n", $1, $2}'
    echo ""
    
    echo "📁 Output files:"
    echo "───────────────────────────────────────────────────────────"
    echo "  ✓ ${OUTPUT_DIR}/chimeric_dsrna.fasta"
    echo "  ✓ ${OUTPUT_DIR}/chimeric_sequence.txt"
    echo "  ✓ ${OUTPUT_DIR}/chimeric_dsrna_stats.csv"
    echo "  ✓ ${OUTPUT_DIR}/chimeric_design.log"
    echo ""
    
else
    echo ""
    echo "═══════════════════════════════════════════════════════════"
    echo "❌ ERROR: Chimeric design failed"
    echo "═══════════════════════════════════════════════════════════"
    echo ""
    echo "Try:"
    echo "  --construct-len 250  (reduce length)"
    echo "  --kmer-len 19        (shorter kmers)"
    echo "  --ot-kmer-len 18     (more permissive)"
    echo ""
    exit 1
fi

echo "═══════════════════════════════════════════════════════════"
echo "PHASE 3A Complete"
echo "═══════════════════════════════════════════════════════════"