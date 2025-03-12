#!/bin/bash

# Directorio de salida
output_dir="/mnt/data2/STORAGE/lucia_jimenez"

# Lista de accesiones a descargar (modifica esta variable con tus accesiones)
accessions=("SRR6924534" "SRR6924535" "SRR6924536")

# Función para obtener la URL correcta de ENA
get_ena_url() {
    local acc=$1
    local prefix=${acc:0:6}
    
    # Tomar el último dígito de la accesión y formatearlo como "00X"
    local last_digit=${acc: -1}
    local subdir="00${last_digit}"

    echo "https://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/${subdir}/${acc}"
}

# Función para verificar si un archivo existe en el servidor
check_file_exists() {
    local url=$1
    wget --spider --quiet "$url"
    return $?  # Devuelve 0 si existe, 1 si no
}

# Función para descargar una sola accesión
download_accession() {
    local acc=$1
    local outdir=$2
    echo "Procesando $acc..."

    # Obtener la URL base
    base_url=$(get_ena_url "$acc")

    # Construir URLs de los archivos
    url_1="${base_url}/${acc}_1.fastq.gz"
    url_2="${base_url}/${acc}_2.fastq.gz"

    # Descargar Read 1
    if check_file_exists "$url_1"; then
        echo "Descargando ${acc}_1.fastq.gz..."
        wget --no-check-certificate "$url_1" -O "${outdir}/${acc}_1.fastq.gz"
    else
        echo "Error: ${acc}_1.fastq.gz no encontrado."
    fi

    # Descargar Read 2
    if check_file_exists "$url_2"; then
        echo "Descargando ${acc}_2.fastq.gz..."
        wget --no-check-certificate "$url_2" -O "${outdir}/${acc}_2.fastq.gz"
    else
        echo "Error: ${acc}_2.fastq.gz no encontrado."
    fi
}

# Crear directorio de salida si no existe
mkdir -p "$output_dir"

# Procesar accesiones definidas en la variable
for acc in "${accessions[@]}"; do
    download_accession "$acc" "$output_dir"
    echo "------------------------------------"
done

echo "Todas las descargas completadas."
