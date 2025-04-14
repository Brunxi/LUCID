# Load necessary libraries
library(dplyr)       # For data manipulation
library(data.table)  # For fast data import and manipulation
library(SummarizedExperiment) # For handling genomic data
library(lattice)     # For visualization
library(Rsubread)    # For RNA-seq alignment and counting
library(ggplot2)     # For visualization
library(GenomicRanges) # For genomic interval manipulation
library(tidyr)       # For data reshaping
library(DT)          # For interactive tables
library(httr)        # For HTTP requests
library(Biostrings)  # For biological sequence manipulation
library(jsonlite)    # For JSON parsing

#' Check if a UniProt ID mapping job is complete
#'
#' This function polls the UniProt API to check if an ID mapping job is ready
#'
#' @param jobId Character string with the job identifier
#' @return Logical indicating if the job is ready (TRUE) or not (FALSE)
isJobReady <- function(jobId) {
    pollingInterval <- 5  # Seconds between API calls
    nTries <- 20         # Maximum number of attempts
    
    for (i in 1:nTries) {
        url <- paste0("https://rest.uniprot.org/idmapping/status/", jobId)
        r <- GET(url = url, accept("application/json"))
        status <- httr::content(r, as = "parsed", type = "application/json")
        
        # Check if results or failed IDs are available
        if (!is.null(status[["results"]]) || !is.null(status[["failedIds"]])) {
            return(TRUE)
        }
        
        # Check for error messages
        if (!is.null(status[["messages"]])) {
            print(status[["messages"]])
            return(FALSE)
        }
        
        # Wait before next attempt
        Sys.sleep(pollingInterval)
    }
    
    return(FALSE)
}

#' Get the results URL for a UniProt ID mapping job
#'
#' This function converts a redirect URL to a stream URL for retrieving results
#'
#' @param redirectURL Character string with the redirect URL
#' @return Character string with the stream URL
getResultsURL <- function(redirectURL) {
    if (grepl("/idmapping/results/", redirectURL, fixed = TRUE)) {
        url <- gsub("/idmapping/results/", "/idmapping/stream/", redirectURL)
    } else {
        url <- gsub("/results/", "/results/stream/", redirectURL)
    }
    return(url)
}

#' Map Ensembl protein IDs to UniProt IDs
#'
#' This function maps Ensembl protein IDs to UniProt IDs using the UniProt API
#' with batch processing for large sets of IDs
#'
#' @param protein_ids Vector of Ensembl protein IDs
#' @param batch_size Integer specifying the batch size for API requests (default: 500)
#' @return Data frame containing the mapping results
map_to_uniprot_ids <- function(protein_ids, batch_size = 500) {
    # Return empty dataframe if no protein IDs are provided
    if (length(protein_ids) == 0) {
        return(data.frame())
    }
    
    total_proteins <- length(protein_ids)
    results <- data.frame()
    
    # Process in batches to avoid API limitations
    for (i in seq(1, total_proteins, by = batch_size)) {
        batch_ids <- protein_ids[i:min(i + batch_size - 1, total_proteins)]
        files <- list(
            from = "Ensembl_Genomes_Protein",
            to = "UniProtKB",
            ids = paste(batch_ids, collapse = ",")
        )
        
        tryCatch({
            # Submit the ID mapping job
            r <- POST(
                url = "https://rest.uniprot.org/idmapping/run",
                body = files,
                encode = "multipart",
                accept("application/json")
            )
            
            submission <- httr::content(r, as = "parsed", type = "application/json")
            
            # Check if the job is complete
            if (isJobReady(submission[["jobId"]])) {
                # Get job details
                url <- paste0("https://rest.uniprot.org/idmapping/details/", submission[["jobId"]])
                r <- GET(url = url, accept("application/json"))
                details <- httr::content(r, as = "parsed", type = "application/json")
                
                # Get results in TSV format
                url <- getResultsURL(details[["redirectURL"]])
                url <- paste0(url, "?format=tsv")
                r <- GET(url = url, accept("text/tsv"))
                resultsTable <- read.delim(text = rawToChar(r$content), sep = "\t")
                
                # Combine with previous results
                results <- rbind(results, resultsTable)
            }
        }, error = function(e) {
            warning("Error processing batch ", i, ": ", e$message)
        })
        
        # Add a small pause to avoid overloading the API
        Sys.sleep(1)
    }
    
    return(results)
}

#' Calculate TPM (Transcripts Per Million) from a count matrix
#'
#' @param counts Matrix of read counts
#' @param gene_lengths Vector of gene lengths
#' @return Matrix of TPM values
calculate_tpm <- function(counts, gene_lengths) {
    # Step 1: Normalize by gene length (get RPK - reads per kilobase)
    rpk <- counts / gene_lengths
    
    # Step 2: Scale by million (get TPM)
    scaling_factor <- colSums(rpk) / 1e6
    tpm <- sweep(rpk, 2, scaling_factor, "/")
    
    return(tpm)
}

#' Process gene expression data using TPM normalization
#'
#' This function aligns reads, counts features, and identifies highly expressed genes
#' based on TPM (Transcripts Per Million) values
#'
#' @param genome_gtf_dir Directory containing genome GTF and FASTA files
#' @param reads_dir Directory containing FASTQ files
#' @param output_dir Directory for output files
#' @param threshold TPM threshold for highly expressed genes
#' @return Vector of highly expressed gene names
process_expression_tpm <- function(genome_gtf_dir, reads_dir, output_dir, threshold) {
    
    # Define the full index name with path to the genome directory
    basename_index <- file.path(genome_gtf_dir, "reference_index")
    
    # Build the index using the full path if it doesn't exist
    if (!file.exists(paste0(basename_index, ".00.b.array"))) {
        buildindex(basename = basename_index, reference = file.path(genome_gtf_dir, "genome.fa"))
    }
    
    # Get paths for R1 and R2 reads
    reads1 <- list.files(reads_dir, pattern = "*_1\\.fastq$", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
    reads2 <- list.files(reads_dir, pattern = "*_2\\.fastq$", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
    
    # Ensure there are matching R1 and R2 files
    if (length(reads1) != length(reads2)) {
        stop("The number of R1 and R2 files does not match.")
    }
    
    # Alignment, skipping if BAM files already exist
    for (i in seq_along(reads1)) {
        bam_file <- file.path(reads_dir, paste0(basename(reads1[i]), ".subread.BAM"))
        if (!file.exists(bam_file)) {
            align(index = file.path(genome_gtf_dir, "reference_index"), 
                useAnnotation = TRUE, 
                annot.ext = file.path(genome_gtf_dir, "genome.gtf"), 
                isGTF = TRUE, 
                readfile1 = reads1[i],
                readfile2 = reads2[i],
                type = "rna", 
                nthreads = 4, 
                output_file = bam_file,
                minFragLength = 50, 
                maxFragLength = 600)
        } else {
            message("BAM file for ", reads1[i], " already exists. Skipping alignment.")
        }
    }
    
    # Get paths for BAM files
    bam_files <- list.files(reads_dir, pattern = "*.subread.BAM$", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
    
    # Count reads per feature with gene length information
    fc <- featureCounts(files = bam_files, 
                      annot.ext = file.path(genome_gtf_dir, "genome.gtf"), 
                      isGTFAnnotationFile = TRUE, 
                      isPairedEnd = TRUE, 
                      useMetaFeatures = TRUE, 
                      GTF.featureType = "exon",
                      GTF.attrType = "gene_id")
    
    # Get gene lengths
    gene_lengths <- fc$annotation$Length
    names(gene_lengths) <- fc$annotation$GeneID
    
     # Calculate TPM
    tpm_matrix <- calculate_tpm(fc$counts, gene_lengths)
    
    # Write the TPM matrix
    write.table(tpm_matrix, file = file.path(output_dir, "tpm_matrix.tsv"), 
                sep = "\t", quote = FALSE, col.names = NA)
    
    # Calculate mean TPM across all samples
    mean_tpm <- rowMeans(tpm_matrix)
    
    # Sort genes by mean TPM
    sorted_tpm <- sort(mean_tpm, decreasing = TRUE)
    
    # Filter genes with TPM greater than the threshold
    highly_expressed_genes <- names(sorted_tpm[sorted_tpm > threshold])
    
    # Save highly expressed genes
    highly_expressed_gene_names <- gsub("^gene:", "", highly_expressed_genes)
    write.table(highly_expressed_gene_names, 
               file = file.path(output_dir, "highly_expressed_genes.tsv"), 
               sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    message(sprintf("Found %d highly expressed genes", length(highly_expressed_genes)))
    
    return(highly_expressed_genes)
}

#' Create a dictionary of core genes conserved in all phytopathogenic fungi
#'
#' @param file_path Path to OrthoFinder TSV results
#' @return List of orthogroups with corresponding genes for each organism
create_core_gene_dict <- function(file_path) {
    # Load data from the TSV file
    data <- read.delim(file_path, na.strings = c("", " "))
    
    # Filter data to remove rows with NA
    filtered_data <- data[complete.cases(data), ]
    
    # Convert comma-separated genes to individual rows
    df_long <- filtered_data %>%
        pivot_longer(-Orthogroup, names_to = "organism", values_to = "genes") %>%
        separate_rows(genes, sep = ",\\s*")
    
    # Create a new column with the corresponding orthogroup for each gene
    df_long <- df_long %>%
        left_join(filtered_data, by = "Orthogroup") %>%
        select(organism, genes, Orthogroup)
    
    # Structure the content of common orthogroups for all organisms
    dict_core <- list()
    
    # Iterate over each row of the data frame
    for (i in 1:nrow(df_long)) {
        # Get the values of organism, genes, and Orthogroup from the current row
        organism <- df_long$organism[i]
        genes <- unlist(strsplit(df_long$genes[i], ", "))
        orthogroup <- df_long$Orthogroup[i]
        
        # Check if the orthogroup exists in the dictionary
        if (!(orthogroup %in% names(dict_core))) {
            dict_core[[orthogroup]] <- list()  # Create a subdictionary for the orthogroup
        }
        
        # Check if the organism exists in the subdictionary of the orthogroup
        if (!(organism %in% names(dict_core[[orthogroup]]))) {
            dict_core[[orthogroup]][[organism]] <- list()  # Create a list for the genes of the organism
        }
        
        # Add genes to the organism's subdictionary
        dict_core[[orthogroup]][[organism]] <- c(dict_core[[orthogroup]][[organism]], genes)
    }
    
    return(dict_core)
}

#' Create a dictionary with protein ID as key and gene as value from a transcriptome
#'
#' @param file_path Path to transcriptome FASTA file
#' @return List mapping protein IDs to gene IDs
create_protein_gene_dict <- function(file_path) {
    # Load the multifasta
    sequences <- readDNAStringSet(file_path)
    
    # Create an empty dictionary
    dictionary <- list()
    
    # Iterate through sequences and extract information from headers
    for (i in 1:length(sequences)) {
        header <- names(sequences)[i]
        
        # Extract gene name and transcript value
        gene <- sub(".*gene:([^ ]+).*", "\\1", header)
        transcript <- sub("^([^ ]+).*", "\\1", header)
        
        # Add the key-value pair to the dictionary
        dictionary[[transcript]] <- gene
    }
    
    return(dictionary)
}

#' Get protein IDs for a list of genes
#'
#' @param dictionary Dictionary mapping protein IDs to gene IDs
#' @param gene_list List of gene IDs
#' @return Data frame with genes and their corresponding protein IDs
get_protein_ids <- function(dictionary, gene_list) {
    # Create an empty data frame to store results
    genes_with_protein <- data.frame(gene = character(), protein = character(), stringsAsFactors = FALSE)
    
    # Iterate through genes and find corresponding transcript identifiers
    for (gene in gene_list) {
        if (gene %in% unlist(dictionary)) {
            # Get all transcript identifiers for the gene
            transcript_values <- names(dictionary[dictionary == gene])
            
            # Create a temporary data frame with the gene and each of its transcripts
            temp_df <- data.frame(gene = rep(gene, length(transcript_values)), 
                               protein = transcript_values, 
                               stringsAsFactors = FALSE)
            
            # Add rows to the main data frame
            genes_with_protein <- rbind(genes_with_protein, temp_df)
        }
    }
    
    return(genes_with_protein)
}

#' Filter a proteome to get protein sequences of interest
#'
#' @param protein_list List of protein IDs to filter for
#' @param proteome_path Path to proteome FASTA file
#' @return DNAStringSet containing filtered protein sequences
filter_protein_sequences <- function(protein_list, proteome_path) {
    # Load protein sequences from the multifasta file
    sequences <- readDNAStringSet(proteome_path)
    
    filtered_sequences <- DNAStringSet()  # Initialize an empty sequence set
    
    # Iterate through all sequence headers
    for (i in 1:length(names(sequences))) {
        header <- names(sequences)[i]
        
        # Extract protein ID from header
        protein_id <- sub("^([^ ]+).*", "\\1", header)
        
        # Check if protein ID is in the list of proteins to filter
        if (!is.na(protein_id) && protein_id %in% protein_list) {
            # Add sequence to the filtered protein sequences set
            filtered_sequences <- c(filtered_sequences, sequences[i])
        }
    }
    
    return(filtered_sequences)
}

#' Check if a protein file is valid
#'
#' @param input_path Path to the protein FASTA file
#' @return Logical indicating if the file is valid (TRUE) or not (FALSE)
clean_protein_file <- function(input_path) {
    tryCatch({
        # Read sequences
        seqs <- readAAStringSet(input_path)
        message(sprintf("Read %d sequences from %s", length(seqs), input_path))
        return(TRUE)
    }, error = function(e) {
        message(sprintf("Error reading %s: %s", input_path, e$message))
        return(FALSE)
    })
}

#' Identify highly expressed Core Effector Proteins (CEPs) from RNA-seq data
#'
#' This function processes RNA-seq data to identify highly expressed core effector proteins
#' and conserved non-annotated proteins (CNAPs) in phytopathogenic fungi
#'
#' @param genome_gtf_dir Directory containing genome GTF and FASTA files
#' @param reads_dir Directory containing RNA-seq FASTQ files
#' @param threshold TPM threshold for highly expressed genes
#' @param essential_fasta_path Path to FASTA file with essential proteins
#' @param orthofinder_tsv_path Path to OrthoFinder TSV results
#' @param output_dir Directory for output files
#' @param proteome_path Path to proteome FASTA file
#' @return List with CEPs and CNAPs sequences, IDs, and summary
get_up_CEPs <- function(genome_gtf_dir, reads_dir, threshold, 
                       essential_fasta_path, orthofinder_tsv_path, output_dir, proteome_path) {
    # 1. Get highly expressed genes from RNA-seq analysis
    message("Processing gene expression...")
    highly_expressed_genes <- process_expression_tpm(genome_gtf_dir, reads_dir, output_dir, threshold)
    
    # Clean gene IDs (remove "gene:" prefix)
    highly_expressed_genes <- gsub("^gene:", "", highly_expressed_genes)
    message(sprintf("Found %d highly expressed genes", length(highly_expressed_genes)))
    
    # 2. Create dictionary of core genes
    message("Creating core gene dictionary...")
    core_dictionary <- create_core_gene_dict(orthofinder_tsv_path)
    message(sprintf("Created dictionary with %d entries", length(core_dictionary)))
    
    # 3. Get highly expressed core proteins
    message("Getting core proteins...")
    proteome_dictionary <- create_protein_gene_dict(proteome_path)
    protein_list <- get_protein_ids(proteome_dictionary, highly_expressed_genes)
    
    # Debug information
    message(sprintf("Number of proteins found: %d", nrow(protein_list)))
    message(sprintf("Highly expressed genes: %s", paste(head(highly_expressed_genes), collapse=", ")))
    message(sprintf("First entries of core dictionary: %s", paste(head(names(unlist(core_dictionary))), collapse=", ")))
    
    core_protein_list <- protein_list$protein[protein_list$protein %in% unlist(core_dictionary)]
    message(sprintf("Found %d highly expressed core proteins", length(core_protein_list)))
    
    # 4. DIAMOND: identify hits against essential proteins
    diamond_db_path <- file.path(output_dir, "essential_db.dmnd")
    system2("diamond", args = c("makedb", "--in", essential_fasta_path, "-d", diamond_db_path))
    
    diamond_output_path <- file.path(output_dir, "diamond_results.txt")
    system2("diamond", args = c("blastp", "-d", diamond_db_path,
                               "-q", proteome_path, "-o", diamond_output_path))
    
    # Process DIAMOND results
    diamond_results <- fread(diamond_output_path, header = FALSE, sep = "\t")
    colnames(diamond_results) <- c("query", "subject", "identity", "length", "mismatch", 
                                 "gapopen", "q.start", "q.end", "s.start", "s.end", 
                                 "evalue", "bitscore")
    
    # Filter significant DIAMOND hits (>40% identity, bitscore >50)
    diamond_hits <- diamond_results %>%
        filter(identity > 40, bitscore > 50) %>%
        pull(query) %>%
        unique()
    
    # 5. Identify CEPs: core proteins with hits in DIAMOND
    ceps_ids <- intersect(core_protein_list, diamond_hits)
    ceps_sequences <- filter_protein_sequences(ceps_ids, proteome_path)
    
    # 6. Identify potential CNAPs: core proteins without annotation
    uniprot_annotations <- map_to_uniprot_ids(core_protein_list)
    
    # Identify CNAPs based on the Protein.names field (uncharacterized proteins)
    cnaps_ids <- character(0)
    if (nrow(uniprot_annotations) > 0) {
        cnaps_ids <- uniprot_annotations$From[
            grepl("uncharacterized protein", tolower(uniprot_annotations$Protein.names), fixed = TRUE)
        ]
    }
    
    # Get CNAP sequences
    cnaps_sequences <- filter_protein_sequences(cnaps_ids, proteome_path)
    
    # 7. Save results
    writeXStringSet(ceps_sequences, filepath = file.path(output_dir, "CEPs.fasta"))
    writeXStringSet(cnaps_sequences, filepath = file.path(output_dir, "CNAPs.fasta"))
    
    write.table(ceps_ids, file = file.path(output_dir, "CEPs_ids.txt"),
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(cnaps_ids, file = file.path(output_dir, "CNAPs_ids.txt"),
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # Save UniProt annotations for both groups
    if (nrow(uniprot_annotations) > 0) {
        # CEPs annotations
        ceps_annotations <- uniprot_annotations[uniprot_annotations$From %in% ceps_ids, ]
        write.table(ceps_annotations,
                   file = file.path(output_dir, "CEPs_uniprot_annotations.tsv"),
                   sep = "\t",
                   row.names = FALSE,
                   quote = FALSE)
        
        # CNAPs annotations
        cnaps_annotations <- uniprot_annotations[uniprot_annotations$From %in% cnaps_ids, ]
        write.table(cnaps_annotations,
                   file = file.path(output_dir, "CNAPs_uniprot_annotations.tsv"),
                   sep = "\t",
                   row.names = FALSE,
                   quote = FALSE)
    }
    
    # Create summary
    summary_df <- data.frame(
        category = c("Total highly expressed core proteins", "CEPs", "CNAPs"),
        count = c(length(core_protein_list),
                 length(ceps_ids),
                 length(cnaps_ids))
    )
    
    write.table(summary_df, file = file.path(output_dir, "analysis_summary.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    return(list(
        ceps_sequences = ceps_sequences,
        ceps_ids = ceps_ids,
        cnaps_sequences = cnaps_sequences,
        cnaps_ids = cnaps_ids,
        summary = summary_df
    ))
}

#' Main function to run the entire analysis pipeline
#'
#' This function processes command line arguments and executes the analysis
#' with robust error handling
main <- function() {
    tryCatch({
        # Collect arguments
        args <- commandArgs(trailingOnly = TRUE)
        
        if (length(args) < 7) {
            stop("Insufficient number of arguments. Please provide: <genome_gtf_dir> <reads_dir> <threshold> <essential_fasta_path> <orthofinder_tsv_path> <output_dir> <proteome_path>")
        }
        
        # Assign arguments to variables with path normalization
        genome_gtf_dir <- normalizePath(args[1], mustWork = TRUE)
        reads_dir <- normalizePath(args[2], mustWork = TRUE)
        threshold <- as.numeric(args[3])
        essential_fasta_path <- normalizePath(args[4], mustWork = TRUE)
        orthofinder_tsv_path <- normalizePath(args[5], mustWork = TRUE)
        output_dir <- args[6]
        proteome_path <- normalizePath(args[7], mustWork = TRUE)
        
        # Create output directory if it doesn't exist
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        
        # Run analysis
        message("Starting analysis...")
        results <- get_up_CEPs(genome_gtf_dir, reads_dir, 
                             threshold, essential_fasta_path, orthofinder_tsv_path, 
                             output_dir, proteome_path)
        
        message("Analysis completed successfully")
        print(results$summary)
        
    }, error = function(e) {
        message(sprintf("Error in execution: %s", e$message))
        quit(status = 1)
    })
}

# Call the main function
main()
