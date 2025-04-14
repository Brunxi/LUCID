# Load necessary libraries
library(dplyr)       # For data manipulation
library(data.table)  # For fast data import and manipulation
library(DESeq2)      # For differential gene expression analysis
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
    total_proteins <- length(protein_ids)
    results <- data.frame()
    
    # Process in batches to avoid API limitations
    for (i in seq(1, total_proteins, by = batch_size)) {
        batch_ids <- protein_ids[i:min(i + batch_size - 1, total_proteins)]
        
        # Prepare the request body
        files <- list(
            from = "Ensembl_Genomes_Protein",
            to = "UniProtKB",
            ids = paste(batch_ids, collapse = ",")
        )
        
        # Use error handling to prevent script failure on API errors
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
    }
    
    return(results)
}

#' Process gene expression data from RNA-seq experiments
#'
#' This function aligns reads, counts features, and performs differential expression analysis
#'
#' @param genome_gtf_dir Directory containing genome GTF and FASTA files
#' @param reads_dir Directory containing FASTQ files
#' @param output_dir Directory for output files
#' @param threshold Log2FoldChange threshold for significant genes
#' @return Vector of significantly upregulated gene names
process_gene_expression <- function(genome_gtf_dir, reads_dir, output_dir, threshold) {
  
    # Define the full index name with the path to the genome directory
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
  
    # Create a vector of clean sample names
    sample_names <- gsub("_1\\.fastq\\.subread\\.BAM$", "", basename(bam_files))
  
    # Count reads per feature
    fc <- featureCounts(files = bam_files, 
                      annot.ext = file.path(genome_gtf_dir, "genome.gtf"), 
                      isGTFAnnotationFile = TRUE, 
                      isPairedEnd = TRUE, 
                      useMetaFeatures = TRUE, 
                      GTF.featureType = "exon",
                      GTF.attrType = "gene_id")
  
    # Replace column names with clean sample names
    colnames(fc$counts) <- sample_names
  
    # Write the count matrix
    write.table(fc$counts, file = file.path(output_dir, "count_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)
  
    # Create coldata based on sample names
    coldata_path <- file.path(genome_gtf_dir, "coldata.csv")
  
    # Read the completed coldata file
    coldata <- read.csv(coldata_path, row.names = 1)
  
    # Verify there are no NAs in the 'treatment' column
    if (any(is.na(coldata$treatment))) {
        stop("Could not determine treatment for one or more samples. Check values in coldata.csv.")
    }
  
    # Verify that names in coldata match sample_names
    if (!all(sample_names %in% rownames(coldata))) {
        stop("Sample names do not match names in coldata.csv. Please verify.")
    }
  
    # Ensure coldata is in the same order as sample names
    coldata <- coldata[sample_names, , drop = FALSE]
  
    # Differential expression analysis using DESeq2
    dds <- DESeqDataSetFromMatrix(countData = fc$counts,
                                colData = coldata,
                                design = ~ treatment)
    dds <- DESeq(dds)
  
    # Data transformation and PCA
    vsdata <- varianceStabilizingTransformation(dds, blind = FALSE)
    pca_plot <- plotPCA(vsdata, intgroup = "treatment")
  
    # Save PCA plot to output directory
    pca_plot_path <- file.path(output_dir, "pca_plot.png")
    png(pca_plot_path, width = 1000, height = 800)
    print(pca_plot)
    dev.off()  
  
    # Get unique treatments excluding the control
    treatments <- unique(coldata$treatment[coldata$treatment != "ctrl"])
  
    # Generate results for each treatment compared to the control
    significant_genes_list <- list()
  
    for (treatment in treatments) {
        results_df <- results(dds, contrast = c("treatment", treatment, "ctrl"))
        results_df <- as.data.frame(results_df)
    
        # Filter genes with log2FoldChange greater than the threshold    
        significant_genes <- results_df[!is.na(results_df$log2FoldChange) & results_df$log2FoldChange > threshold, ]
    
        significant_genes_list[[treatment]] <- significant_genes
    
        # Save significant genes in TSV with gene names only
        significant_gene_names <- gsub("^gene:", "", rownames(significant_genes))
        write.table(significant_gene_names, 
                   file = file.path(output_dir, paste0("significant_genes_", treatment, "_vs_ctrl.tsv")), 
                   sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    }
  
    return(significant_gene_names)
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

#' Identify upregulated Core Effector Proteins (CEPs) from RNA-seq data
#'
#' This function processes RNA-seq data to identify upregulated core effector proteins
#' and conserved non-annotated proteins (CNAPs) in phytopathogenic fungi
#'
#' @param genome_gtf_dir Directory containing genome GTF and FASTA files
#' @param reads_dir Directory containing RNA-seq FASTQ files
#' @param threshold Log2FoldChange threshold for significant genes
#' @param essential_fasta_path Path to FASTA file with essential proteins
#' @param orthofinder_tsv_path Path to OrthoFinder TSV results
#' @param output_dir Directory for output files
#' @param proteome_path Path to proteome FASTA file
#' @return List with CEPs and CNAPs sequences, IDs, and summary
get_up_CEPs <- function(genome_gtf_dir, reads_dir, threshold, 
                       essential_fasta_path, orthofinder_tsv_path, output_dir, proteome_path) {
    # 1. Get upregulated genes from RNA-seq analysis
    significant_genes <- process_gene_expression(genome_gtf_dir, reads_dir, output_dir, 1.0)
    
    # 2. Create dictionary of core genes from OrthoFinder results
    core_dictionary <- create_core_gene_dict(orthofinder_tsv_path)
    
    # 3. Get upregulated core proteins
    proteome_dictionary <- create_protein_gene_dict(proteome_path)
    protein_list <- get_protein_ids(proteome_dictionary, significant_genes)
    core_protein_list <- protein_list$protein[protein_list$protein %in% unlist(core_dictionary)]
    
    # 4. DIAMOND: identify hits against essential proteins
    # Create DIAMOND database from essential proteins
    diamond_db_path <- file.path(output_dir, "essential_db.dmnd")
    system2("diamond", args = c("makedb", "--in", essential_fasta_path, "-d", diamond_db_path))
    
    # Run DIAMOND blastp
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
    
    # 6. Get UniProt annotations for all core proteins
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
    # Save FASTA sequences
    writeXStringSet(ceps_sequences, filepath = file.path(output_dir, "CEPs.fasta"))
    writeXStringSet(cnaps_sequences, filepath = file.path(output_dir, "CNAPs.fasta"))
    
    # Save IDs
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
        category = c("Total upregulated core proteins", "CEPs", "CNAPs"),
        count = c(length(core_protein_list),
                 length(ceps_ids),
                 length(cnaps_ids))
    )
    
    write.table(summary_df, file = file.path(output_dir, "analysis_summary.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Return results as a list
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
main <- function() {
    # Collect arguments passed from the command line
    args <- commandArgs(trailingOnly = TRUE)
  
    # Verify if the correct number of arguments has been passed
    if (length(args) < 7) {
        stop("Insufficient number of arguments. Provide: <genome_gtf_dir> <reads_dir> <threshold> <essential_fasta_path> <orthofinder_tsv_path> <output_dir> <proteome_path>")
    }
  
    # Assign arguments to variables
    genome_gtf_dir <- args[1]
    reads_dir <- args[2]
    threshold <- as.numeric(args[3])
    essential_fasta_path <- args[4]
    orthofinder_tsv_path <- args[5]
    output_dir <- args[6]
    proteome_path <- args[7]
  
    # Execute the desired function with the arguments
    results <- get_up_CEPs(genome_gtf_dir, reads_dir, threshold, essential_fasta_path, orthofinder_tsv_path, output_dir, proteome_path)
  
    # Print the final results to be seen in the console
    print(results)
}

# Call the main function
main()
