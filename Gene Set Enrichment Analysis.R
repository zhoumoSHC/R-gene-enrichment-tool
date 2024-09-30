#!/usr/bin/env Rscript

# Title: Gene Set Enrichment Analysis and Pathway Annotation
# Author: Æ¤Ë¿, improved by ¿ËÀÏÊ¦
# Date: 2024-09-30
# Description: This script performs gene set enrichment analysis and annotates genes with their associated pathways.

# ===== USER CONFIGURATION =====
# Set these parameters before running the script

# Base directory for input/output files
base_dir <- ''

# Path to the gene list file
genelist_file <- 'genelist'

# Path to the GMT file (download from https://www.gsea-msigdb.org/gsea/downloads.jsp )
gmt_file <- 'h.all.v2023.2.Hs.symbols.gmt'

# Column number in the gene list file that contains the gene symbols
symbol_loc <- 17

# ===== END OF USER CONFIGURATION =====

log_message <- function(message) {
  cat(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - ", message, "\n"))
}

read_and_process_pathways <- function(file_path) {
  log_message("Reading and processing pathways")
  pathways <- scan(file_path, sep = '\n', what = "character", quiet = TRUE)
  pathways <- strsplit(pathways, "\t")
  names(pathways) <- sapply(pathways, `[[`, 1)
  lapply(pathways, `[`, -c(1, 2))
}

perform_enrichment <- function(gene_list, pathways, hallmark_genes, symbol_loc) {
  log_message("Performing enrichment analysis")
  result_df <- data.frame(stringsAsFactors = FALSE)
  
  for (i in seq_along(pathways)) {
    pathway_genes <- unlist(pathways[i])
    overlap_genes <- gene_list[gene_list[, symbol_loc] %in% pathway_genes, symbol_loc]
    
    count <- length(overlap_genes)
    pophit <- length(pathway_genes)
    poptotal <- length(hallmark_genes)
    listtotal <- sum(gene_list[, symbol_loc] %in% hallmark_genes)
    
    FE <- (count / listtotal) / (pophit / poptotal)
    term <- names(pathways)[i]
    
    contingency_table <- matrix(c(count, listtotal - count, pophit - count, poptotal - listtotal - pophit + count), nrow = 2)
    fisher_test <- fisher.test(contingency_table, alternative = "greater")
    p_value <- fisher_test$p.value
    
    result_df <- rbind(result_df, data.frame(
      Category = term,
      term = term,
      Count = count,
      percent = count * 100 / listtotal,
      PValue = p_value,
      genes = paste(overlap_genes, collapse = ", "),
      listtotal = listtotal,
      pophit = pophit,
      poptotal = poptotal,
      FoldEnrichment = FE,
      stringsAsFactors = FALSE
    ))
  }
  
  result_df[order(result_df$PValue), ]
}

annotate_genes_with_pathways <- function(gene_list, pathways, symbol_loc) {
  log_message("Annotating genes with pathways")
  
  unique_genes <- unique(gene_list[, symbol_loc])
  unique_pathways <- names(pathways)
  
  annotation_matrix <- matrix(0, 
                              nrow = length(unique_genes), 
                              ncol = length(unique_pathways),
                              dimnames = list(unique_genes, unique_pathways))
  
  for (pathway_name in unique_pathways) {
    pathway_genes <- pathways[[pathway_name]]
    common_genes <- intersect(pathway_genes, unique_genes)
    if (length(common_genes) > 0) {
      annotation_matrix[common_genes, pathway_name] <- 1
    }
  }
  
  annotation_df <- as.data.frame(annotation_matrix)
  annotation_df$Gene <- rownames(annotation_df)
  
  annotation_df <- annotation_df[, c(ncol(annotation_df), 1:(ncol(annotation_df)-1))]
  
  annotation_df$Pathways_Count <- rowSums(annotation_df[, -1])
  annotation_df$Pathways_List <- apply(annotation_df[, -1], 1, function(x) {
    paste(names(x)[x == 1], collapse = "; ")
  })
  
  col_order <- c("Gene", "Pathways_Count", "Pathways_List", 
                 setdiff(colnames(annotation_df), c("Gene", "Pathways_Count", "Pathways_List")))
  annotation_df[, col_order]
}

main <- function() {
  log_message("Analysis started")
  
  # Set working directory
  setwd(base_dir)
  genelist_path <- file.path(base_dir, genelist_file)
  inflam.path <- read_and_process_pathways(gmt_file)
  
  # Read gene list
  log_message("Reading gene list")
  gene_list <- read.table(genelist_path, header = TRUE, stringsAsFactors = FALSE)
  gene_list <- gene_list[!duplicated(gene_list[, symbol_loc]), ]
  
  # Process Hallmark genes
  hallmark_genes <- unique(unlist(inflam.path))
  
  # Perform enrichment analysis
  result <- perform_enrichment(gene_list, inflam.path, hallmark_genes, symbol_loc)
  
  # Multiple testing correction
  log_message("Performing multiple testing correction")
  result$Bonferroni <- p.adjust(result$PValue, method = "bonferroni")
  result$Benjamini <- p.adjust(result$PValue, method = "BH")
  result$FDR <- p.adjust(result$PValue, method = "fdr")
  
  # Write enrichment analysis results
  output_file <- file.path(base_dir, paste0('enrichment_results_', format(Sys.time(), "%Y%m%d_%H%M%S"), '.tsv'))
  log_message(paste("Writing enrichment analysis results to", output_file))
  write.table(result, output_file, sep = '\t', quote = FALSE, row.names = FALSE)
  
  # Perform gene pathway annotation
  gene_pathway_annotations <- annotate_genes_with_pathways(gene_list, inflam.path, symbol_loc)
  
  # Write the gene pathway annotations to a file
  annotation_output_file <- file.path(base_dir, paste0("gene_pathway_annotations_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv"))
  log_message(paste("Writing gene pathway annotations to", annotation_output_file))
  write.table(gene_pathway_annotations, file = annotation_output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Log final statistics
  log_message(paste("Number of genes analyzed:", nrow(gene_list)))
  log_message(paste("Number of pathways analyzed:", length(inflam.path)))
  log_message("Analysis completed")
}

main()
