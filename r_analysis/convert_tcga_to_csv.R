# R Script: convert_tcga_to_csv.R
# Author: Gemini Assistant
# Description: This script reads a complex, multi-header TCGA Level 3 RSEM gene expression file (.txt),
#              extracts only the raw counts, cleans up gene and sample identifiers,
#              and saves the result as a clean, human-readable CSV file.

# --- Configuration ---
# Input file should be in the same directory as the script
input_file <- "OV.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt"
# Output file will be created in the same directory
output_file <- "TCGA-OV-counts_cleaned.csv"

cat("Starting TCGA data conversion process...\n")
cat("Input file:", input_file, "\n")

# --- 1. Parse Header to get Sample IDs ---
cat("Step 1: Parsing header to extract sample IDs...\n")
header_line <- read.delim(input_file, nrows = 1, header = FALSE, stringsAsFactors = FALSE)
# The sample IDs are in columns starting from the 2nd, and repeat every 3 columns
sample_ids_raw <- as.character(header_line[1, seq(2, ncol(header_line), by = 3)])
# Per user request, use the raw, unmodified sample IDs
sample_ids <- sample_ids_raw


# --- 2. Read and Extract Raw Count Data ---
cat("Step 2: Reading main data table and extracting raw_count columns...\n")
# Read the main data, skipping the first (header) line. The first column will be treated as row names.
full_data <- read.delim(input_file, skip = 1, header = TRUE, stringsAsFactors = FALSE, row.names = 1)

# The 'raw_count' columns are the 1st, 4th, 7th, etc.
raw_count_cols <- seq(from = 1, to = ncol(full_data), by = 3)
counts_matrix <- full_data[, raw_count_cols]

# Assign the raw sample IDs as column names
colnames(counts_matrix) <- sample_ids

# --- 3. Handle Duplicates (using original gene IDs) ---
cat("Step 3: Handling potential duplicate gene IDs...\n")

# Use the original, unmodified gene IDs from the row names
gene_ids_raw <- rownames(counts_matrix)

# Remove genes with unknown symbols ('?')
unknown_genes_mask <- grepl("^\\?", gene_ids_raw)
counts_matrix <- counts_matrix[!unknown_genes_mask, ]
gene_ids_raw <- gene_ids_raw[!unknown_genes_mask]
cat(sum(unknown_genes_mask), "genes with unknown symbols ('?') were removed.\n")


# Handle duplicate gene symbols by taking the mean expression
if (any(duplicated(gene_ids_raw))) {
    cat("Found duplicated gene IDs. Merging them by taking the mean expression...\n")
    # Use the aggregate function to group by the raw gene ID and calculate the mean for each row
    counts_matrix_agg <- aggregate(counts_matrix, by = list(gene_id = gene_ids_raw), FUN = mean)
    # The gene IDs become the first column after aggregation, set them as row names
    rownames(counts_matrix_agg) <- counts_matrix_agg$gene_id
    # Remove the temporary gene_id column
    counts_matrix <- counts_matrix_agg[, -1]
} else {
    # If no duplicates, just set the raw gene IDs as row names
    rownames(counts_matrix) <- gene_ids_raw
}

# --- 4. Final Touches and Save to CSV ---
cat("Step 4: Rounding values and saving to CSV...\n")
# Ensure all values are integers, as DESeq2 prefers
counts_matrix_final <- round(counts_matrix)

# Write the clean data frame to a CSV file
write.csv(counts_matrix_final, file = output_file, row.names = TRUE)

cat("\nConversion complete!\n")
cat("Clean data has been saved to:", output_file, "\n") 