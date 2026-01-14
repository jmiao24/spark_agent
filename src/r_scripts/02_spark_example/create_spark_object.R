#!/usr/bin/env Rscript
library(optparse)
library(data.table)
library(SPARK)

option_list <- list(
  make_option("--counts", type = "character",
              help = "Count matrix CSV file (genes × spots). Rows are genes, columns are spots/cells."),
  make_option("--location", type = "character",
              help = "Location CSV file with x,y coordinates. Two columns (x, y), rows match count matrix columns."),
  make_option("--percentage", type = "double", default = 0.1,
              help = "Gene filtering threshold: keep genes expressed in at least this fraction of spots [default: %default]"),
  make_option("--min_total_counts", type = "integer", default = 10,
              help = "Minimum total counts per spot to keep [default: %default]"),
  make_option("--seed", type = "integer", default = 42,
              help = "Random seed for reproducibility [default: %default]"),
  make_option("--output", type = "character",
              help = "Output RDS file for SPARK object (required)")
)

args <- parse_args(OptionParser(option_list = option_list))

# Set random seed for reproducibility
set.seed(args$seed)

# Read input data
counts_dt <- fread(args$counts)
counts <- as.matrix(counts_dt, rownames = 1)

location_dt <- fread(args$location)
location <- as.data.frame(location_dt[, -1])  # Remove first column (row names)
rownames(location) <- location_dt[[1]]  # Set row names from first column

# Create SPARK object
spark_obj <- CreateSPARKObject(
  counts = counts,
  location = location,
  percentage = args$percentage,
  min_total_counts = args$min_total_counts
)

# Calculate library size
spark_obj@lib_size <- apply(spark_obj@counts, 2, sum)

# Save SPARK object
saveRDS(spark_obj, args$output)

# Output summary statistics
result <- data.table(
  n_genes = nrow(spark_obj@counts),
  n_spots = ncol(spark_obj@counts),
  total_counts = sum(spark_obj@counts)
)
fwrite(result, sub("\\.rds$", "_summary.csv", args$output))
