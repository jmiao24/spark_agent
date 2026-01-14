#!/usr/bin/env Rscript
library(optparse)
library(data.table)
library(SPARK)

option_list <- list(
  make_option("--spark_object", type = "character",
              help = "Path to fitted SPARK object RDS file (from spark_vc output)"),
  make_option("--check_positive", type = "logical", default = TRUE,
              help = "Check if kernel matrix is positive definite [default: %default]"),
  make_option("--verbose", type = "logical", default = FALSE,
              help = "Print detailed progress messages [default: %default]"),
  make_option("--seed", type = "integer", default = 42,
              help = "Random seed for reproducibility [default: %default]"),
  make_option("--output", type = "character",
              help = "Output CSV file for test results (required)")
)

args <- parse_args(OptionParser(option_list = option_list))

# Set random seed for reproducibility
set.seed(args$seed)

# Load fitted SPARK object
spark_obj <- readRDS(args$spark_object)

# Test for spatially expressed genes
spark_obj <- spark.test(
  spark_obj,
  check_positive = args$check_positive,
  verbose = args$verbose
)

# Extract results
results <- as.data.table(spark_obj@res_mtest, keep.rownames = "gene")

# Save results
fwrite(results, args$output)

# Also save updated SPARK object with test results
saveRDS(spark_obj, sub("\\.csv$", ".rds", args$output))
