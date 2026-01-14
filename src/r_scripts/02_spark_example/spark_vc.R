#!/usr/bin/env Rscript
library(optparse)
library(data.table)
library(SPARK)

option_list <- list(
  make_option("--spark_object", type = "character",
              help = "Path to SPARK object RDS file (from create_spark_object output)"),
  make_option("--covariates", type = "character", default = "",
              help = "Optional CSV file with covariates matrix (spots × covariates). Leave empty for no covariates."),
  make_option("--num_core", type = "integer", default = 1,
              help = "Number of CPU cores for parallel processing [default: %default]"),
  make_option("--verbose", type = "logical", default = FALSE,
              help = "Print detailed progress messages [default: %default]"),
  make_option("--seed", type = "integer", default = 42,
              help = "Random seed for reproducibility [default: %default]"),
  make_option("--output", type = "character",
              help = "Output RDS file for fitted SPARK object (required)")
)

args <- parse_args(OptionParser(option_list = option_list))

# Set random seed for reproducibility
set.seed(args$seed)

# Load SPARK object
spark_obj <- readRDS(args$spark_object)

# Load covariates if provided
covariates <- NULL
if (args$covariates != "") {
  covariates <- as.matrix(fread(args$covariates))
}

# Estimate parameters under null hypothesis
spark_obj <- spark.vc(
  spark_obj,
  covariates = covariates,
  lib_size = spark_obj@lib_size,
  num_core = args$num_core,
  verbose = args$verbose
)

# Save fitted SPARK object
saveRDS(spark_obj, args$output)

# Output summary
result <- data.table(
  n_genes_fitted = length(spark_obj@res_vc),
  status = "Model parameters estimated under null hypothesis"
)
fwrite(result, sub("\\.rds$", "_summary.csv", args$output))
