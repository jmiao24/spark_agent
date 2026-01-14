"""MCP tools for SPARK - Spatial Pattern Recognition via Kernels

SPARK detects genes with spatial expression patterns in spatially resolved
transcriptomic data using generalized linear spatial models.
"""

import subprocess
import tempfile
from pathlib import Path
from typing import Annotated, Optional
import pandas as pd
from mcp.server.fastmcp import FastMCP

mcp = FastMCP("spark-example")

# Point to the R scripts directory for this tutorial
R_SCRIPT_DIR = Path(__file__).parent.parent / "r_scripts" / "02_spark_example"


@mcp.tool()
def create_spark_object(
    counts_csv: Annotated[str,
        "Path to CSV file with expression count matrix (genes × spots). "
        "First column should be gene names, remaining columns are spot IDs. "
        "Each cell contains the raw count for that gene in that spot."],
    location_csv: Annotated[str,
        "Path to CSV file with spatial coordinates (spots × 2). "
        "Two columns: x and y coordinates. Row order must match column order in counts_csv. "
        "No header required, just numeric coordinates."],
    percentage: Annotated[float,
        "Gene filtering threshold: retain genes expressed in at least this fraction of spots. "
        "Range: 0.0 to 1.0. Default 0.1 means keep genes expressed in ≥10% of spots."] = 0.1,
    min_total_counts: Annotated[int,
        "Minimum total counts per spot to retain. Spots with fewer total counts are filtered out. "
        "Typical value: 10-100 depending on data sparsity."] = 10,
    seed: Annotated[int,
        "Random seed for reproducibility. Use same seed to get identical results."] = 42,
) -> dict:
    """Create SPARK object from count matrix and spatial coordinates.

    This initializes a SPARK analysis by filtering lowly expressed genes
    and low-quality spots based on expression thresholds. The SPARK object
    contains the filtered data and calculated library sizes for normalization.
    """
    with tempfile.NamedTemporaryFile(suffix=".rds", delete=False) as f:
        output_rds = f.name

    subprocess.run([
        "Rscript", str(R_SCRIPT_DIR / "create_spark_object.R"),
        "--counts", counts_csv,
        "--location", location_csv,
        "--percentage", str(percentage),
        "--min_total_counts", str(min_total_counts),
        "--seed", str(seed),
        "--output", output_rds
    ], check=True)

    summary_csv = output_rds.replace(".rds", "_summary.csv")
    summary_df = pd.read_csv(summary_csv)
    Path(summary_csv).unlink()

    return {
        "message": f"SPARK object created with {summary_df['n_genes'][0]} genes and {summary_df['n_spots'][0]} spots",
        "reference": "https://github.com/xzhoulab/SPARK/blob/master/docs/pages/02_SPARK_Example.md",
        "spark_object_path": output_rds,
        "n_genes": int(summary_df["n_genes"][0]),
        "n_spots": int(summary_df["n_spots"][0]),
        "total_counts": int(summary_df["total_counts"][0]),
    }


@mcp.tool()
def spark_vc(
    spark_object_rds: Annotated[str,
        "Path to SPARK object RDS file (from create_spark_object output). "
        "This object contains the filtered count matrix and spatial coordinates."],
    covariates_csv: Annotated[Optional[str],
        "Path to CSV file with covariates matrix (spots × covariates). "
        "Each row corresponds to a spot, columns are covariates (e.g., batch effects, cell types). "
        "Leave empty (None) to fit model without covariates."] = None,
    num_core: Annotated[int,
        "Number of CPU cores for parallel processing. Higher values speed up computation "
        "for datasets with many genes. Typical range: 1-8."] = 1,
    verbose: Annotated[bool,
        "Print detailed progress messages during model fitting. "
        "Set to True for debugging or monitoring long-running jobs."] = False,
    seed: Annotated[int,
        "Random seed for reproducibility. Use same seed to get identical results."] = 42,
) -> dict:
    """Estimate SPARK model parameters under null hypothesis.

    This function fits a count-based spatial model under the null hypothesis
    of no spatial pattern. It estimates variance components and prepares the
    model for hypothesis testing in the next step (spark_test).
    """
    with tempfile.NamedTemporaryFile(suffix=".rds", delete=False) as f:
        output_rds = f.name

    cmd = [
        "Rscript", str(R_SCRIPT_DIR / "spark_vc.R"),
        "--spark_object", spark_object_rds,
        "--num_core", str(num_core),
        "--verbose", str(verbose).upper(),
        "--seed", str(seed),
        "--output", output_rds
    ]

    if covariates_csv:
        cmd.extend(["--covariates", covariates_csv])

    subprocess.run(cmd, check=True)

    summary_csv = output_rds.replace(".rds", "_summary.csv")
    summary_df = pd.read_csv(summary_csv)
    Path(summary_csv).unlink()

    return {
        "message": f"Model parameters estimated for {summary_df['n_genes_fitted'][0]} genes",
        "reference": "https://github.com/xzhoulab/SPARK/blob/master/docs/pages/02_SPARK_Example.md",
        "fitted_spark_object_path": output_rds,
        "n_genes_fitted": int(summary_df["n_genes_fitted"][0]),
        "status": str(summary_df["status"][0]),
    }


@mcp.tool()
def spark_test(
    fitted_spark_object_rds: Annotated[str,
        "Path to fitted SPARK object RDS file (from spark_vc output). "
        "This object contains the estimated model parameters under null hypothesis."],
    check_positive: Annotated[bool,
        "Validate that kernel matrices are positive definite. "
        "Set to True for robust analysis (recommended). False may speed up testing but risks numerical issues."] = True,
    verbose: Annotated[bool,
        "Print detailed progress messages during hypothesis testing. "
        "Shows progress for each gene and kernel tested."] = False,
    seed: Annotated[int,
        "Random seed for reproducibility. Use same seed to get identical results."] = 42,
) -> dict:
    """Test genes for spatial expression patterns.

    This function performs hypothesis testing to identify genes with significant
    spatial patterns. It tests multiple kernel functions (Gaussian, Periodic)
    and combines results. Returns p-values for each gene.
    """
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
        output_csv = f.name

    subprocess.run([
        "Rscript", str(R_SCRIPT_DIR / "spark_test.R"),
        "--spark_object", fitted_spark_object_rds,
        "--check_positive", str(check_positive).upper(),
        "--verbose", str(verbose).upper(),
        "--seed", str(seed),
        "--output", output_csv
    ], check=True)

    results_df = pd.read_csv(output_csv)
    tested_spark_rds = output_csv.replace(".csv", ".rds")
    Path(output_csv).unlink()

    # Count significant genes (adjusted p-value < 0.05)
    n_significant = (results_df["adjusted_pvalue"] < 0.05).sum()

    return {
        "message": f"Spatial pattern testing completed. Found {n_significant} significant genes (FDR < 0.05).",
        "reference": "https://github.com/xzhoulab/SPARK/blob/master/docs/pages/02_SPARK_Example.md",
        "tested_spark_object_path": tested_spark_rds,
        "n_genes_tested": len(results_df),
        "n_significant_genes": int(n_significant),
        "results_preview": results_df.nsmallest(10, "adjusted_pvalue")[
            ["gene", "combined_pvalue", "adjusted_pvalue"]
        ].to_dict(orient="records"),
    }


if __name__ == "__main__":
    mcp.run()
