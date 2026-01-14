# [Paper2Agent](https://github.com/jmiao24/Paper2Agent): SPARK Demo

A demonstration of turning the [SPARK paper](https://www.nature.com/articles/s41592-019-0701-7) into an interactive AI agent. This project transforms Sun et al.'s SPARK (Spatial Pattern Recognition via Kernels) method for identifying spatially variable genes in spatially resolved transcriptomic data into a conversational agent that can create SPARK objects, fit spatial models, and test for spatial expression patterns through natural language.

## Folder Structure

```
spark_agent/
├── src/
│   ├── SPARK_mcp.py              # MCP server entry point
│   ├── requirements.txt          # Python dependencies
│   ├── tools/                    # Python wrappers for R scripts
│   │   └── 02_spark_example.py   # SPARK tool implementations
│   └── r_scripts/                # R scripts for each tool
│       └── 02_spark_example/
│           ├── create_spark_object.R  # Initialize SPARK object
│           ├── spark_vc.R             # Estimate null model parameters
│           └── spark_test.R           # Test for spatial patterns
└── README.md
```

## Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/jmiao24/spark_agent.git
cd spark_agent
```

### 2. Install Gemini CLI

Install the [Google Gemini CLI](https://github.com/google-gemini/gemini-cli):

```bash
brew install gemini-cli
```

### 3. Install R Dependencies

The agent requires R and the `SPARK` package. Install R from [CRAN](https://cran.r-project.org/), then install the required packages:

```r
# Install dependencies
install.packages(c("optparse", "data.table"))

# Install SPARK from GitHub
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("xzhoulab/SPARK")
```

### 4. Install FastMCP

```bash
pip install fastmcp
```

### 5. Install MCP Server

Install the SPARK MCP server using fastmcp:

```bash
fastmcp install gemini-cli ./src/SPARK_mcp.py --with-requirements ./src/requirements.txt
```

### 6. Start the Agent

Start Gemini CLI in the repository folder:

```bash
gemini
```

You will now have access to the SPARK agent with all available tools.

## Example Query

```
Create a SPARK object from my count matrix and location data, then fit the spatial
model and test for spatially variable genes. Show me the top significant genes.
```

## Available Agent Tools

The agent provides the following capabilities through natural language:

### Data Initialization
- `create_spark_object`: Create a SPARK object from count matrix and spatial coordinates
  - Filters genes based on expression percentage threshold
  - Filters spots based on minimum total counts
  - Calculates library sizes for normalization

### Model Fitting
- `spark_vc`: Estimate SPARK model parameters under null hypothesis
  - Fits count-based spatial model
  - Estimates variance components
  - Supports covariates for batch effects or cell types
  - Supports parallel processing with multiple cores

### Hypothesis Testing
- `spark_test`: Test genes for spatial expression patterns
  - Tests multiple kernel functions (Gaussian, Periodic)
  - Computes p-values for each gene
  - Identifies spatially variable genes (FDR < 0.05)
  - Returns results preview with top significant genes

## Analysis Workflow

The typical SPARK analysis follows a three-stage pipeline:

1. **Data Initialization**: Prepare count matrix (genes × spots) and location data, then call `create_spark_object` to filter and initialize the analysis object

2. **Model Fitting**: Call `spark_vc` to estimate model parameters under the null hypothesis of no spatial pattern

3. **Hypothesis Testing**: Call `spark_test` to identify genes with significant spatial expression patterns

## About SPARK

SPARK (Spatial Pattern Recognition via Kernels) is a statistical method for identifying genes with spatial expression patterns in spatially resolved transcriptomic data. Key features include:

- Generalized linear spatial models for count-based expression data
- Multiple kernel functions to capture diverse spatial patterns
- Variance component estimation for spatial correlation
- Hypothesis testing with combined p-values across kernels
- Support for covariates to control for batch effects

For more details, see the [original paper](https://www.nature.com/articles/s41592-019-0701-7) and the [SPARK documentation](https://xzhoulab.github.io/SPARK/).
