# Fix multiprocessing deadlock in asyncio context
# Must be called before any other imports
import multiprocessing
multiprocessing.set_start_method("spawn", force=True)

"""
Model Context Protocol (MCP) for SPARK

SPARK (Spatial Pattern Recognition via Kernels) detects genes with spatial expression
patterns in spatially resolved transcriptomic data using generalized linear spatial models.
This package enables identification of spatially variable genes in tissue samples.

This MCP Server provides Python interfaces to R tools extracted from the following tutorial files:
1. 02_spark_example (Example Analysis with SPARK: Breast Cancer Data)
    - create_spark_object: Create SPARK object from count matrix and spatial coordinates (calls R via Rscript)
    - spark_vc: Estimate SPARK model parameters under null hypothesis (calls R via Rscript)
    - spark_test: Test genes for spatial expression patterns (calls R via Rscript)

Note: All tools execute R code via Rscript subprocess calls. Ensure R is installed
and the package dependencies are available in the renv environment at repo/SPARK/.
"""

from mcp.server.fastmcp import FastMCP
import importlib

# Import statements (using importlib for numeric-prefixed modules)
spark_example_02 = importlib.import_module('tools.02_spark_example')

# Server definition
mcp = FastMCP(name="SPARK")

# Register tools from 02_spark_example module
create_spark_object = spark_example_02.create_spark_object
spark_vc = spark_example_02.spark_vc
spark_test = spark_example_02.spark_test

mcp.add_tool(create_spark_object)
mcp.add_tool(spark_vc)
mcp.add_tool(spark_test)

if __name__ == "__main__":
    mcp.run()
