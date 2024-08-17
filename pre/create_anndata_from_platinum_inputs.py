# This script builds and saves an AnnData object from a set of standard inputs
# It is meant to be run in situations where it takes a long time to load in the mtx file

# Note that this can also be accomplished usually by just calling sc.read_10x_mtx

# Required inputs
# matrix.mtx.gz - raw counts in sparse MEX format
# barcodes.tsv.gz - list of barcodes corresponding to column indices of the matrix
# features.tsv.gz - list of features corresponding to row indices of the matrix
# metadata.csv.gz - csv file with rows corresponding to barcodes.tsv.gz and any covariates of interest for this data


import argparse
import os
import sys
import loompy as lp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

sc.settings.verbosity = 3


def main(args):

    # Define files
    MTX_FILE = args.mtx_file
    OBS_FILE = args.obs_file
    VAR_FILE = args.var_file
    OUT_DIR = args.out_dir

    # Load in matrix
    print(f"Reading in {MTX_FILE} using scanpy read_mtx...")
    adata = sc.read_mtx(MTX_FILE)

    # Read in the barcode information
    print(f"Reading in {OBS_FILE}...")
    obs_df = pd.read_csv(OBS_FILE, delimiter="\t", index_col=0)
    obs_df.head()
    adata.obs = obs_df

    # Add feature data
    print(f"Reading in {VAR_FILE}...")
    var_df = pd.read_csv(VAR_FILE, delimiter="\t", header=None)
    var_df[["gene_name", "gene_id"]] = [row for row in var_df[0].str.split(":")]
    var_df = var_df.drop(0, axis=1).set_index("gene_name")
    var_df.head()
    adata.var = var_df
    adata.var_names_make_unique()

    # Prelim filtters, takes a bout a minute
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create ScanPy compatible h5ad and loom files from mtx, obs, and var files")
    parser.add_argument("--mtx_file", help="matrix market file to read in", required=True)
    parser.add_argument("--obs_file", help="tsv file with cell barcodes and metadata", required=True)
    parser.add_argument("--var_file", help="tsv file with gene names and metadata", required=True)
    parser.add_argument("--out_dir", help="directory to save output files", required=True)
    args = parser.parse_args()
    main(args)
    