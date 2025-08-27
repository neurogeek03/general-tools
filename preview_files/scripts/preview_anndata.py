import scanpy as sc
import argparse
import pandas as pd

def preview_obs_var(h5ad_file, n=5):
    # Load the AnnData object
    adata = sc.read_h5ad(h5ad_file)

    # Print file info
    print(f"\nLoaded {h5ad_file} with shape {adata.shape}")
    print(adata)

    # Preview .obs
    print(f"\nPreview of .obs (first {n} rows, all columns):")
    with pd.option_context('display.max_columns', None):
        print(adata.obs.head(n))

    # Preview .var
    print(f"\nPreview of .var (first {n} rows, all columns):")
    with pd.option_context('display.max_columns', None):
        print(adata.var.head(n))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preview the .obs, .var of an AnnData (.h5ad) file.")
    parser.add_argument("file", help="Path to the .h5ad file.")
    parser.add_argument("-n", type=int, default=5, help="Number of rows to preview (default: 5)")
    args = parser.parse_args()

    preview_obs_var(args.file, args.n)
