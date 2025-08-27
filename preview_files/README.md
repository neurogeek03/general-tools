## Preview Seurat / Anndata objects ( `.rds` or `.h5ad` ) & large `.csv` files

<aside>

> **Features:**
> 
> - Requests an interactive compute node (`n`) for each file
> - Activates the correct conda environment based on file type & runs the appropriate preview script:
> `.h5ad` → anndata_env → `preview_anndata.py`
> `.rds` → seurat_env → `preview_rds.R`
> `.csv` → anndata_env → prints first 10 rows with all columns using pandas
> 
> **Usage:**
> `preview <file1> [file2 ... fileN]`
> 
</aside>

I wrote 2 very simple scripts (one for `.rds` and one for `.h5ad` files) that can be used to preview objects. I have pasted them below: 

### R script

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

# Get the first command-line argument (the .rds file path)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Please provide a .rds file as the first argument.", call. = FALSE)
}

file_path <- args[1]

# Load the Seurat object
seurat_obj <- readRDS(file_path)

cat(sprintf("\nLoaded %s\n", file_path))
cat(sprintf("Object class: %s\n", class(seurat_obj)))

# Print summary info
print(seurat_obj)

# Preview @meta.data if available
if ("meta.data" %in% slotNames(seurat_obj)) {
  cat("\nPreview of @meta.data (first 5 rows):\n")
  print(head(seurat_obj@meta.data, 5))
} else {
  cat("No @meta.data slot found in this object.\n")
}

# Preview gene names
if ("assays" %in% slotNames(seurat_obj)) {
  default_assay <- DefaultAssay(seurat_obj)
  cat(sprintf("\nDefault assay: %s\n", default_assay))
  
  gene_names <- rownames(seurat_obj[[default_assay]])
  
  cat("\nGene names (first 10):\n")
  print(head(gene_names, 10))
} else {
  cat("No assays found in this object.\n")
}
```

### Python script

```python
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

```

Normally, one can execute these by: 

1. Requesting a compute node `n` (Instructions in my setup tutorial)
2. Activating the right environment (personally, I have 2 separate environments, one for each object type)
3. Finding the path of the script and running something like this to execute either script:
    
    ```bash
    <python or R> <path/to/preview/script> <path/to/file/of/choice>
    ```
    

But all this is a lot of work each time you just want to preview a file, so I set up the following: 

### 1. Place your preview scripts in 1 folder and name it `preview`

### 2. Move that folder to your home dir

```bash
mv ./preview ~
```

### 3. Create a binaries folder in your home dir

```bash
cd ~
mkdir -p ~/bin
```

### 4. Create a binary file called `preview` within the bin folder

Then paste this & save:

```bash
#!/bin/bash -i
# Usage: preview <file1> [file2 ... fileN]
# Works for .h5ad, .rds, and .csv files

if [ $# -lt 1 ]; then
    echo "Usage: preview <file1> [file2 ... fileN]"
    exit 1
fi

for FILE in "$@"; do
    EXT="${FILE##*.}"
    echo "Previewing $FILE (.$EXT)..."

    n <<EOF
    source ~/.bashrc

    case "$EXT" in
        h5ad)
            conda activate anndata_env
            python ~/preview/preview_anndata.py "$FILE"
            ;;
        rds)
            conda activate seurat_env
            Rscript ~/preview/preview_rds.R "$FILE"
            ;;
        csv)
            conda activate anndata_env
            python - <<PYTHON_EOF
import pandas as pd
df = pd.read_csv("$FILE")
with pd.option_context('display.max_columns', None):
    print(df.head(10))
PYTHON_EOF
            ;;
        *)
            echo "Unsupported file type: $EXT"
            ;;
    esac
EOF
done
```

### 5. Make it executable

```bash
chmod +x ~/bin/preview
```

### 6. Add ~/bin to $PATH

We add `~/bin` to your `$PATH` so the shell knows to look there when you type a command.

Here’s the idea:

- When you type a command like `ls` or `python`, the shell searches a list of directories (stored in `$PATH`) to find the executable.
- By default, your home folder `~/bin` isn’t usually in `$PATH`, so if you put a script there, the shell won’t find it.
- By adding it:

```bash
echo 'export PATH=$HOME/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

you tell the shell: “Look in my personal `~/bin` folder first when searching for commands.”

You can verify if you added it like so:

```bash
echo $PATH | tr ':' '\n'
```

### 7. Test it!

From anywhere in the HPC, you should be able to preview any `.rds` or `.h5ad` file, even if it’s large.

```bash
preview <path/to/your/file>
```