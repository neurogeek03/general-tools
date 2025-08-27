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
