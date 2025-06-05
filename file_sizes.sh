#!/bin/bash

# Check if base directory was passed
if [ -z "$1" ]; then
  echo "Usage: $0 <base_directory>"
  exit 1
fi

BASE_DIR="$1"

# Loop over all immediate subdirectories
for dir in "$BASE_DIR"/*/; do
  [ -d "$dir" ] || continue  # Skip if not a directory
  echo "===== Directory: $dir ====="
  find "$dir" -type f -print0 | while IFS= read -r -d '' file; do
    size=$(du -h "$file" | cut -f1)
    echo "$file: $size"
  done
  echo ""
done