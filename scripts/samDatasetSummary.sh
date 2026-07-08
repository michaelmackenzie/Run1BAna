#!/bin/bash
# usage
if [[ -z "$1" ]]; then
  echo "Usage: $0 <dataset>" >&2
  exit 1
fi

# Get the dataset name
dataset="$1"
# Obtain summary once
summary_txt=$(samweb list-files --summary "dh.dataset $dataset and availability:anylocation" 2>/dev/null)
nfiles=$(echo "$summary_txt" | awk '/File count:/ {print $3}')
size=$(echo "$summary_txt" | awk '/Total size:/ {print $3}')
triggered=$(echo "$summary_txt" | awk '/Event count:/ {print $3}')

# Ensure numeric defaults
nfiles=${nfiles:-0}
size=${size:-0}
triggered=${triggered:-0}

# Exit if no files found
if (( nfiles == 0 )); then
    echo "Error: No files found in dataset" >&2
    exit 1
fi

# Calculate total generated events from all files
echo "Calculating the N(generated) information"
generated=$(samweb list-files "dh.dataset $dataset and availability: anylocation" 2>/dev/null | \
  xargs -n1 samweb get-metadata | awk '/dh.gencount/ { sum += $2;} END { print sum+0 }')

printf "N(events): %s\nN(generated): %s\nN(files): %s\nDataset size: %s\n" "$triggered" "$generated" "$nfiles" "$size"
