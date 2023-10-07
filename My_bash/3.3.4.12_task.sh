#!/usr/bin/env bash

grep -e $'protein_coding' gencode.v25.primary_assembly.annotation.gtf | awk -F'\t' '$3 == "exon"' > protein_coding_exons.gtf

awk 'BEGIN { max_diff = 0; gene_name = "" } { diff = $5 - $4 + 1; if (diff > max_diff) { max_diff = diff; gene_name = $18 } } END { print gene_name,max_diff }' protein_coding_exons.gtf

