#!/usr/bin/env bash

cut -f9 gencode.v25.primary_assembly.annotation.gtf | \
  awk -F '[ ;]' '{for(i=1;i<=NF;i++) if($i == "gene_type") print $(i+2)}' | \
  sort | uniq -c | wc -l
