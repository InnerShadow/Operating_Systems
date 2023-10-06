#!/usr/bin/env bash

# Извлечение кодирующих интервалов для протеин-кодирующих генов и сохранение в unmerged_CDS.bed
grep 'gene_type "protein_coding"' gencode.v25.primary_assembly.annotation.gtf | \
    awk -F '\t' '{if ($3 == "CDS") print $1, $4-1, $5}' OFS='\t' > unmerged_CDS.bed


