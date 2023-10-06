#!/usr/bin/env bash

# Извлечение кодирующих интервалов для протеин-кодирующих генов и сохранение в unmerged_CDS.bed
grep 'gene_type "protein_coding"' gencode.v25.primary_assembly.annotation.gtf | \
    awk -F '\t' '{if ($3 == "CDS") print $1, $4-1, $5}' OFS='\t' > unmerged_CDS.bed

# Сортируем
sort -k1,1 -k2,2n unmerged_CDS.bed > sorted_CDS.bed

# Мерджим
bedtools merge -i sorted_CDS.bed > merged_CDS.bed

# Подсчет суммарного количества покрытых нуклеотидов
total_nucleotides=$(awk '{sum += $3 - $2} END {print sum}' merged_CDS.bed)

echo "total_nucleotides = $total_nucleotides"
