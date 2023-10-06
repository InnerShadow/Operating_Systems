#!/usr/bin/env bash

# Обработка CDS (кодирующих последовательностей)
grep 'gene_type "protein_coding"' gencode.v25.primary_assembly.annotation.gtf | \
  awk -F'\t' '$3 == "CDS" {OFS="\t"; split($9, a, /[";]+/); print $1, $4-1, $5, a[2], "0", $7}' > unmerged_CDS.bed

bedtools sort -i unmerged_CDS.bed > sorted_CDS.bed
bedtools merge -i sorted_CDS.bed > merged_CDS.bed

# Подсчет количества нуклеотидов для CDS
total_nucleotides=$(awk '{sum += $3 - $2} END {print sum}' merged_CDS.bed)
echo "Total nucleotides in CDS: $total_nucleotides"

# Обработка экзонов
grep 'gene_type "protein_coding"' gencode.v25.primary_assembly.annotation.gtf | \
  awk -F'\t' '$3 == "exon" {OFS="\t"; split($9, a, /[";]+/); print $1, $4-1, $5, a[2], "0", $7}' > unmerged_exon.bed

bedtools sort -i unmerged_exon.bed > sorted_exon.bed
bedtools merge -i sorted_exon.bed > merged_exon.bed

# Подсчет количества нуклеотидов для экзонов
total_exon_nucleotides=$(awk '{sum += $3 - $2} END {print sum}' merged_exon.bed)
echo "Total nucleotides in exons: $total_exon_nucleotides"

# Подсчет процента покрытия кодирующими последовательностями
coverage_percentage=$(echo "scale=3; ($total_nucleotides / $total_exon_nucleotides) * 100" | bc)
echo "Coverage by coding sequences: $coverage_percentage%"
