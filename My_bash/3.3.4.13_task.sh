#!/usr/bin/env bash

grep -w "exon" gencode.v25.primary_assembly.annotation.gtf | grep -w "gene_type \"protein_coding\"" | grep -w "transcript_type \"protein_coding\"" > mRNA.gtf

cut -f9 mRNA.gtf | sed 's/;/\n/g' | grep "transcript_id" | sed 's/.*transcript_id "\([^"]*\)".*/\1/' | sort -u > transcript_ids.txt

unique_transcript_count=$(wc -l < transcript_ids.txt)
echo "Unique Transcript Count: $unique_transcript_count"

total_mRNA_length=$(awk '{sum += $5 - $4 + 1} END {print sum}' mRNA.gtf)
echo "Total mRNA Length: $total_mRNA_length"

average_mRNA_length=$(echo "scale=2; $total_mRNA_length / $unique_transcript_count" | bc)
echo "Average mRNA Length: $average_mRNA_length"

