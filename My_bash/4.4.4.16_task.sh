#!/usr/bin/env bash

samtools view -b Dm_rnaseq.sam -o Dm_rnaseq.bam
samtools sort Dm_rnaseq.bam -o Dm_rnaseq.sorted.bam

samtools index Dm_rnaseq.sorted.bam

samtools view Dm_rnaseq.sorted.bam | awk '{ split($6, arr, "N"); for (i = 2; i <= NF; i++) { if (arr[i] ~ /^[0-9]+M$/) print arr[i] } }' > intron_lengths.txt

top_10_introns=$(sort -nr intron_lengths.txt | head -n 10)

echo "Top 10:"
echo "$top_10_introns"
