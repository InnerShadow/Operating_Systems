#!/usr/bin/env bash

samtools view -bS -o Dm_rnaseq.bam Dm_rnaseq.sam
samtools sort -o Dm_rnaseq.sorted.bam Dm_rnaseq.bam

samtools index Dm_rnaseq.sorted.bam

samtools view -h Dm_rnaseq.sorted.bam | \
  awk '{ if($3=="2L" && $7=="2L" && $4>$8) print $4, $8 }' | \
  awk '{ intron_length = $4 - $8; if (intron_length > max_intron_length) max_intron_length = intron_length } END { print max_intron_length }'

