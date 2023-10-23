#!/usr/bin/env bash

samtools view -b Dm_rnaseq.sam -o Dm_rnaseq.bam
samtools sort Dm_rnaseq.bam -o Dm_rnaseq.sorted.bam

samtools index Dm_rnaseq.sorted.bam

longest_intron_length=$(samtools view Dm_rnaseq.sorted.bam | awk '{ print $6 }' | cut -d 'N' -f 2 | sort -nr | head -n 1)

echo "Len: $longest_intron_length"
