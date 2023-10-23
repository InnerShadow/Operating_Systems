#!/usr/bin/env bash

fastq-dump --split-3 SRR3146648
hisat2 -S Dm_rnaseq.sam -x genome -U Dm_rnaseq.fastq

samtools view -b Dm_rnaseq.sam -o Dm_rnaseq.bam

samtools sort Dm_rnaseq.bam -o Dm_rnaseq.sorted.bam

samtools index Dm_rnaseq.sorted.bam

max_intron_length=$(samtools view Dm_rnaseq.sorted.bam | awk '{ print $6 }' | grep -o ".*[^*]*" | cut -c 3- | tr ',' '\n' | sort -nr | head -n 1)

echo $max_intron_length

