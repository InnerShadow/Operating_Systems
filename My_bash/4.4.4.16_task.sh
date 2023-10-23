#!/usr/bin/env bash

fastq-dump --split-3 SRR3146648

samtools view -b Dm_rnaseq.sam -o Dm_rnaseq.bam

samtools sort Dm_rnaseq.bam -o Dm_rnaseq.sorted.bam

samtools index Dm_rnaseq.sorted.bam

samtools view Dm_rnaseq.sorted.bam | awk '{if ($6 ~ /N/) print $6}' | cut -f 2 -d 'N' | sort -n | tail -n 1

