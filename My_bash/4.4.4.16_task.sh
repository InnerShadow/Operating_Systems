#!/usr/bin/env bash

samtools view -bS Dm_rnaseq.sam > Dm_rnaseq.bam
samtools sort Dm_rnaseq.bam -o Dm_rnaseq.sorted.bam
samtools index Dm_rnaseq.sorted.bam

samtools view Dm_rnaseq.sorted.bam | cut -f 6 | grep -oE '(\d+)N' | tr -d 'N' | sort -nr | head -n 1
