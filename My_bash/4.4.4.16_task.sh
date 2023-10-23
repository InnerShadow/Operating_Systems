#!/usr/bin/env bash

samtools view -b Dm_rnaseq.sam -o Dm_rnaseq.bam
samtools sort Dm_rnaseq.bam -o Dm_rnaseq.sorted.bam

samtools index Dm_rnaseq.sorted.bam
samtools view Dm_rnaseq.sorted.bam | awk $6 '{ print $6 }' | cat | grep -o ".M[^N]*" | cut -c 3- | sort -nr | less > 1212.txt