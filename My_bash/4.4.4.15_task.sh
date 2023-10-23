#!/usr/bin/env bash

samtools view -bS -o Dm_rnaseq.bam Dm_rnaseq.sam
samtools sort -o Dm_rnaseq.sorted.bam Dm_rnaseq.bam
samtools index Dm_rnaseq.sorted.bam

samtools view Dm_rnaseq.sorted.bam | cut -f 6 | grep "N" | tr -d "N" | tr "," "\n" | sort -n | tail -1
