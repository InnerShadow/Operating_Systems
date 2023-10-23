#!/usr/bin/env bash

amtools view -F 4 Dm_rnaseq.sorted.bam | cut -f 1 | sort | uniq -c | awk -F " " '{print $1}' | sort -n| uniq -c > counts.txt 
head -1 counts.txt 
cat counts.txt | awk -F " " '{sum+=$1} END {print sum}'

