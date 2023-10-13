#!/usr/bin/env bash

samtools view -b Ce_Pol2.sam -o Ce_Pol2.bam

samtools sort Ce_Pol2.bam -o Ce_Pol2.sorted.bam

samtools index Ce_Pol2.sorted.bam

samtools view Ce_Pol2.sorted.bam | cut -f 5 | sort | uniq

samtools view Ce_Pol2.sorted.bam | awk '{c += 1; if ($5==42) i +=1} END { print i/c}'

