#!/bin/bash

java -jar picard.jar MarkDuplicates I=Ce_Pol2.sorted.bam O=Ce_Pol2.markdup.bam M=file.metrics

samtools rmdup -sS Ce_Pol2.markdup.bam rmdup.bam

samtools view -q 10 -b rmdup.bam -o Ce_Pol2.dedup_q10.bam

count_original=$(samtools view Ce_Pol2.sorted.bam | grep "MtDNA" | wc -l)
count_filtered=$(samtools view Ce_Pol2.dedup_q10.bam | grep "MtDNA" | wc -l)

echo "Количество прочтений на митохондриальном геноме в исходном файле: $count_original"
echo "Количество прочтений на митохондриальном геноме в файле после фильтрации: $count_filtered"
