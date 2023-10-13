#!/usr/bin/env bash

wget https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec

gunzip ERR426367_1.fastq.gz

mv ERR426367_1.fastq Legionella_R1.fastq

sed -n '1~4s/^@/>/p;2~4p' Legionella_R1.fastq > Legionella_R1.fa

makeblastdb -in UniVec -parse_seqids -dbtype nucl

blastn -db UniVec -reward 1 -penalty -5 -gapopen 3 \
-gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp \
1750000000000 -query Legionella_R1.fa -outfmt 6 -out file.format6.out

cut -f2 file.format6.out | sort > file.txt

cat file.txt | grep "gnl|uv|NGB01096.*:*-**" | wc -l

grep "gnl|uv|NGB00726.*:*-**" UniVec

echo "Top 3 most frequent sequences:"
sort file.txt | uniq -c | sort -nr | head -n 3
