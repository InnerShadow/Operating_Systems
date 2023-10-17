#!/usr/bin/env bash

wget http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz

gunzip dm6.fa.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR314/008/SRR3146648/SRR3146648.fastq.gz

gunzip SRR3146648.fastq.gz

hisat2-build -p 2 dm6.fa genome

hisat2 -S Dm_rnaseq.sam -x genome -U SRR3146648.fastq

