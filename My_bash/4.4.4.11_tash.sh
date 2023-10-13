#!/usr/bin/env bash

wget -O ce_ens_86.fa.gz https://ftp.ensembl.org/pub/release-86/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz

gunzip ce_ens_86.fa.gz

samtools faidx ce_ens_86.fa > ce_ens_86.fa.fai

prefetch -v SRR3535767

fastq-dump --split-3 SRR3535767.sra > Ce_Pol2.fastq

mv SRR3535767.fastq Ce_Pol2.fastq

bowtie2-build ce_ens_86.fa genome

bowtie2 -S Ce_Pol2.sam -x genome -U Ce_Pol2.fastq

