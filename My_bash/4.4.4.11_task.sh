#!/usr/bin/env bash

wget ftp://ftp.ensembl.org/pub/release-86/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz

gzip -d Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz

mv Caenorhabditis_elegans.WBcel235.dna.toplevel.fa ce_ens_86.fa

samtools faidx ce_ens_86.fa

prefetch SRR3535767
fastq-dump --split-3 SRR3535767

mv SRR3535767.fastq Ce_Pol2.fastq

bowtie2-build ce_ens_86.fa ce_ens_86_genome

bowtie2 -S Ce_Pol2.sam -x ce_ens_86_genome -U Ce_Pol2.fastq

