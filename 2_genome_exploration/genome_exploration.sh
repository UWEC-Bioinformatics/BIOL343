#!bin/bash
#Create directory for for genome and annotation file
mkdir -p genome
#Retrieve and decompress genome 
wget -nc -O genome/genome.fa.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/706/575/GCF_015706575.1_ASM1570657v1/GCF_015706575.1_ASM1570657v1_genomic.fna.gz"
gzip -df genome/genome.fa.gz
#Retrieve and decompress the annotations
wget -nc -O genome/annotations.gff.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/706/575/GCF_015706575.1_ASM1570657v1/GCF_015706575.1_ASM1570657v1_genomic.gff.gz"
gzip -df genome/annotations.gff.gz
#Count the number of contigs/chromosomes
grep -c "^>" genome/genome.fa > contig_count.txt
#List the name of each contig/chromosom in the header of the sequence
grep "^>" genome/genome.fa > contig_headers.txt
