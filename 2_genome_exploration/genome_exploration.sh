cd 2_genome_exploration
mkdir placenta_genome
cd placenta_genome
wget -nc -O placenta_genome.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz
gzip -d -f placenta_genome.fa.gz
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.gff.gz | gzip -f -d > placenta_annotations.gff
gffread placenta_annotations.gff -T -o placenta_annotations.gtf

samtools faidx placenta_genome.fa

grep -c ">" placenta_genome.fa > placenta_gene_count.txt
grep ">" placenta_genome.fa > placenta_genes.txt
cut -f 2 placenta_annotations.gtf | sort | uniq -c > placent_sort_cut_uniq.txt
grep "rna39844" placenta_annotations.gtf| grep -c "exon" > placenta_exon_count_gene_rna39844.txt