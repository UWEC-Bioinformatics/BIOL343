mkdir hw3_genome
cd hw3_genome
wget -nc -O zebrafish_genome.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.5_GRCz10/GCF_000002035.5_GRCz10_genomic.fna.gz
gzip -d -f zebrafish_genome.fa.gz
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.5_GRCz10/GCF_000002035.5_GRCz10_genomic.gtf.gz | gzip -f -d > zebrafish_annotations.gtf
grep -c ">" zebrafish_genome.fa > num_chromosomes.txt
grep -F ">" zebrafish_genome.fa > chromosome_headers.txt
cut -f 3 zebrafish_annotations.gtf | sort | uniq -c > annotation_features.txt
grep "cep97" zebrafish_annotations.gtf | grep "XR_001800919.1" | cut -f 3 | grep -c "exon" > num_exons_gene.txt