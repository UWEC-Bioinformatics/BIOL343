
cd beegenome
wget -nc -O genome.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz
gzip -d -f genome.fa.gz
grep ">" --count genome.fa > week_3_hw_genome_count.txt
grep ">" genome.fa > week_3_hw_genome_lines.txt
grep 'XR_001705490.2' genome_annotations.gtf | cut -f 3,4 | grep 'exon' --count > week_3_hw_exon_count.txt
grep '"XR_001705490.2' genome_annotations.gtf > week_3_hw_gene_of_intrest.txt
