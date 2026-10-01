#mkdir genome_for_script/
#wget -nc -O genome_for_script/genome.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz


#gzip -d -f genome_for_script/genome.fna.gz
#wget -nc -O genome_for_script/annotations.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf.gz
#gzip -d -f genome_for_script/annotations.gtf.gz
#grep '>' genome_for_script/genome.fna >contig_headers.txt
#grep -c '>' genome_for_script/genome.fna > contig_count.txt