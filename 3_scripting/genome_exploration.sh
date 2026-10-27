mkdir cowgenes
wget -nc -O cowgenes/genome.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.3_ARS-UCD2.0/GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz
gzip -d -f cowgenes/genome.fna.gz
wget -nc -O cowgenes/annotations.gtf https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9913/GCF_002263795.3-RS_2023_09/GCF_002263795.3_ARS-UCD2.0_genomic.gtf.gz
gzip -d -f cowgenes/annotations.gtf
grep -c '>' cowgenes/genome.fna > chromosome_count.txt
grep '>' cowgenes/genome.fna > chromosome_name.txt