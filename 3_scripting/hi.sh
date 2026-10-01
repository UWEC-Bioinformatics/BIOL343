mkdir genome_for_script/
wget -nc -O genome/genome.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/052/040/795/GCF_052040795.1_GRCz12ab/GCF_052040795.1_GRCz12ab_genomic.fna.gz
gzip -d -f genome/genome.fna.gz