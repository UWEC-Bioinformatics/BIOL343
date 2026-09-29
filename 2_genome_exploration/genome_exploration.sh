mkdir genome1
wget -nc -O genome1/genome.fna.gz https://www.axolotl-omics.org/dl/AmexG_v6.0-DD.fa.gz
gzip -d -f genome1/genome.fna.gz

wget -nc -O genome1/annotations.gtf.gz https://www.axolotl-omics.org/dl/AmexT_v47-AmexG_v6.0-DD.gtf.gz
gzip -f -d genome1/annotations.gtf.gz

grep -c '>' genome1/genome.fna > contigs_chromosomes_count.txt

grep '>' genome1/genome.fna > contigs_chromosomes_names.txt