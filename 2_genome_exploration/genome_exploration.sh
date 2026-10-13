# this script doesn't run because the genome is too big, but files already exist in directory genome
# output of grep commands into txt file added to homework 2
mkdir genome1
wget -nc -O genome1/genome_axolotl.fna.gz https://www.axolotl-omics.org/dl/AmexG_v6.0-DD.fa.gz
gzip -d -f genome1/genome_axolotl.fna.gz

wget -nc -O genome1/annotations_axolotl.gtf.gz https://www.axolotl-omics.org/dl/AmexT_v47-AmexG_v6.0-DD.gtf.gz
gzip -f -d genome1/annotations_axolotl.gtf.gz

grep -c '>' genome1/genome_axolotl.fna > contigs_chromosomes_count.txt

grep '>' genome1/genome_axolotl.fna > contigs_chromosomes_names.txt