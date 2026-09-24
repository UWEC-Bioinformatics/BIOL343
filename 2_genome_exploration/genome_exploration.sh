#Make new directory for downloaded files
mkdir shell_genome

#Download and unzip (decompress) reference genome
wget -nc -O shell_genome/genome.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/237/465/GCF_015237465.1_rCheMyd1.pri/GCF_015237465.1_rCheMyd1.pri_genomic.fna.gz
gzip -d -f shell_genome/genome.fna.gz

#Download and unzip (decompress) reference genome annotations
wget -nc -O shell_genome/annotations.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/237/465/GCF_015237465.1_rCheMyd1.pri/GCF_015237465.1_rCheMyd1.pri_genomic.gtf.gz
gzip -f -d shell_genome/annotations.gtf.gz

#returns text file containing the number of chromosomes in the reference genome
grep -c '>' shell_genome/genome.fna > shell_genome/chromosome_count.txt

#returns a text file containing the name of every chromosome, contig, and scaffolding in the genome
grep '>' shell_genome/genome.fna > shell_genome/chromosome_names.txt