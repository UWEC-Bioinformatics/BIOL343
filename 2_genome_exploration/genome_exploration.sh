#create a directory 
mkdir shell_files

#Download and save the reference genome
wget -nc -O shell_files/my_genome.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.20_GRCm38/GCF_000001635.20_GRCm38_genomic.fna.gz
gzip -d -f shell_files/my_genome.fna.gz

#Download and save the annotations
wget -nc -O shell_files/my_annotations.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.20_GRCm38/GCF_000001635.20_GRCm38_genomic.gtf.gz
gzip -f -d shell_files/my_annotations.gtf.gz

#Chromosome count as a text file
grep -c '>' shell_files/my_genome.fna > shell_files/chromosome_count.txt

#Chromosome list as a text file
grep '>' shell_files/my_genome.fna > shell_files/chromosome_list.txt
 

