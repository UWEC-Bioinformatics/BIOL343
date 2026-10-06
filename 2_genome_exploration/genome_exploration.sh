# mkdir genome_for_script
# wget -nc -O genome_for_script/genome.fna.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/GRCm38.primary_assembly.genome.fa.gz
# gzip -d -f genome_for_script/genome.fna.gz
# head -20 genome_for_script/genome.fna
# tail genome_for_script/genome.fna
# wget -nc -O genome_for_script/annotations.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.annotation.gtf.gz
# gzip -f -d genome_for_script/annotations.gtf.gz
# grep -c '>' genome_for_script/genome.fna > contig_count.txt
# grep '>' genome_for_script/genome.fna > contig_header.txt
# head -10 genome_for_script/annotations.gtf 
# tail -10 genome_for_script/annotations.gtf 
