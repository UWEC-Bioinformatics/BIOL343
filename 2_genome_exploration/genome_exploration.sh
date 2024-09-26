cd 2_genome_exploration
mkdir lemons
cd lemons
wget -nc -O lemon_genome.fa https://www.citrusgenomedb.org/citrus_downloads/Citrus_limon/C.limon_EMF-UC_v1-Alternative_genome/assembly/Climon_v1_alternative-scaffolds.fa
wget -O - https://www.citrusgenomedb.org/citrus_downloads/Citrus_limon/C.limon_EMF-UC_v1-Alternative_genome/annotation/Climon_v1_alternative-annotation.gff | \
    gzip -f -d > lemon_annotations.gff
gffread lemon_annotations.gff -T -o lemon_annotations.gtf
grep -c ">" lemon_genome.fa > num_chromosomes.txt
grep ">" lemon_genome.fa > chromosome_names.txt
cut -f 3 lemon_annotations.gtf | sort | uniq -c > feature_info.txt
grep 'CL5G054352012' lemon_annotations.gtf | grep -c "CDS" > GOI_exons.txt
samtools faidx lemon_genome.fa