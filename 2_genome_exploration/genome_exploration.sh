mkdir c_elegans_genome
cd c_elegans_genome

wget -nc -O c_elegans_genome.fa.gz https://downloads.wormbase.org/species/c_elegans/sequence/genomic/c_elegans.PRJNA275000.current.genomic.fa.gz
gzip -d -f c_elegans_genome.fa.gz

wget -O - https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS19.canonical_geneset.gtf.gz | \
    gzip -f -d > c_elegans_annotations.gtf

grep -c ">" c_elegans_genome.fa > num_chromosomes.txt

grep ">" c_elegans_genome.fa > name_chromosomes.txt

cut -f 3 c_elegans_annotations.gtf | sort | uniq -c > feature_types.txt

grep 'Y38C1AB.8a.1' c_elegans_annotations.gtf | cut -f 3 | grep 'exon' | uniq -c > num_exons_gene.txt

samtools faidx c_elegans_genome.fa
