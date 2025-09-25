#!/bin/bash

OUTDIR="${BASE_DIR}/2_genome_exploration"
WORKDIR="${OUTDIR}/genome_for_homework"
mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

echo "[INFO] working in: $WORKDIR"

# download the reference and the annotation files

wget -nc -O genome_for_homework.fa.gz "https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz" 
gzip -d -f genome_for_homework.fa.gz

wget -O - https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz | \
   gzip -f -d > annotations_for_homework.gtf
    
# output of the number of contigs/chromosomes
echo "[INFO] Writing number_of_sequences.txt"
grep -c ">" genome_for_homework.fa > number_of_sequences.txt

# output of the header of each sequence
echo "[INFO] Writing headers_list.txt"
grep -r ">" genome_for_homework.fa > headers_list.txt

# output of the feature type of the 3rd column of the annotation file
echo "[INFO] Writing feature_type_counts.txt"
grep -v "#" annotations_for_homework.gtf | cut -f3 | sort | uniq -c > feature_type_counts.txt

# output of the count of exons of my interested gene
echo "[INFO] Writing Abcb11_ENSMUST00000102710_exon_count.txt"
grep "Abcb11" genome_for_homework/annotations_for_homework.gtf | grep 'ENSMUST00000102710'|cut -f 3| grep -c "exon" > Abcb11_ENSMUST00000102710_exon_count.txt

echo "[DONE] All files created in: ${WORKDIR}"