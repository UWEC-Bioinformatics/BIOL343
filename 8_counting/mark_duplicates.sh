#!/bin/bash
#SBATCH --job-name=hw8_markduplicates
#SBATCH --partition=week
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=1-00:00:00
#SBATCH --output=markduplicates1.out


mkdir logs
mkdir dedup

picard MarkDuplicates \
    -I /data/users/carltoga2473/BIOL343/6_alignment/alignment/star/axolotl/Aligned.sorted.bam \
    -M /data/users/carltoga2473/BIOL343/8_counting/star_duplicates \
    -O /data/users/carltoga2473/BIOL343/8_counting/star.bam \
    --VALIDATION_STRINGENCY SILENT