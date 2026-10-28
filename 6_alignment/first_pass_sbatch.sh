#!/bin/bash
#SBATCH --job-name=star_first_pass
#SBATCH --partition=highmemory
#SBATCH --cpus-per-task=32
#SBATCH --mem=700G
#SBATCH --time=10:00:00
#SBATCH --output=first_pass.out

cd /data/users/carltoga2473/BIOL343/6_alignment
mkdir -p alignment/star/axolotl/

STAR \
--runThreadN 32 \
--runMode alignReads \
--genomeDir ../2_genome_exploration/genome/star \
--readFilesManifest manifest_axolotl.tsv \
--readFilesCommand zcat \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--outFileNamePrefix alignment/star/axolotl/

mv alignment/star/axolotl/Log.final.out alignment/star/axolotl/first-pass.final.out