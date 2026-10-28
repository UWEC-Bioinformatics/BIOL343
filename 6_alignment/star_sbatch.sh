#!/bin/bash
#SBATCH --job-name=star_sbatch
#SBATCH --partition=highmemory
#SBATCH --cpus-per-task=32
#SBATCH --mem=700G
#SBATCH --time=10:00:00
#SBATCH --output=star_sbatch.out

cd /data/users/carltoga2473/BIOL343/6_alignment
rm -rf ../2_genome_exploration/genome/star/_STARtmp

STAR --runThreadN 32 --runMode genomeGenerate \
  --genomeDir ../2_genome_exploration/genome/star \
  --genomeFastaFiles ../2_genome_exploration/genome/genome_axolotl.fna \
  --sjdbGTFfile ../2_genome_exploration/genome/annotations_axolotl.gtf \
  --sjdbOverhang 149 --genomeSAindexNbases 14 \
  --limitGenomeGenerateRAM 600000000000