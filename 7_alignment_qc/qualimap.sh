#!/bin/bash
#SBATCH --job-name=hw7_qualimap
#SBATCH --partition=week
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=qualimap_rnaseq.out

samtools sort -@ 32 -o /data/users/carltoga2473/BIOL343/6_alignment/alignment/star/axolotl/Aligned.sorted.bam /data/users/carltoga2473/BIOL343/6_alignment/alignment/star/axolotl/Aligned.out.bam

mkdir -p qualimap_axolotl/star

qualimap bamqc -nt 32 --java-mem-size=100G -outdir qualimap_axolotl/star/bam -bam /data/users/carltoga2473/BIOL343/6_alignment/alignment/star/axolotl/Aligned.sorted.bam  --feature-file /data/users/carltoga2473/BIOL343/2_genome_exploration/genome/annotations_axolotl.gtf
qualimap rnaseq --java-mem-size=100G -pe -outdir qualimap_axolotl/star/rnaseq -bam /data/users/carltoga2473/BIOL343/6_alignment/alignment/star/axolotl/Aligned.sorted.bam -gtf /data/users/carltoga2473/BIOL343/2_genome_exploration/genome/annotations_axolotl.gtf