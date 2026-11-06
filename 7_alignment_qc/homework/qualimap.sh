#!/bin/bash
#SBATCH --job-name=hw7_qualimap
#SBATCH --partition=week
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=qualimap_rnaseq.out

qualimap rnaseq --java-mem-size=64G -outdir /qualimap/star/rnaseq -bam ../../6_alignment/my_alignment/star/Aligned.sorted.bam -gtf ../../2_genome_exploration/my_genome/annotations.gtf