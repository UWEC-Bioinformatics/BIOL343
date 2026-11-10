#!/bin/bash
#SBATCH --job-name=hw8_feature_counts_1
#SBATCH --partition=week
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=feature_counts_slurm_1.out

featureCounts /data/users/oswaldbn7080/BIOL343/8_counting/homework/star.bam -T 32 --byReadGroup --ignoreDup -M --fraction -s 1 -a /data/users/oswaldbn7080/BIOL343/2_genome_exploration/my_genome/annotations.gtf -o star_counts.tsv --verbose
