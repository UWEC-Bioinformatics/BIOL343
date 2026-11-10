#!/bin/bash
#SBATCH --job-name=hw8_marking_duplicates
#SBATCH --partition=week
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=marking_duplicates_slurm.out


picard MarkDuplicates -I /data/users/oswaldbn7080/BIOL343/6_alignment/my_alignment/star/Aligned.sorted.bam -M /data/users/oswaldbn7080/BIOL343/8_counting/homework/star_duplicates -O /data/users/oswaldbn7080/BIOL343/8_counting/homework/star.bam --VALIDATION_STRINGENCY SILENT
