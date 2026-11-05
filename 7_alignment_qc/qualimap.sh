#!/bin/bash
#SBATCH --job-name=hw7_qualimap
#SBATCH --partition=week
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=1-00:00:00
#SBATCH --output=qualimap_rnaseq5.out

align=/data/users/carltoga2473/BIOL343/6_alignment/alignment/star/axolotl/Aligned.sorted.bam
annotation=/data/users/carltoga2473/BIOL343/2_genome_exploration/genome/annotations_axolotl.gtf

# samtools sort -@ 32 -o $align /data/users/carltoga2473/BIOL343/6_alignment/alignment/star/axolotl/Aligned.out.bam

mkdir -p qualimap_axolotl/star

#qualimap bamqc -nt 32 --java-mem-size=100G -outdir qualimap_axolotl/star/bam -bam $align  --feature-file $annotation


qualimap rnaseq --java-mem-size=240G -pe -outdir qualimap_axolotl/star/rnaseq -bam $align -gtf $annotation 

