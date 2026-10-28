#!/bin/bash
#SBATCH --job-name=star_second_pass
#SBATCH --partition=highmemory
#SBATCH --cpus-per-task=32
#SBATCH --mem=700G
#SBATCH --time=10:00:00
#SBATCH --output=second_pass.out

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
--outSAMattributes NH HI AS nM RG \
--outFileNamePrefix alignment/star/axolotl/ \
--sjdbFileChrStartEnd alignment/star/axolotl/SJ.out.tab

mv alignment/star/axolotl/Log.final.out alignment/star/axolotl/second-pass.final.out