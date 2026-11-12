#!/bin/bash

# ---- SLURM SETTINGS ---- #

# -- Job Specific -- #
#SBATCH --job-name="two_pass_align"	# What is your job called?
#SBATCH --output=%j-two_pass_align-out.txt	# Output file - Use %j to inject job id, like output-%j.txt
#SBATCH --error=%j-two_pass_align-err.txt	# Error file - Use %j to inject job id, like error-%j.txt

#SBATCH --partition=week	# Which group of nodes do you want to use? Use "GPU" for graphics card support
#SBATCH --time=7-00:00:00	# What is the max time you expect the job to finish by? DD-HH:MM:SS
##SBATCH --account=

# -- Resource Requirements -- #
#SBATCH --mem=64G		# How much memory do you need?
#SBATCH --ntasks-per-node=64	# How many CPU cores do you want to use per node (max 64)?
#SBATCH --nodes=1		# How many nodes do you need to use at once?
##SBATCH --gpus=0		# Do you require a graphics card? How many (up to 3 per node)? Remove the first "#" to activate.

# -- Email Support -- #
#SBATCH --mail-type=ALL		# What notifications should be emailed about? (Options: NONE, ALL, BEGIN, END, FAIL, QUEUE) 

# ---- YOUR SCRIPT ---- #
STAR \
--runThreadN 64 \
--runMode alignReads \
--genomeDir ../2_genome_exploration/genome_hw_script/star_hmwk \
--readFilesManifest hmwk_manifest.tsv \
--readFilesCommand zcat \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--outFileNamePrefix hmwk_alignment/hmwk_star/
mv hmwk_alignment/hmwk_star/Log.final.out hmwk_alignment/hmwk_star/first-pass.final.out
STAR \
--runThreadN 64 \
--runMode alignReads \
--genomeDir ../2_genome_exploration/genome_hw_script/star_hmwk \
--readFilesManifest hmwk_manifest.tsv \
--readFilesCommand zcat \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--outSAMattributes NH HI AS nM RG \
--outFileNamePrefix hmwk_alignment/hmwk_star/ \
--sjdbFileChrStartEnd hmwk_alignment/hmwk_star/SJ.out.tab
mv hmwk_alignment/hmwk_star/Log.final.out hmwk_alignment/hmwk_star/second-pass.final.out