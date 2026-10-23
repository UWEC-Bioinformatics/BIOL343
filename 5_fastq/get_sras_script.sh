#!/bin/bash

# ---- SLURM SETTINGS ---- #

# -- Job Specific -- #
#SBATCH --job-name="get_sras"	# What is your job called?
#SBATCH --output=%j-get_sras.txt	# Output file - Use %j to inject job id, like output-%j.txt
#SBATCH --error=%j-get_sras.txt	# Error file - Use %j to inject job id, like error-%j.txt

#SBATCH --partition=week	# Which group of nodes do you want to use? Use "GPU" for graphics card support
#SBATCH --time=7-00:00:00	# What is the max time you expect the job to finish by? DD-HH:MM:SS
##SBATCH --account=

# -- Resource Requirements -- #
#SBATCH --mem=64G		# How much memory do you need?
#SBATCH --ntasks-per-node=32	# How many CPU cores do you want to use per node (max 64)?
#SBATCH --nodes=1		# How many nodes do you need to use at once?
##SBATCH --gpus=0		# Do you require a graphics card? How many (up to 3 per node)? Remove the first "#" to activate.

# -- Email Support -- #
#SBATCH --mail-type=ALL		# What notifications should be emailed about? (Options: NONE, ALL, BEGIN, END, FAIL, QUEUE) 

# ---- YOUR SCRIPT ---- #

module load python-libs
source activate biol343-something
while IFS= read -r line; do \
    echo "Getting $line from NCBI SRA"; \
    parallel-fastq-dump --split-files --sra-id $line --threads 32 --outdir hmwk_fastq --gzip; \
    done < hmwk_sra_list.txt