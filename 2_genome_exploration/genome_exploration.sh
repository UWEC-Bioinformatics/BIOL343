#!/bin/bash

# ---- SLURM SETTINGS ---- #

# -- Job Specific -- #
#SBATCH --job-name="genome"	# What is your job called?
#SBATCH --output=%j-genome-out.txt	# Output file - Use %j to inject job id, like output-%j.txt
#SBATCH --error=%j-genome-err.txt	# Error file - Use %j to inject job id, like error-%j.txt

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

## # Week 3 Homework

## ***Due (pushed to your GitHub branch) on 10/3 by 11:59 pm***

## ## Convert your Week 2 Homework notebook to a Shell script

## - Create a new file called `genome_exploration.sh` and save it in the `2_genome_exploration/` directory
## - All downloaded files should be written to a new directory in `2_genome_exploration`
## - Save the output of `grep` commands to `.txt` files using the `>` operator 
##     - These files should include the ***exact*** same information as was written out in `2_genome_exploration.ipynb`
## - Run the file in BOSE with the command `bash genome_exploration.sh`
##     - Check to ensure the proper files were created
##     - Commit/push to your branch
##     - Dr. Wheeler should be able to re-run `genome_exploration.sh` to regenerate all the files associated with Week 2 Homework

mkdir genome_hw_script
wget -nc -O genome_hw_script/tintling_genome.fna.gz https://mushroomdb.brc.hu/files/genome_and_genes/CopciAB_new_jgi_20220113.fasta.gz
gzip -d -f genome_hw_script/tintling_genome.fna.gz

wget -nc -O genome_hw_script/tintling_annotation.gtf.gz https://mushroomdb.brc.hu/files/genome_and_genes/CopciAB_new_jgi_20220113.gtf.gz
gzip -d -f genome_hw_script/tintling_annotation.gtf.gz

grep -c '>' genome_hw_script/tintling_genome.fna > genome_hw_script/num_contigs.txt
grep '>' genome_hw_script/tintling_genome.fna > genome_hw_script/names_contigs.txt
