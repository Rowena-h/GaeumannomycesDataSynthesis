#!/bin/bash
#SBATCH -p ei-short                      	# queue
#SBATCH -N 1                            	# number of nodes
#SBATCH -c 1                            	# number of cores
#SBATCH --mem 1000	                     	# memory pool for all cores
#SBATCH --output logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=end,fail            	# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk 	# send-to address

out_dir=alignments

source mafft-7.271

#Align each gene
for marker in $(cat markers)
do

	#Create alignment
	mafft ${out_dir}/${marker}.fasta > ${out_dir}/${marker}_aln.fasta

done
