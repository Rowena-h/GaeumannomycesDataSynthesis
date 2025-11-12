#!/bin/bash
#SBATCH -p ei-short                             # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 1                                    # number of cores
#SBATCH --mem 1G                                # memory pool for all cores
#SBATCH --output logs/%x.o%j                    # STDOUT and STDERR
#SBATCH --mail-type=end,fail                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

out_dir=alignments

#Concatenate
python ~/scripts/AMAS.py concat -f fasta \
				-d dna \
				-i $(find ${out_dir} -maxdepth 1 -name '*_aln_checked_trimmed.fa' ! -name '*genetree*') \
				-p ${out_dir}/gaeumannomyces_partition.txt \
				-t ${out_dir}/gaeumannomyces_concat.fa \
				-u fasta

#Add gene models
sed -i 's/^/GTR+G, /' ${out_dir}/gaeumannomyces_partition.txt
