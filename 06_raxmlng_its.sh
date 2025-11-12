#!/bin/bash
#SBATCH -p ei-short                             # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 1                                    # number of cores
#SBATCH --mem 10G                               # memory pool for all cores
#SBATCH --output logs/%x.o%j                    # STDOUT and STDERR
#SBATCH --mail-type=end,fail                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

alignment_dir=alignments
out_dir=raxmlng

#Make output directory
mkdir -p ${out_dir}

for lineage in ITS ITS1 ITS2
do

	#Convert file for RAxML-NG
	singularity exec ~/programmes/raxml-ng/raxml-ng-1.1.0.img raxml-ng \
		--parse \
		--msa ${alignment_dir}/${lineage}_aln_checked_trimmed.fa \
		--model GTR+G \
		--prefix ${out_dir}/gaeumannomyces_${lineage}

	#Run ML tree search and bootstrapping for <=1000 iterations
	singularity exec ~/programmes/raxml-ng/raxml-ng-1.1.0.img raxml-ng \
		--all \
		--msa ${alignment_dir}/${lineage}_aln_checked_trimmed.fa \
	        --model GTR+G \
	        --prefix ${out_dir}/gaeumannomyces_${lineage} \
	        --seed 2 \
	        --threads ${SLURM_CPUS_PER_TASK} \
	        --bs-trees autoMRE{1000}

	#Check convergence
	singularity exec ~/programmes/raxml-ng/raxml-ng-1.1.0.img raxml-ng \
		--bsconverge \
		--bs-trees ${out_dir}/gaeumannomyces_${lineage}.raxml.bootstraps \
	       	--prefix ${out_dir}/gaeumannomyces_${lineage}_convergence_test \
	        --seed 2

done
