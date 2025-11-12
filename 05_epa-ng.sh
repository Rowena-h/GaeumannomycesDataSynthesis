#!/bin/bash
#SBATCH -p ei-short                             # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 3                                    # number of cores
#SBATCH --mem 10G                               # memory pool for all cores
#SBATCH --output logs/%x.o%j                    # STDOUT and STDERR
#SBATCH --mail-type=end,fail                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

lineage=gaeumannomyces
alignment_dir=alignments
raxml_dir=raxmlng
out_dir=epang/its

#Make output directory
mkdir -p ${out_dir}

#Convert leading or trailing ?s to gaps
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ${alignment_dir}/ITS_aln_checked_trimmed.fa | \
	awk -F '[^?](.*[^?]|$)' '{s=$0; h=gsub(/./,"-",$1); t=gsub(/./,"-",$2); print $1 substr(s,h+1, length(s)-h-t) $2}' > tmp && mv tmp ${alignment_dir}/ITS_aln_checked_trimmed.fa

source mafft-7.271

#Add GlobalFungi sequences to ITS alignment
mafft   --auto \
        --addfragments genus_Gaeumannomyces.fas \
        --reorder --keeplength \
        --thread ${SLURM_CPUS_PER_TASK} \
        ${alignment_dir}/ITS_aln_checked_trimmed.fa > ${alignment_dir}/genus_Gaeumannomyces_ITS_aln.fa

#Run RAxML-NG on on just ITS to get best model estimations
singularity exec ~/programmes/raxml-ng/raxml-ng-1.1.0.img raxml-ng \
	--evaluate \
	--msa ${alignment_dir}/ITS_aln_checked_trimmed.fa \
	--tree ${raxml_dir}/${lineage}_concat.raxml.support_rooted \
	--prefix ${raxml_dir}/constrained/${lineage}_its \
	--model GTR+G

#Run EPA
singularity exec ~/programmes/epa-ng/epa-ng.img epa-ng \
	-s ${alignment_dir}/ITS_aln_checked_trimmed.fa \
	-t ${raxml_dir}/${lineage}_concat.raxml.support_rooted \
	-q ${alignment_dir}/genus_Gaeumannomyces_ITS_aln.fa \
	--model ${raxml_dir}/constrained/${lineage}_its.raxml.bestModel \
	-w ${out_dir} \
	-T ${SLURM_CPUS_PER_TASK}
