#!/bin/bash
#SBATCH -p ei-short                             # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 1                                    # number of cores
#SBATCH --mem 1G                                # memory pool for all cores
#SBATCH --output logs/%x.o%j                    # STDOUT and STDERR
#SBATCH --mail-type=end,fail                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

out_dir=alignments

source mafft-7.271

#Trim LSU
for marker in LSU
do

	#Trim
        singularity exec ~/programmes/trimAl/trimAl.img trimal \
                -in ${out_dir}/${marker}_aln_checked.fa \
                -fasta -gappyout > ${out_dir}/${marker}_aln_checked_trimmed.fa

done

for marker in $(cat markers)
do

        #Convert leading or trailing gaps to ?s
        awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ${out_dir}/${marker}_aln_checked_trimmed.fa | \
                awk -F '[^-](.*[^-]|$)' '{s=$0; h=gsub(/./,"?",$1); t=gsub(/./,"?",$2); print $1 substr(s,h+1, length(s)-h-t) $2}' > tmp.fa && mv tmp.fa ${out_dir}/${marker}_aln_checked_trimmed.fa

done
