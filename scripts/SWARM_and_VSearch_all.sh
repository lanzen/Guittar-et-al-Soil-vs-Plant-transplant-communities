#!/bin/bash

scriptdir=$HOME/script/drive5

#Assumes that all reads are using 1 line for the seq and are in the current folder ending with "_read_fixed.fasta"

#Concatenate all files
cat *_fixed.fa* > All_reads.fasta
sed 's/ 1:N:0:.*;/;/g' All_reads.fasta > All_reads_fixed.fasta
rm All_reads.fasta

#Sort and dereplicate reads
vsearch -derep_fulllength All_reads_fixed.fasta -output All_reads_U.fa -sizeout

#Abundance sort and retain singletons 
vsearch -sortbysize All_reads_U.fa -output All_reads_U_sorted.fa -minsize 0
rm All_reads_U.fa

#Linearize for SWARM
awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}' All_reads_U_sorted.fa  > All_reads_U_fixed.fa
rm All_reads_U_sorted.fa

# Run SWARM
swarm -f -z -t 6 -w SWARM_OTUs.fasta All_reads_U_fixed.fa >/dev/null

#Chimera filtering
vsearch --uchime_ref SWARM_OTUs.fasta -db ~/db/rdp_gold.fa -strand plus -nonchimeras SWARM_OTUs_uchime_ref.fa

vsearch --uchime_denovo SWARM_OTUs_uchime_ref.fa --nonchimeras SWARM_OTUs_uchime_ref_denovo.fa

#OTU Clustering
vsearch --cluster_size SWARM_OTUs_uchime_ref_denovo.fa --id 0.97 --consout All_OTUs.fa

#Relabel OTUs
python $scriptdir/fasta_number.py All_OTUs.fa OTU_ > All_OTUs_Final.fa
rm All_OTUs.fa SWARM_OTUs_uchime_ref.fa

#Linearize
awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}' All_OTUs_Final.fa > All_OTUs_Final_L.fa
rm All_OTUs_Final.fa

#Map reads
vsearch -usearch_global All_reads_fixed.fasta -db All_OTUs_Final_L.fa -strand plus -id 0.97 -uc All_map.uc

python $scriptdir/uc2otutab.py All_map.uc > All_OTU_table.csv

# Taxonomic classification
blastn -task megablast -query All_OTUs_Final_L.fa -db ~/projects/PhyloRefDB/silvamod106/silvamod.fasta -num_alignments 100 -outfmt 5 -out OTUs_Silvamod.xml -num_threads 4
gzip OTUs_Silvamod.xml

classify -i All_OTUs_Final_L.fa -t  All_OTU_table.csv -k OTUs_Silvamod.xml.gz





