# Adapted version to FASTQ files with removed primers, Anders Lanzen, Dec 2016

#Length of amplicon minus primers (787-519-15=253)
#(the rev. primer is 20 nt but starts at pos. 787)
crop_length=252 

#Make directory where the merged read files will be placed (fastq and fasta)
mkdir merged_reads

#Iterate over all subdirectories, assume they have the sample name, enter and
#uncompress fastq files
for d in $*; do
    echo $d
    cd $d
    pwd
    i=${d//\/}
    echo $i

    gunzip -kf *.fastq.gz 
    
    # Merge the reads (non-staggering since adaptors removed) max 5 differences

    echo $d >> ../readprep.log

    vsearch --fastq_mergepairs *_R1*.fastq --reverse *_R2*.fastq  --fastq_allowmergestagger --fastqout ${i}_m.fastq 2>> ../readprep.log

    vsearch --fastq_filter ${i}_m --fastaout ${i}_merged_QF.fasta --fastq_maxee 1 --fastq_trunclen $crop_length 2>>../readprep.log    
    
    #Fix read names to include "barcode" label
   python ~/script/drive5/fixreads.py ${i}_merged_QF.fasta $i > ../merged_reads/${i}_reads_fixed.fasta

   rm ${i}_merged_QF.fasta

   cd ..
done

echo "Preparing FASTQC reports"

cat */*_m.fastq > all_merged.fastq
fastqc all_merged.fastq
rm all_merged.fastq

cat */*R1*.fastq > all_R1.fastq
fastqc all_R1.fastq
rm all_R1.fastq

cat */*R2*.fastq > all_R2.fastq
fastqc all_R2.fastq
rm all_R2.fastq

rm */*.fastq
rm *fastqc.zip


