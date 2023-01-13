#!/bin/bash

#Code adapted from code on Taylor lab github site

export NA=$1
mkdir ${NA}
cd ${NA}
ln -s ../${NA}.fastq

# initial alignment and processing of files
OUTPUT="$(bowtie2 -x /Users/johnwyrick/bin/bowtie2-2.2.6/saccer3 -U ${NA}.fastq -S ${NA}.sam)"
echo "${OUTPUT}"
samtools view -b -S ${NA}.sam >${NA}.bam
bamToBed -i ${NA}.bam >${NA}.bed

# get damage coordinate and sequence of damage
perl ../finddamagecoord.pl <${NA}.bed >${NA}_damage.bed

# get nucleotides of damage sites
fastaFromBed -s -name -fi ../saccer3_genome.fa -bed ${NA}_damage.bed -fo ${NA}_damage.fa

# count damage nucleotides
perl ../count_dinucleotides.pl <${NA}_damage.fa >${NA}_nucount.txt

# extract Dipyrimidine reads
printf "${NA}_damage.bed\n${NA}_damage.fa\n" | perl ../extract_dipyrimidine_reads.pl >${NA}_dipy.bed

# sort dipyrimidine reads 
cat ${NA}_dipy.bed | sort -k1,1 -k2,2n -k 6 >${NA}_dipy_sorted.bed

# split strands
printf "${NA}_dipy_sorted.bed\n" | perl ../split_strands.pl 

# rest of processing should be done mannually --> count function of IGV tools, withi window size 1, zoom 10, to generate .wig files
# then run perl set_background.pl >[bkgd set wig file]

