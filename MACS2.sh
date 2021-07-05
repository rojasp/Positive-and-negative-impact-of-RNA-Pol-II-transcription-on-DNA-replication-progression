#!/bin/bash


#Run MACS2
date

# Peaks were called in BrdU sorted bam files derived from Bowtie2 (Galaxy). 
# Create sample names (prefix) in to a file: 
# ls -1 *bam | sort | sed -r 's/-[exp-ctr]{1}-brdu.sorted.bam//g' | sort | uniq > brdu_sample_names.txt
# exp is IP. 
# ctr is input control. 

# Create an array with the unique sample names 
sample_files_brdu=($(cut -f 1 brdu_sample_names.txt))

# print out all the names stored in the array
echo "${sample_files_brdu[@]}"

echo --------------------------START MACS2-----------------
# call peaks with macs2  

for file in "${sample_files_brdu[@]}"
do
	IP_bam="${file}"-exp-brdu.sorted.bam
	Input_bam="${file}"-ctr-brdu.sorted.bam

	macs2 callpeak -t "$IP_bam" -c "$Input_bam" -g hs -s 36 -m 8 30 -p 0.00001 --to-large 
	--keep-dup 2 --nomodel --extsize 200 -B -n "${file}"-exp-ctr-callpeaks 
	
done




