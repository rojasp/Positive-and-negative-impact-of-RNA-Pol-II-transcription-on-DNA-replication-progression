#!/bin/bash  


#Run samtools
date
echo --------------------------START SAMtools-----------------


for file in *.bam
do name=$(basename $file .bam) 
#sort bam files 
samtools sort -o ${name}_sorted.bam ${name}.bam;

#remove replicates
samtools rmdup -s ${name}.bam ${name}_NoDup.bam;

#generate index 
samtools index ${name}_NoDup.bam;

#get some statistics
samtools flagstat ${name}_NoDup.bam;


done

