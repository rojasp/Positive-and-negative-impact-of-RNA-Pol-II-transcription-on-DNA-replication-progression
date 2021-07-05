#!/bin/bash  

#run fastqc
date
echo --------------------------START fastq -----------------

#create a folder to store the results
mkdir -p qa/raw
  
#creat a directory to hold the raw data
mkdir raw

fastqc -o qa/raw -t 16 raw/*.fastq.gz

# run multiqc to summarize all data in one plot
cd qa/raw
multiqc .


