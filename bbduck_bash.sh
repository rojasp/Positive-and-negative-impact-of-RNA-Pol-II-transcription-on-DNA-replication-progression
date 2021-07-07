#!/bin/bash  


#Run bbduck
date
echo --------------------------START bbduck-----------------


for i in *.fastq.gz; 
do 

./bbduk.sh -Xmx4g in=$i  out=clean_$i ktrim=r k=21 ref=/path_to_bbmap/resources/adapters.fa

done


