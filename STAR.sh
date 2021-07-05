#!/bin/bash  


#Run STAR
date
echo --------------------------START STAR-----------------

# download the human genome GRCh38
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.p10.genome.fa
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz

GTF=$(gencode.v28.annotation.gtf)  
echo "$GTF"  
GRC_fa=$(GRCh38.p10.genome.fa)  
echo "$GRC_fa" 

#Generate genome index for first pass 
STAR --runMode genomeGenerate --genomeDir . --genomeFastaFiles $GRC_fa --sjdbOverhang 100 --sjdbGTFfile $GTF;

#Generate alignments for first pass
for i in *.fastq.gz; 

do 

STAR --runThreadN 16 --runMode  alignReads --genomeDir . alignReads --readFilesIn $i --readFilesCommand zcat
--outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --alignIntronMax 500000 --alignMatesGapMax 1000000 
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outSAMtype BAM SortedByCoordinate 
--outWigType bedGraph --quantMode GeneCounts --outFileNamePrefix ./star_out/$i;
done


