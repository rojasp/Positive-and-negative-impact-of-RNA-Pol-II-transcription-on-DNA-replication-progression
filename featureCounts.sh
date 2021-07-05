#!/bin/bash  


#Run featureCounts
date

# download the human genome GRCh38 annotation
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz

echo --------------------------START featureCounts -----------------
GTF=$(gencode.v28.annotation.gtf)  
echo "$GTF"  
featureCounts -a $GTF -t gene -g gene_id -o counts.txt *.NoDup.bam

done


