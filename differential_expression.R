

##############################
# 0 - Load libraries
##############################
library(DESeq2)
library(easyRNASeq)
library(RColorBrewer)
library(vsn)
library(scatterplot3d)
library(VennDiagram)
library(LSD)
library(gplots)
library(arrayQualityMetrics)
library(parallel)

############################## 
# 1 - Source file 
##############################
## Use the counts.txt data generate by featureCounts.sh

count <- read.csv("counts.txt", sep = "\t", row.names = 1)
count <-as.matrix(count)
# class(count)

# meta: is a meta file with sample names and treatment (DRB or a-Amanitin)  

meta<- read.csv("meta", sep = "\t")
meta$treatment <- as.factor(meta$treatment)

# get some frequencies 
meta %>% group_by(treatment) %>% summarize(n=n()) %>% mutate(freq = n / sum(n))




###########################################
# 2- Differential expression analysis 
##########################################
#run DESeq2
dds <- DESeqDataSetFromMatrix(countData=count , 
                              colData=meta, 
                              design=~treatment,
                              tidy = FALSE)

dds <- DESeq(dds)
res <- results(dds)

#remove remove rows of the DESeqDataSet that have no counts
dds <- dds[rowSums(counts(dds)) > 1, ]
nrow(dds)
class(dds)

#estimate size factor, dispersion and plot
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- names(res)
sizes
boxplot(sizes,main="relative library sizes",ylab="scaling factor")
cds = estimateDispersions(dds)
plotDispEsts(cds)

#Transfor data on the log 2 scale 
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vsd.fast <- vst(dds, blind = TRUE)
vst <- assay(vsd.fast)
colnames(vst) <- colnames(countData)
#rownames(vst) <- rownames(countData)
write.table(vst, "gene_matrix_RNA_matrix" , sep="\t", col.names=T, row.names = T)


#Remove all values so non expressed genes have a vst value of 0
vst <- vst - min(vst)
class(vst)
class(vsd)
head(colData(vsd),3)
write.table(vst, "gene_matrix_RNA_matrix_over0", sep = "\t")

#PCA analysis
DESeq2::plotPCA(vsd, intgroup = "treatment")
norm_count <- counts(dds,normalized=T)
counts(dds)

#Get differential expressed genes 
dds2 <- DESeq(dds)
res2 <- results(dds2)
# write.table(resSig, file="DE_results")

#Order and filter data according with padj values 
resOrdered <- res2[order(res2$padj),]
summary(res2)
resSig <- subset(resOrdered, padj < 0.1)
resSig
res_lfc <- subset(resSig, abs(log2FoldChange) > 2) 
# write.table(resSig, file="DE_by_padj<0.1")

#Plot the 50 most expressed genes 

mat <- assay(vsd)[head(order(resSig$padj)),50] 
pheatmap(mat) 
res.05 <- results(dds2, alpha = 0.05)
resOrdered_pvalue <- res.05[order(res.05$pvalue),]
mat <- assay(vsd)[head(order(res.05$padj)), ]
pheatmap(mat) 


