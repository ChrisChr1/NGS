---
title: "Exploring Sample Relationships in R"
output:
  pdf_document: default
date: "2024-03-19"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

In this tutorial we are working with mouse mammary tissue at three developmental stages: virgin, pregnant and lactating.

```{r, install packages and read in data}
library(pheatmap)
library(RColorBrewer)
library(BiocManager)
BiocManager::install("org.Mm.eg.db")
BiocManager::install("DESeq2")
BiocManager::install("genefilter")
library(org.Mm.eg.db)
library(DESeq2)
library(genefilter)
```

## Show counts, exp_counts and metadata object

```{r pressure, echo=FALSE}
exp_counts <- read.table(file="~/Downloads/mammary_counts_v1.txt", sep = "\t", header = T)
sample_data <- read.table(file="~/Downloads/sample_data.txt", sep = "\t", header = T, row.names=1)
head(counts)
nrow(exp_counts)
sample_data
```

```{r}
#use org.Mm.eg.db to create a two-column table called “mapping” that takes the Ensembl identifiers in exp_counts and looks up their gene symbol counterparts.
mapping <-select(org.Mm.eg.db, as.character(exp_counts$ENSEMBL), keytype = "ENSEMBL", column="SYMBOL")
head(mapping)
nrow(mapping)
```
```{r}
#Print the number of duplicate rows
d <- duplicated(mapping$ENSEMBL)
sum(d)
mapping <- mapping[!d,]
nrow(mapping)
```

```{r}
#merge expression counts and mapping
exp_counts <- merge(exp_counts, mapping, by = "ENSEMBL")
head(exp_counts)
```
```{r}
#Remove rows where ENSEMBL did not map to a gene symbol
missing <- is.na(exp_counts$SYMBOL)
exp_counts <- exp_counts[!missing,]
nrow(exp_counts)
``` 

```{r}
#Remove duplicate gene symbols
o <- order(rowSums(exp_counts[,c(2:13)]), decreasing=TRUE)
exp_counts <- exp_counts[o,]
d2 <- duplicated(exp_counts$SYMBOL)
exp_counts <- exp_counts[!d2,]
```

```{r}
sum(d2)
nrow(exp_counts)
```

```{r}
row.names(exp_counts) <- exp_counts$SYMBOL
exp_counts$SYMBOL <- NULL
exp_counts$ENSEMBL <- NULL
head(exp_counts)
```

```{r}
#Create deseq dataset
data_deseq <- DESeqDataSetFromMatrix(countData = exp_counts, colData = sample_data, design = ~ 1)
head(counts(data_deseq))
```

```{r}
nrow(data_deseq)
data_deseq <- data_deseq[ rowSums(counts(data_deseq)) > 1, ]
rld <- rlog(data_deseq, blind=FALSE)
```

```{r}
sampleDists <- dist(t(assay(rld)))
sampleDists
```

```{r}
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( rld$cell_type, rld$dev_stage, rld$replicate, sep="-" )
DistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- paste( rld$cell_type, rld$dev_stage, rld$replicate, sep="-" )
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
plotPCA(rld, intgroup = c("cell_type", "dev_stage"))
```

```{r}
geneVars <- rowVars(assay(rld))
geneVarsOrdered <- order(geneVars, decreasing = TRUE)
topVarGenes <- head(geneVarsOrdered, 50)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("cell_type","dev_stage")])
clear_col_names <- paste( rld$cell_type, rld$dev_stage, rld$replicate, sep=".")
topGenesHeatmap <- pheatmap(mat, annotation_col=df, labels_col = clear_col_names)
```

1.	ID Mapping
a.	What was the overall goal of the ID mapping section?
      In this step, we are cleaning up the data. By converting each gene into from one type to another, we maintain readability and consistency of ID.

b.	Describe at least two steps of the ID mapping where information was lost.
      Two steps that removed information after mapping: 1) removing duplicates and removing duplicates after the first duplicated item 2) ENSEMBL failed to map and we removed rows with missing values.
      
2.	In this tutorial, we filtered the counts table to remove rows where the total counts were less than 1. However, there are many ways to filter the counts to retain “more interesting” genes. Name a different filtering criterion, explain why you chose it, and show the R command that you would apply to data_deseq to filter according to this criterion.
    data_deseq_1 <- data_deseq[ rowSums(counts(data_deseq)) < 1, ]
    data_deseq <- data_deseq_1[!data_deseq_1,]
    
3.	Sample distance heatmap
a.	Show an image of your sample distance heatmap. (look above)
b.	Discuss the relationships among the samples as seen in the heatmap. Which samples are more closely related? Which are less similar? Does this correspond to what you would expect based on the biology of the experiment? Explain.
      Luminal cell tissue is more closely related to luminal cell tissue than basal cell tissue samples. Basal cell tissue is more closely related to basal cell tissue than luminal cell tissue. Each tissue is closely related by developmental stage. You can identify these relationships through chared cohorts in the dendrogram.
      
4.	PCA plot
a.	Show an image of your PCA plot. 
b.	What percent of the variance is captured by:
    PC1 =  81%
    PC2 = 12%
c.	PC1 separates the samples into two major groups. What sample characteristic distinguishes these groups? 
      It looks as though PC1 separates the samples into basal and luminal groups. 
      
d.	What sample characteristic is mainly distinguished by PC2?
      PC2 distinguishes developmental stages.
      
e.	Do the replicates agree well with each other?
      I believe that the samples agree well with each other. The pregnant basal cells seems to have a bit more variance between the replicates compared to other sample replicates. 
      
5.	Gene heatmap
a.	Provide an image of your heatmap.
b.	Look at the clustering of the samples (columns of the heatmap). The tree diagram has four levels of branching. Which samples are separated from each other at the first branching level of the tree? The second branching level? The third branching level? The fourth level?
      The fourth level separates into different cell types: basal and luminal. The third level separates into virgin as one group and pregnant and lactating as another group. The second level separates lactating ang pregnant cells. Developmental stages are separated at the first level of the tree.
      
c.	Look at the clustering of the genes (rows of the heatmap). The first branch of the tree diagram separates the genes into two clusters--a larger cluster of ~35 genes and a smaller cluster of ~15 genes. For the genes in the larger cluster, what samples show high expression of these genes? Low expression?
      For the genes in the larger cluster, we mainly see luminal cells showing low expression (with the exception of virgin luminal cells with high expression) and basal cells showing higher expression.
      
For the genes in the smaller cluster, what samples show the highest expression? The lowest expression? Intermediate levels of expression?
    For the genes in the smaller cluster, we mainly see virgin luminal cells and all basal cells showing low expression. Lactating luminal cells and show high expression and pregnant luminal cells show moderate expression.
    
d.	Pick one gene from each of the two clusters in the heatmap and look up their function in bioinformatics resources such as UniProtKB or GeneCards. Based on the functional information you find, is the expression pattern of these genes in the different cell types and developmental stages what you expect? Explain. 
      I chose Dmbt1 (Deleted In Malignant Brain Tumors 1) because it is expressed more in virgin luminal cells and downregulated in every other cell sample. According to UniProtKB, this gene is involved in providing mucosal and cellular immunity and is expressed highly from 18.5 dpc (days post coitum) to birth and gradually decreases as the mouse enters into adulthood. Therefore, it makes sense that this gene is upregulated in virgin mice. Yet, we do not see the same expression pattern in virgin basal cells. I expect this is due to difference in the cellular functions. The next gene I chose was ACTG2 (Actin, gamma-enteric smooth muscle) which is mainly expressed in the basal cells. According to UniProt, this gene is expressed in smooth muscle.  
      
Note that mammary basal cells are also called myoepithelial cells. They resemble smooth muscle cells and secrete proteins that make up the basement membrane, a type of extracellular matrix. The luminal cells are responsible for making and secreting milk during lactation. A good paper discussing the role of mammary cell types can be found here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3193434/.
