---
title: "edgeR"
output: pdf_document
date: "2024-03-26"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## edgeR Tutorial
We are using raw read counts, not normalized counts as input to edgeR. 

```{r, echo=FALSE, include=FALSE}
library(edgeR)
library(org.Hs.eg.db)
library(dplyr)
library(statmod)
```

```{r}
cancer_counts <- read.table(file="oral_carcinoma_counts_2018.txt", sep =
"\t", header = T)
head(cancer_counts)
```

## ID Mapping

Map the ENSEMBL gene identifiers to Entrez Gene IDs

```{r}
mapping <-AnnotationDbi::select(org.Hs.eg.db, as.character(cancer_counts$ENSEMBL), keytype = "ENSEMBL", column="ENTREZID")

head(mapping)
```
## Remove duplicates

Use the duplicated function to deduplicate the rows in the mapping table

```{r}
d <- duplicated(mapping$ENSEMBL)
sum(d)
mapping <- mapping[!d,]
```
## Merge mapping and cancer counts


```{r}
cancer_counts <- merge(cancer_counts, mapping, by = "ENSEMBL")
head(cancer_counts)
```

## Remove missing data and duplicate ENSEMBL IDs

```{r}
missing <- is.na(cancer_counts$ENTREZID)
sum(missing)
cancer_counts <- cancer_counts[!missing,]
o <- order(rowSums(cancer_counts[,c(2:7)]), decreasing=TRUE)
cancer_counts <- cancer_counts[o,]
d2 <- duplicated(cancer_counts$ENTREZID)
sum(d2)
cancer_counts <- cancer_counts[!d2,]
```

##Create new mapping table

```{r}
mapping2 <- AnnotationDbi::select(org.Hs.eg.db, as.character(cancer_counts$ENTREZID),
keytype = "ENTREZID", column="SYMBOL")
d3 <- duplicated(mapping2$ENTREZID)
sum(d3)
#Remove duplicates
mapping2 <- mapping2[!d3,]

#Merge the mapping2 and cancer counts table
cancer_counts <- merge(cancer_counts, mapping2, by= "ENTREZID")
head(cancer_counts)
# Column 1 = ENTREZID
# Column 2= ENSEMBL
# Column 9 = Symbol
```

## MDS plot

```{r}
y <- DGEList(counts=cancer_counts[,3:8], genes=cancer_counts[,c(1:2,9)])
head(y$genes)
head(y$samples)
head(y$counts)
rownames(y$counts) <- rownames(y$genes) <- y$genes$ENTREZID

y$genes$ENTREZID <- NULL
y <- calcNormFactors(y)
y$samples$group = c("N", "T", "N", "T", "N", "T")
plotMDS(y)
```

## Dispersion and BCV

```{r}
y <- estimateDisp(y)
y$common.dispersion
#There is little variability between replicates.
plotBCV(y)
```

## Differentially expressed genes—exact test

```{r}
et <- exactTest(y, pair=c("N","T"))
summary(de<-decideTestsDGE(et))
```

```{r}
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h=c(-1, 1), col="blue")
```

```{r}
diffExpGenes <- topTags(et, n=1000, p.value = 0.05)
head(diffExpGenes)
write.table(diffExpGenes$table, file="tumor_v_normal_exactTest.txt", sep =
"\t", row.names=TRUE, col.names=NA)
```

## Differentially expressed genes—generalized linear model

```{r}
Patient <- factor(c(8,8,33,33,51,51))
Tissue <- factor(c("N","T","N","T","N","T"))
design <- model.matrix(~Patient+Tissue)
rownames(design) <- colnames(y)
design
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
#When design matrix was taken into account, the common dispersion decreased. This may be due to incorporating sample type information into the estimation.
```

```{r}
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=4)
summary(de2 <- decideTestsDGE(lrt))
diffExpGenes2 <- topTags(lrt, n=1000, p.value = 0.05)
head(diffExpGenes2$table)

```


1. ID Mapping. You have a data frame in R called “counts” that contains gene symbols
in the first column called “SYMBOL” and integer read counts from six human
samples in the subsequent columns. You want to map the gene symbols to Ensembl
IDs.
a. Show the R command you would use to create a table that maps between the
symbols and the Ensembl IDs.
      
      mapping <- select(org.Hs.eg.db, as.character(counts$SYMBOL), keytype = "SYMBOL", column="ENSEMBL")

b. Show the R commands you would use to determine whether there are any
Ensembl IDs in your list that map to multiple gene symbols.
      
      d <- duplicated(mapping$SYMBOL)
    
      sum(d)
c. Show the R commands you would use to determine whether there are gene
symbols that do not have a corresponding Ensembl ID.
      
      counts <- merge(counts, mapping, by = "SYMBOL")
      
      missing <- is.na(counts$ENSEMBL)
      
      sum(missing)
2. MDS plot.
a. Provide a screenshot of your MDS plot.
    
    Provided above.
b. What sample characteristic is distinguished by dimension 1 (horizontal
axis)?
   
    Tumor
c. What sample characteristic is distinguished by dimension 2 (vertical axis)?
    
    Patient
3. Dispersion and BCV.
a. What was the common dispersion for the exact test? For the generalized
linear model?
      
      The common dispersion was 0.21 for the exact test and 0.16 for the linear model.
b. Provide a screenshot of your BCV plot for the exact test.
     
      Provided above.
c. Is the common BCV in the range of what you would expect given the nature
of the experiment? Why or why not?
      
      Considering these are human samples, I would expect that the common BCV to be near or above 0.4, indicating that there is on average 40% variability in the expression of genes across diseased and healthy groups. We also see an increase of BVC as counts increase which I would not expect.
      
4. Differentially expressed genes—exact test
a. How many significantly up-regulated, significantly down-regulated, and non-
significant genes did you find?

      Down-regulated =     571
      Not Significant =  9650
      Up-regulated =      222

b. Provide a screenshot of your log-fold change vs. average log CPM plot.
      
      Provided above.
c. In what range(s) of fold-change do most of the significantly differentially
expressed genes lie?
      1:6 for up-regulated and -1:-7 for doqn-regulated
d. What were the top five most differentially regulated genes? Were they up- or
down-regulated?
      ATP2A1 - Down-regulated
      SH3BGRL2 - Down-regulated
      PYGM - Down-regulated
      PTHLH - Up-regulated
      SASH1 - Down-regulated
      
5. Differentially expressed genes—generalized linear model
a. Provide a screenshot of your design matrix.
      Provided above.
b. How many significantly up-regulated, significantly down-regulated, and non-
significant genes did you find?
             
    Down-regulated = 936
    Not Significant = 9175
    Up-regulated = 332
    
c. What were the top five most differentially regulated genes? Were they up- or
down-regulated?
    PTGFR - Down-regulated
    PTHLH - Up-regulated
    IGF1 - Down-regulated
    COL4A6 - Up-regulated
    ABCA8 - Down-regulated
    
d. Were there any genes in the top five that did not appear in the top five of the
exact test? If so, what p-value and rank did these genes have in the exact test?
    ABCA8 is ranked 8th with pval 2.04E-11, PRGFR is ranked 6th with pval 3.27E-13, IGF1 is ranked 262nd with pval            9.51E-05 and COL4A6 is ranked 16th with pval 1.59E-10.

e. Show how you would modify the R command lrt <- glmLRT(fit, coef=4), so
that you are testing the patient-specific effect of Patient 51 rather than the
tumor-specific effect. (You don’t have to carry out this analysis—just show
    >Patient <- factor(c(8,8,33,33,51,51))
   
    >design <- model.matrix(~Patient)
   
    >rownames(design) <- colnames(y)
   
    >design
   
    >y <- estimateDisp(y, design, robust=TRUE)
   
    >y$common.dispersion
   
    >fit <- glmFit(y, design)
   
    >lrt <- glmLRT(fit, coef=3)
   
    >summary(de2 <- decideTestsDGE(lrt))
   
    >diffExpGenes2 <- topTags(lrt, n=1000, p.value = 0.05)
    
