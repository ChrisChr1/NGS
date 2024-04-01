---
title: "BRCA"
output: word_document
date: "2024-03-31"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Enrichment analysis of Triple Negative Breast Cancer Patients vs Normal tissue. 

###Install packages
```{r, include=FALSE}
BiocManager::install("clusterProfile")
BiocManager::install("enrichplot")

library(clusterProfile)
library(enrichplot)
```

###Load in the data


```{r}
DEresults <- read.table(file = "tnbc_vs_normal_up_FDR05.txt", sep = "\t", header = T)
```
###Use EnrichGO to perform an enrichment analysis from the DE results. Create a table and write it to a file.

```{r}
enrich_results <- enrichGO(gene = DEresults$GeneSymbol, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL", readable = TRUE) 
head(enrich_results)
write.table(enrich_results, file="GOBP_enrichment_results_R.txt", sep = "\t", row.names=TRUE, col.names=NA)
```
### Visualization of enrichment with a bar plot
```{r}
barplot(enrich_results, x = "GeneRatio", color = "p.adjust", showCategory = 20)
```

### Visualization of enrichment with a dot plot
```{r}
dotplot(enrich_results, x = "Count", color = "qvalue", size = "GeneRatio",showCategory=30)
```


### Visualization of relationships between input genes with Gene-GO Term Network
```{r}
fold_change <- DEresults$logFC
names(fold_change) <- DEresults$GeneSymbol
cnetplot(enrich_results, color.params = list(foldChange = fold_change, edge = TRUE), cex.param = list(gene_label = 0.5), showCategory = 5)
```

###Visualization of relationships between input genes with a heatmap
```{r}
heatplot(enrich_results, foldChange=fold_change, showCategory = 10)
```
