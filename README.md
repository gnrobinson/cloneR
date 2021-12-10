<img align="right" width="140" height="150" src="https://github.com/image.png">
# cloneR - An R package to help clone correction in partially clonal populations
Authors: Guy N.M. Robinson and Douglas R. Cook
Artwork: Kaitlyn Coleman
Department of Plant Pathology, University of California - Davis

---

## Introduction

cloneR is an R package to identify isolates with the same multilocus genotypes from genome-wide data

Analysis of sexual recombination is an important first step to discerning the "evolvability" of organisms to novel environments. This can have practical implications in agriculture where resistant crop cultivars are deployed against pathogen populations that are highly structured, leading to inconsistent practical success. In partially clonal species, asexual cycles can obfuscate signal from sexual cycles because it significantly violates Hardy-Weinberg assumptions. As a result, a process of "clonal correction" is undertaken prior to genetic analyses to the estimated effective population size.

Previous efforts to correct for clonal cycles have entailed either A) identifying multilocus genotypes, or B) identifying recombination events during phylogenic construction. Identifying multilocus genotypes from sequence data from conserved genes has historically been the fastest method to determine clonal relationships. With the advent of genome-wide sequencing, this method has been significantly challenged since much sexual recombination can occur outside of these conserved genes. Furthermore, this technique could not be applied to a genome-wide scale because of the accumulation of mutations during the evolutionary process. 

Identifying recombination events during phylogenetic construction is another common method to detect and subsequently correct for clonality in population datasets. While far more advanced that simply looking at sequence diversity, it was not built with clonal correction in mind. There are many reasons why a lack of recombination may not be simply asexual cycles, most notably barriers to genetic flow such as geographic factors. Furthermore, it does not differentiate between sexual mechanisms and horizontal gene transfer.

cloneR exploits advances in genome-wide sequencing to determine multilocus sequences including information about the standard deviation of sequence. It aims to determine whether two individuals are likely to have a shared ancestry that cannot be further subdivided at different K populations. It builds upon an R package (LEA) that predicts population subdivision by predicting the population's ancestry. cloneR utilizes a non-negative matrix factorization algorithm (PCA method to be implemented soon[12/2021]) to calculate the ancestry matrix. cloneR will then calculate the euclidean distances between isolates to deteremine their ancestry distance. Multiple randomly sub-sampled SNP sets are used to account for SNP choice bias.

## Pipeline Overview

![](/github/pipeline_figure.png)
### 1. Installation
Install package via Github, using "devtools"
```{r}
# install dependencies
install.packages(c("LEA", "ggplot2", "ggraph", "igraph", "tools"))

# install cloneR
devtools::install_github("gnrobinson/cloneR")
```

### 2. Running cloneR
cloneR uses an unzipped VCF file as input. It can be run in full or sequentially:
```{r}
library(cloneR)

cloneR(example.vcf, snps = 1000, subsets = 100, K = 2:20)
```

```{r}
library(cloneR)

# creates labels for ancestry
extract_labels(example.vcf)

# convert VCF file to geno file
LEA::vcf2geno(example.vcf, force = TRUE)

# make subsets of vcf file
make_subsets(example.geno, snps = 1000, subsets = 100)

# calculate Q matrices for all subsets
make_Q(K=2:20, CPU = 1, alpha = 10, iterations = 200, ploidy = 2)

# find euclidean distances for all subsets Q matrices
calc_dist(K, plot = TRUE)
```


