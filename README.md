<img align="right" width="140" height="150" src="https://github.com/gnrobinson/cloneR/blob/main/cloneR.png">

# cloneR - An R package to help clone correction in partially clonal populations
Authors: Guy N.M. Robinson and Douglas R. Cook

Artwork: Kaity Coleman

Department of Plant Pathology, University of California - Davis

---

## Introduction

cloneR is an R package that assigns individuals into clonal multilocus genotypes (MLGs) from genome-wide data by estimating each individual's ancestry.

Analysis of sexual recombination is an important first step to discerning the "evolvability" of organisms to novel environments. This can have practical implications for the design of vaccines, the breeding of crop resistance, the development of diagnostics, and legal ramifications of quarantine efforts. In partially clonal species, asexual cycles can obfuscate signal from sexual cycles and complicate genetic analyses. Clonal correction has to be undertaken before using many population analyses, since clonal cycles violate the Hardy-Weinberg assumption of random mating. Clonal correction is a method of compressing MLGs into a single representative isolate and estimate the effective population size. Recent efforts to gain a deeper understanding of the evolution of partially clonal organisms have utilized clonal correction as the basis of their evolutionary hypotheses.

Previous efforts to correct for clonal cycles have entailed either A) identifying multilocus genotypes from conserved sequence data (typically genes), or B) identifying recombination events during phylogenic construction. Identifying multilocus genotypes from sequence data from conserved genes has historically been the fastest method to determine clonal relationships. With the advent of genome-wide sequencing, this method has been significantly challenged because of the accumulation of erroneous mutations. Furthermore, conserved genes, such as the BUSCO genes, are under high purifying selection and therefore may under represent the total number of MLGs.

Identifying recombination events during phylogenetic construction is another common method to detect and subsequently correct for clonality in population datasets. While it is a conservative approach to ensure panmixia, it was not built to understand the complexity of evolutionary factors surrounding specie's diversity. There are many reasons why a lack of recombination may not be simply asexual cycles, most notably barriers to genetic flow such as geographic factors. Furthermore, it does not differentiate between sexual mechanisms and horizontal gene transfer.

## Theory behind CloneR

The frequency of clonal to sexual cycles is an important factor in determining the evolutionary trajectory of a species. Due to the hierarchical nature of phylogenetic trees, the clonal correction process can bias the frequency of sexual cycles in a population during linkage disequilibrium analysis. Choosing the hierarchical level to clonally correct to will change the linkage disequilibrium observed, since it is expected that lower linkage disequilibrium would be observed deeper into phylogenetic trees (Figure 1). Clonal correction methods have the possibility to both "under-correct" or "over-correct", each of which could change the evolutionary model used by researchers.

![](/phylogenetic_figure.png)

CloneR exploits advances in genome-wide sequencing to determine multilocus sequences including information about the standard deviation of sequence. It aims to determine whether two individuals are likely to have a shared ancestry. It builds upon another R package (LEA) that predicts population subdivision by predicting the population's ancestry. LEA utilizes a non-negative matrix factorization algorithm to calculate the ancestry matrix. By continually subdividing each individual's genome into K clusters, this separates isolates that are related clonally versus sexually because at continually higher K values clonal isolates do not continue to subdivide (Figure 2). 

![](/subdivision_figure.png)

Using the ancestry coefficients from LEA, CloneR will calculate the euclidean distances between isolates to determine their ancestry distance. Multiple randomly sub-sampled SNP sets are used to account for SNP choice bias and variance in ancestry coefficients based on SNP choice (Figure 3)

![](/variance_figure.png)

Once the grouping of clonal isolates is done, two output files can be generated. The first is a membership table that contains all of the group assignments for each individual. The second is a network plot showing the population structure between the different isolates with nodes representing the isolates and edges representing clonal relations.

![](/output_figure.png)

## Pipeline Overview

![](/cloneR_pipeline.png)

### 1. Installation
Install package via Github, using "devtools"
```{r}
# install dependencies
install.packages(c("LEA", "ggplot2", "ggraph", "igraph", "tools", "dplyr"))

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


