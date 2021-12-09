# cloneR - An R package to help clone correction in partially clonal populations

## Authors: Guy N.M. Robinson and Douglas R. Cook

---
## Introduction
cloneR is an R package to identify isolates with the same multilocus genotypes. 
---
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

