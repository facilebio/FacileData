
<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
========

The `FacileData` package was written to facilitate easier analysis of large, multi-assay high-throughput genomics datasets. To this end, the `FacileData` package provides two things:

1.  A *FacileData Access API* that defines a fluent interface over multi-assay genomics datasets that fits into the [tidyverse](https://www.tidyverse.org/). This enables analysts to more naturally query and retrieve data for general exploratory data analysis; and
2.  A reference implementation of a datastore that implements the *FacileData Access API* called a *FacileDataSet*. The `FacileDataSet` provides efficient storage and retrieval of arbitrarily large high-throughput genomics datasets. For example, a single `FacileDataSet` can be used to store *all* of the RNA-seq, microarray, RPPA, etc. data from the [The Cancer Genome Atlas](https://cancergenome.nih.gov/). This singular `FacileDataSet` allows analysts easy access to arbitrary subsets of these data without having to load all of it into memory.

Installation
============

The FacileData suite of packages is only available from github from now. You will want to install three `FacileData*` packages to appreciate the its utility:

``` r
# install.packages("devtools")
devtools::install_github("faciledata/FacileData")
devtools::install_github("faciledata/FacileExplorer")
devtools::install_github("faciledata/FacileTCGADataSet")
```

Example Usage
=============

As a teaser, we provide code snippets that show how to plot HER2 copy number vs expression across the TCGA "BLCA" and "BRCA" indications using the= `FacileDataSet`. We'll then compare that to how the same code might be written using more traditional bioconductor containers.

``` r
library(FacileData)
library(FacileTCGADataSet)
tcga <- FacileTCGADataSet()

features <- filter_features(tcga, name == "ERBB2")

fdat <- tcga %>% 
  filter_samples(indication %in% c("BLCA", "BRCA")) %>% 
  with_assay_data(features, assay_name = "rnaseq", normalized = TRUE) %>%
  with_assay_data(features, assay_name = "cnv_score") %>% 
  with_sample_covariates(c("indication", "sex"))

ggplot(fdat, aes(cnv_score_ERBB2, ERBB2, color=sex)) +
  geom_point() +
  facet_wrap(~ indication)
```

<img src="vignettes/images/her2_cnv_vs_expression.png" />
