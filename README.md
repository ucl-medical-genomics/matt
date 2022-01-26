# MATT: Methylation Analysis ToolkiT
[![R-CMD-check](https://github.com/IsmailM/matt/workflows/R-CMD-check/badge.svg)](https://github.com/IsmailM/matt/actions)

> Pre-release (things may not fully work yet)

## Introduction


Methylation Analysis ToolkiT (MATT), was developed to deal with large numbers of Whole Genome Bisulfite Sequencing (WGBS) datasets in a computationally and memory-efficient manner.

A single typical WGBS dataset has approximately 28 million rows, and a typical analysis would involve multiple datasets. Hence it is important to ensure that all analysis code is efficient and scalable enough to be able to deal with large amounts of data. Common methylation analysis tasks include collapsing strand information, merging data from numerous samples, and pai-wise comparisons. These tasks involve creating union joins between large datasets based on the chromosome and position number. These processes are extremely computationally intensive and are not speed or memory-efficient even when using base R functions or the purpose-built functions from published libraries (e.g. bsseq, Methylkit and RnBeads).

Hence MATT was developed specifically to handle these large datasets.

MATT does not import any existing methylation libraries but is designed from the ground up for the vast amounts of data included in large cohort studies involving WGBS datasets.

`data.table` backend and an optional caching design to achieve high performance and memory efficiency.

MATT makes it feasible to work and analyse WGBS datasets on a large scale.


## Installation 

```r
remotes::install_github("IsmailM/matt")
```

> Note: XML install error may occur if R<=v4.0. Resolved by installing XML package version 3.9-0.


## Citation

> Manuscript in preparation



## Tutorial

```r
test_data <- system.file("extdata", "bismark_coverage_CpG.bedgraph", package="matt")
dt <- read_cpg_data(infile, align_to_reference=F)

# with local cache
dt <- read_cpg_data(infile, "sample_1", align_to_reference=F)

```


## Developer Quickstart

```bash
# Get repository
git clone https://github.com/ucl-medical-genomics/matt.git
cd matt

# Set Bioconductor sources
Rscript -e 'install.packages("BiocManager"); options(repos = BiocManager::repositories())'

# Install packages in DESCRIPTION file
Rscript -e 'install.packages()'
```

### Run tests

```bash
# Load `matt` 
devtools::load_all()
# Run all tests in tests/testthat
devtools::test()
```
