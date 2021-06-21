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


## Tutorial




## Installation 

```r
remotes::install_github("IsmailM/matt")
```

## Citation

Manuscript in preparation




# TODO
Add specific version of data.table (e.g. if we use fcase etc.)

options("matt.placeholder_mode")
options("matt.cache_directory") -> default ./matt_cache_store


external: 
* read_cpg_data()

* annotate()
* annotate.nearest.neighbours()
* data_statistics()

* as.matt.table(, cols = ())
  * data.table
  * data.frame
  * BSSEq object
  * Genomic Ranges object
* as.data.frame()
* as.genomic.range()

* pair_wise_comparison
* all_vs_all
* collapse_all --> into a matrix

* plot

-> run:
- `document()`, `check()` `install()` `test()`


