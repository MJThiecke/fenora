<pre>
███████ ███████ ███    ██  ██████  ██████   █████  
██      ██      ████   ██ ██    ██ ██   ██ ██   ██ 
█████   █████   ██ ██  ██ ██    ██ ██████  ███████ 
██      ██      ██  ██ ██ ██    ██ ██   ██ ██   ██ 
██      ███████ ██   ████  ██████  ██   ██ ██   ██ 
<pre>
**Feature normalisation for regression analysis**

Fenora provides functions for preprocessing and normalising sequencing-based 
feature count data, preparing it for downstream regression analyses. It 
supports a variety of assays such as ChIP-seq, CUT&Tag, ATAC-seq, RNA-seq, and 
more.

NOTE: Fenora does not make any assumptions on replicates. It is recommended to 
perform replicate-comparisons and merging before running Fenora.

---

## Overview

Fenora is designed to take **feature count matrices** (bed-like format with 
genomic coordinates) and perform:

- **Filtering** intervals by length  
- **Anscombe variance-stabilising transformation**  
- **Quantile normalisation** across chromosomes (optional)  
- **Diagnostic plots** to visualise normalisation effects  
- **Export** of processed counts (scores) tables for downstream regression analysis  

---

## Input format

Fenora expects a **bed-like file** with the following columns:

| chr | start | end | feature1 | feature2 | ... | featureN |
|-----|-------|-----|----------|----------|-----|----------|
| 1   | 1     | 16007 | 12345    | 3456     | ... | ...      |
| 1   | 16008 | 24571 | 5678     | 1234     | ... | ...      |

- **chr, start, end** → genomic interval (restriction fragment, peak, or window)  
- **feature1 ... featureN** → counts from quantitative sequencing assays 
(e.g., ChIP-seq, CUT&Tag, ATAC-seq, RNA-seq, etc.)  

---

## Installation

You can install Fenora from GitHub:

install.packages("devtools")
devtools::install_github("MJThiecke/fenora")

Dependencies:
- optparse
- data.table
- ggplot2
- ggrain
- aroma.light
- MASS
- fitdistrplus




