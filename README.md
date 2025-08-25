<pre>
             ███████ ███████ ███    ██  ██████  ██████   █████  
             ██      ██      ████   ██ ██    ██ ██   ██ ██   ██ 
             █████   █████   ██ ██  ██ ██    ██ ██████  ███████ 
             ██      ██      ██  ██ ██ ██    ██ ██   ██ ██   ██ 
             ██      ███████ ██   ████  ██████  ██   ██ ██   ██ 
<pre>
             ~~Feature normalisation for regression analysis~~

Fenora provides functions for preprocessing and normalising sequencing-based 
feature count data, preparing it for downstream regression analyses. It 
can be used for data from assays such as ChIP-seq, CUT&Tag, ATAC-seq, RNA-seq, 
etc.
Fenora is particularly well suited to pre-process count data in intervals of
varying sizes (such as restriction fragments).
Fenora is also capable of mitigating read-count effects of aneuploidy, using
the argument --qNormBetweenChr.

NOTE: Fenora does not make any assumptions on replicates. It is recommended to 
perform replicate-comparisons and merging before running Fenora.


###############################################################################
# Overview

Fenora is designed to take a feature count matrix (bed-like format with genomic
coordinates) and perform:

- Filtering intervals by length  
- Anscombe variance-stabilising transformation
- Quantile normalisation between features
- Quantile normalisation between chromosomes (optional) 
- Correcting for interval-size (makes intervals of different sizes comparable) 
- Diagnostic plots to visualise normalisation effects  
- Export of processed counts (scores) for downstream regression analysis  


###############################################################################
# Input format

Fenora expects a bed-like file with the following columns:

| chr | start | end | feature1 | feature2 | ... | featureN |
|-----|-------|-----|----------|----------|-----|----------|
| 1   | 1     | 16007 | 12345    | 3456     | ... | ...      |
| 1   | 16008 | 24571 | 5678     | 1234     | ... | ...      |

- chr, start, end: genomic interval (restriction fragment, peak, or window)  
- feature1 ... featureN: counts from quantitative sequencing assays (e.g., 
  ChIP-seq, CUT&Tag, ATAC-seq, RNA-seq, etc.)  


###############################################################################
# Installation

You can install Fenora from GitHub:

install.packages("devtools")
devtools::install_github("MJThiecke/fenora")

A good starting point is to run scripts/fenora.R (See 'Test run' below)

Dependencies:
- optparse
- data.table
- ggplot2
- ggrain
- aroma.light (install with BiocManager::install('aroma.light'))
- MASS
- fitdistrplus
- edgeR (install with BiocManager::install('edgeR'))


###############################################################################
## Test run

Please test the package using this mock dataset: fenora/data/test_counts.tsv

Example command for the script fenora.R:
Rscript ./scripts/fenora.R -o ./out -f test_run -c ./data/test_counts.tsv

For more info, run:
Rscript ./scripts/fenora.R --help


###############################################################################
# References

Thiecke, M.J., Wutz, G., Muhar, M., Tang, W., Bevan, S., Malysheva, V., 
Stocsits, R., Neumann, T., Zuber, J., Fraser, P. and Schoenfelder, S., 2020. 
Cohesin-dependent and-independent mechanisms mediate chromosomal contacts 
between promoters and enhancers. Cell reports, 32(3).

Paul F. Harrison (2017). Varistran: Anscombe's variance stabilizing 
transformation for RNA-seq gene expression data. The Journal of Open Source 
Software 2 (16). doi:10.21105/joss.00257

