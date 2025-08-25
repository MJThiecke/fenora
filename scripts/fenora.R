#==============================================================================
#
# Pre-processing feature counts in genomic intervals of arbitrary length.
#
# Author: Michiel J. Thiecke
#
#==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(optparse)
  library(fenora)
})

debug <- FALSE

# Define options
option_list <- list(
  make_option(c("-o", "--outDir"), type = "character", default = ".",
              help = "Output directory [default: current dir]"),
  make_option(c("-f", "--fnOut"), type = "character", default = "fenora",
              help = "Output file name for transformation figures"),
  make_option(c("-c", "--cMatr"), type = "character",
              help = "Counts per interval table (tab-separated, with header)"),
  make_option("--minLen", type = "integer", default = 100,
              help = "Lower limit of interval length (bp) [default: %default]"),
  make_option("--maxLen", type = "integer", default = 50000,
              help = "Upper limit of interval length length (bp) [default: %default]"),
  make_option("--qNormBetweenChr", type = "logical", default = FALSE,
              help = "Perform quantile normalisation between chromosomes. This is usefull when dealing with data from aneuploid origin. [default: %default]"),
  make_option("--buildPlots", type = "logical", default = FALSE,
              help = "Build figures showing effect of transformations [default: %default]"),
  make_option("--transfToFile", type = "logical", default = FALSE,
              help = "Write transformed data to file [default: %default]"),
  make_option("--fixDisp", type = "double", default = -1,
              help = "Fixed dispersion parameter for Anscombe transformation; set -1 to estimate automatically [default: %default]"),
  make_option("--featureIDs", type = "character", default = 'all',
              help = "The features to select from the counts matrix (comma-separated). If 'all', all features are used. [default: all]"),
  make_option("--seed", type = "integer", default = -1,
              help = "Random seed; set -1 for no fixed seed [default: %default]")
)

if (!debug) {
  parser <- OptionParser(option_list = option_list,
    description = "Pre-processing of ChIP-Seq counts per restriction fragment (RF). Uses Anscombe variance stabilising transformation and can optionally perform quantile normalisation to mitigate variance caused by aneuploidy."
  )
  args <- parse_args(parser)
} else {
  args <- list(
    outDir = "[your_dir]",
    cMatr = "[your_dir]/fenora/data/test_counts.tsv",
    fnOut = "test",
    minLen = 100,
    maxLen = 50000,
    qNormBetweenChr = TRUE,
    buildPlots = TRUE,
    transfToFile = TRUE,
    fixDisp = -1,
    featureIDs = 'all',
    seed = 1337
  )
}

# Run pre-processing
preproc_dat <- preprocess_counts_per_rf(
  fn_counts = args$cMatr,
  dir_out = args$outDir,
  fn_stub = args$fnOut,
  thresh_len_min = args$minLen,
  thresh_len_max = args$maxLen,
  perform_qnorm_btchr = args$qNormBetweenChr,
  plot_transformations = args$buildPlots,
  write_transf_to_file = args$transfToFile,
  feature_ids = args$featureIDs,
  seed = args$seed
)

# Regress out chromosome length effects
regressDistance(
  score_per_frag = preproc_dat$preproc, 
  feature_ids = preproc_dat$feature_ids,
  fn_stub = args$fnOut,
  regrType = "ols",
  plotResids = args$buildPlots,
  outDir = args$outDir
  )
  