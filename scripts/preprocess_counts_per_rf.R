#==============================================================================
#
# Pre-processing and plotting of feature counts for OLS regression
# This script is meant to be executed manually. It should not be executed by
# another script.
#
# Author: Michiel J. Thiecke
# Date: 01/11/2019
# Updated: 14/08/2025 - Converted to optparse
#
#==============================================================================

rm(list = ls())

devtools::load_all('/hpc/compgen/users/mthiecke/workspace/fenora')
devtools::document('/hpc/compgen/users/mthiecke/workspace/fenora')

suppressPackageStartupMessages({
  library(optparse)
})

debug <- TRUE

# Define options
option_list <- list(
  make_option(c("-w", "--wd"), type = "character", default = ".",
              help = "Working directory [default: current dir]"),
  make_option(c("-o", "--outDir"), type = "character", default = ".",
              help = "Output directory [default: current dir]"),
  make_option(c("-f", "--fnOut"), type = "character", default = "transf_fig",
              help = "Output file name for transformation figures"),
  make_option(c("-r", "--intervalMap"), type = "character", default = "GRCh37_HindIII.rmap",
              help = "Bed-like file containing the intervals (required columns: chr, start, end)"),
  make_option(c("-c", "--cMatr"), type = "character",
              help = "Counts per interval table (tab-separated, with header)"),
  make_option(c("-i", "--chrInfo"), type = "character",
              help = "Chromosome info file (chr, len, ploidy). Determines which chromosomes are included"),
  make_option("--minLen", type = "integer", default = 100,
              help = "Lower limit of interval length (bp) [default: %default]"),
  make_option("--maxLen", type = "integer", default = 50000,
              help = "Upper limit of interval length length (bp) [default: %default]"),
  make_option("--qNorm", type = "logical", default = FALSE,
              help = "Perform quantile normalisation between chromosomes [default: %default]"),
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
    wd = "/hpc/compgen/users/mthiecke/workspace/fenora/",
    outDir = "/hpc/compgen/users/mthiecke/workspace/fenora_test/out/",
    intervalMap = "data/GRCh37_HindIII.rmap",
    cMatr = "data/test_counts.tsv",
    chrInfo = "data/chromosome_info_default.txt",
    fnOut = "test",
    minLen = 100L,
    maxLen = 50000L,
    qNorm = TRUE,
    buildPlots = TRUE,
    transfToFile = TRUE,
    fixDisp = -1,
    featureIDs = 'all',
    seed = 1337L
  )
}

# Set working directory
setwd(args$wd)

# Run pre-processing
preproc_dat <- preprocess_counts_per_rf(
  fn_counts = args$cMatr,
  fn_chr_info = args$chrInfo,
  dir_out = args$outDir,
  fn_stub = args$fnOut,
  thresh_len_min = args$minLen,
  thresh_len_max = args$maxLen,
  perform_qnorm = args$qNorm,
  plot_transformations = args$buildPlots,
  write_transf_to_file = args$transfToFile,
  fixed_dispersion = args$fixDisp,
  feature_ids = args$featureIDs,
  seed = args$seed
)


# Regress out chromosome length effects
# TODO
test <- regressDistance(
  score_per_frag = preproc_dat$preproc, 
  feature_ids = preproc_dat$feature_ids,
  rmap = args$intervalMap, 
  fn_stub = args$fnOut,
  regrType = "ols",
  plotResids = TRUE,
  outDir = args$outDir
  )
  
  
  
  