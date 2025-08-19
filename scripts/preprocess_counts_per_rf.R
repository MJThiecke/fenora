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
              help = "Working directory [default: %default]"),
  make_option(c("-o", "--outDir"), type = "character", default = ".",
              help = "Output directory [default: %default]"),
  make_option(c("-r", "--rMap"), type = "character", default = "GRCh37_HindIII.rmap",
              help = "Restriction digestion map file"),
  make_option(c("-c", "--cMatr"), type = "character", default = "aggregated_ip_counts_per_rf.txt",
              help = "Counts per RF matrix file"),
  make_option(c("-i", "--chrInfo"), type = "character",
              help = "Chromosome info file (chr, len, optional ploidy)"),
  make_option("--minLen", type = "integer", default = 100,
              help = "Lower limit of RF length (bp) [default: %default]"),
  make_option("--maxLen", type = "integer", default = 50000,
              help = "Upper limit of RF length (bp) [default: %default]"),
  make_option("--qNorm", type = "logical", default = FALSE,
              help = "Perform quantile normalisation between chromosomes [default: %default]"),
  make_option("--buildPlots", type = "logical", default = FALSE,
              help = "Build figures showing effect of between-chromosome quantile normalisation [default: %default]"),
  make_option("--chrLenQ", type = "integer", default = 5,
              help = "Number of chromosome length quantiles for plotting [default: %default]"),
  make_option("--fixDisp", type = "double", default = -1,
              help = "Fixed dispersion parameter for Anscombe transformation; set -1 to estimate automatically [default: %default]"),
  make_option("--exIDs", type = "character", default = "3",
              help = "Examples for plotting: sample IDs or number of samples (randomly selected) [default: %default]"),
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
    rMap = "data/GRCh37_HindIII.rmap",
    cMatr = "data/test_counts.tsv",
    chrInfo = "data/chromosome_info_default.txt",
    minLen = 100L,
    maxLen = 50000L,
    qNorm = TRUE,
    buildPlots = TRUE,
    chrLenQ = 5L,
    fixDisp = -1,
    exIDs = "h3k4me1,h3k27ac",
    seed = 1337L
  )
}

# Set working directory
setwd(args$wd)

# Run pre-processing
preprocess_counts_per_rf(
  rmap_file     = args$rMap,
  cMatr_file    = args$cMatr,
  chrInfo_file  = args$chrInfo,
  outDir        = args$outDir,
  minLen        = args$minLen,
  maxLen        = args$maxLen,
  qNorm         = args$qNorm,
  buildPlots    = args$buildPlots,
  chrLenQ       = args$chrLenQ,
  fixDisp       = args$fixDisp,
  exIDs         = args$exIDs,
  seed          = args$seed
)

message("Preprocessing complete")
