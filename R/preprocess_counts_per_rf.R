#' Preprocess ChIP-Seq counts per restriction fragment (RF)
#'
#' Applies Anscombe variance stabilising transformation, optional quantile normalisation,
#' and generates optional diagnostic plots.
#'
#' @param rMap_file Path to restriction digestion map (columns: chr, start, end, rfID)
#' @param cMatr_file Path to counts per RF matrix (columns: rfID, counts_target1, ...)
#' @param chrInfo_file Path to chromosome info file (columns: chr, len, optional ploidy)
#' @param outDir Output directory for processed data and plots
#' @param minLen Minimum RF length in bp (default: 100)
#' @param maxLen Maximum RF length in bp (default: 50000)
#' @param qNorm Logical; perform quantile normalisation between chromosomes
#' @param buildPlots Logical; generate plots showing normalisation effects
#' @param chrLenQ Number of chromosome length quantiles for plotting
#' @param fixDisp Fixed dispersion parameter for Anscombe transformation; -1 to estimate automatically
#' @param exIDs Comma-separated sample IDs or integer number of samples for example plots
#' @param seed Random seed; -1 for no fixed seed
#'
#' @return Invisibly returns processed counts matrix
#' @export
preprocess_counts_per_rf <- function(
  rMap_file,
  cMatr_file,
  chrInfo_file,
  outDir,
  minLen = 100L,
  maxLen = 50000L,
  qNorm = FALSE,
  buildPlots = FALSE,
  chrLenQ = 5L,
  fixDisp = -1,
  exIDs = "3",
  seed = -1
) {
  # Ensure output directory exists
  if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = TRUE)
  }

  if (seed != -1) {
    set.seed(as.integer(seed))
  }

  message("Reading input files...")
  rMap <- fread(rMap_file, header = TRUE)
  rMap[, fid := paste0(chr, '_', start, '_', end)]
  cMatr <- fread(cMatr_file, header = TRUE, sep = "\t")
  # Get the feature names from the cMatr header
  # The trailing columns must contain the feature counts
  feature_ids <- colnames(cMatr)[4:length(colnames(cMatr))]
  
  cMatr[, fid := paste0(chr, '_', start, '_', end)]
  chrInfo <- fread(chrInfo_file, header = TRUE, sep = "\t")

  message("Filtering restriction fragments by length...")
  rf_lengths <- rMap$end - rMap$start + 1
  keep_idx <- rf_lengths >= minLen & rf_lengths <= maxLen
  rMap <- rMap[keep_idx, ]
  
  # Select only features in the rmap
  cMatr <- cMatr[(fid %in% rMap$fid), ]
  
  message("Performing Anscombe transformation...")
  counts_only <- cMatr[, feature_ids, with = FALSE]
  print(counts_only)
  if (fixDisp == -1) {
    # Estimate dispersion (placeholder - replace with your actual estimation method)
    disp <- edgeR::estimateDisp(counts_only)$common.dispersion
  } else {
    disp <- fixDisp
  }
  
  counts_transformed <- vst(counts_only, method = "anscombe.nb", dispersion = disp)
  
  print(counts_transformed)
  stop('intentional')

  if (qNorm) {
    message("Performing between-chromosome quantile normalisation...")
    counts_transformed <- normalizeQuantile.betweenChr(
      counts_transformed,
      chrInfo = chrInfo,
      rMap = rMap
    )
  }

  if (buildPlots) {
    message("Generating diagnostic plots...")
    plot_between_chr_norm(
      counts_transformed,
      chrInfo = chrInfo,
      rMap = rMap,
      chrLenQ = chrLenQ,
      exIDs = exIDs
    )
  }

  message("Saving processed counts...")
  processed <- cbind(rfID = cMatr$rfID, counts_transformed)
  output_file <- file.path(outDir, "processed_counts.txt")
  write.table(processed, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

  message("Preprocessing complete: ", output_file)
  invisible(processed)
}

