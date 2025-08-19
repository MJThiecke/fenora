#' Preprocess ChIP-Seq counts per restriction fragment (RF)
#'
#' Applies Anscombe variance stabilising transformation, optional quantile normalisation,
#' and generates optional diagnostic plots.
#'
#' @param rmap_file Path to restriction digestion map (columns: chr, start, end, rfID)
#' @param cMatr_file Path to counts per RF matrix (columns: rfID, counts_target1, ...)
#' @param chrInfo_file Path to chromosome info file (columns: chr, len, optional ploidy)
#' @param outDir Output directory for processed data and plots
#' @param minLen Minimum RF length in bp (default: 100)
#' @param maxLen Maximum RF length in bp (default: 50000)
#' @param qNorm Logical; perform quantile normalisation between chromosomes
#' @param buildPlots Logical; generate plots showing normalisation effects
#' @param chrLenQ Number of chromosome length quantiles for plotting
#' @param fixDisp Fixed dispersion parameter for Anscombe transformation; -1 to estimate automatically
#' @param feature_ids Comma-separated sample IDs
#' @param seed Random seed; -1 for no fixed seed
#'
#' @return Invisibly returns processed counts matrix
#' @export
preprocess_counts_per_rf <- function(
  rmap_file,
  cMatr_file,
  chrInfo_file,
  outDir,
  minLen = 100L,
  maxLen = 50000L,
  qNorm = FALSE,
  buildPlots = FALSE,
  chrLenQ = 5L,
  fixDisp = -1,
  feature_ids = NA,
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
  rmap <- fread(rmap_file, header = TRUE)
  rmap[, fid := paste0(chr, '_', start, '_', end)]
  cMatr <- fread(cMatr_file, header = TRUE, sep = "\t")
  
  # If no feature IDs were supplied, attempt to get them from the cMatr header:
  # the trailing columns should contain the feature counts
  if(is.na(feature_ids)){
    feature_ids <- colnames(cMatr)[4:length(colnames(cMatr))]
  } else {
    feature_ids <- unlist(strsplit(feature_ids, ","))
    if(!all(feature_ids %in% colnames(cMatr))){
      stop("Feature IDs provided must match columns in counts matrix")
    }
  }
  
  cMatr[, fid := paste0(chr, '_', start, '_', end)]
  chrInfo <- fread(chrInfo_file, header = TRUE, sep = "\t")

  message("Filtering restriction fragments by length...")
  rf_lengths <- rmap$end - rmap$start + 1
  keep_idx <- rf_lengths >= minLen & rf_lengths <= maxLen
  rmap <- rmap[keep_idx, ]
  
  # Select only features in the rmap
  cMatr <- cMatr[(fid %in% rmap$fid), ]
  
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
  colnames(counts_transformed) <- paste0(colnames(counts_transformed), '_ansc')
  
  # Add transformed counts to cMatr
  cMatr <- cbind(cMatr, counts_transformed)
  
  # Normalise quantiles between chromosomes
  # This is relevant when dealing with data with aneuploid origin
  if (qNorm) {
    message("Performing between-chromosome quantile normalisation...")
    for(current_id in feature_ids){
      
      current_fdat <- cMatr[, c('fid', paste0(current_id, '_ansc')), with = FALSE]
      
      counts_transformed <- normalizeQuantile.betweenChr(
        feature_dat = current_fdat,
        rmap = rmap
      )
      
      # Add the quantile-normalised data to cMatr
      cMatr <- merge.data.table(x = cMatr, y = counts_transformed, by = 'fid')
    }
  }
  
  # Build plots that show the effect of anscombe and quantile transformations
  if (buildPlots) {
    message('Plotting transformation effects')
    
    # The design determines the features to include
    design <- data.table(
      feature = feature_ids,
      anscombe = paste0(feature_ids, '_ansc')
    )
    if(qNorm){
      design[, q_norm := paste0(anscombe, '_qnorm')]
    }
    
    # Build test plot
    plot_transform_effects(
      feature_dat = cMatr,
      design = design,
      fn_out = file.path(outDir, "'example_plot'"))
    
    message('Done')
  }

  message("Saving processed counts...")
  
  output_file <- file.path(outDir, "processed_counts.txt")
  cols_out <- c(
    'chr', 'start', 'end', 
    feature_ids, paste0(feature_ids, '_ansc'))
  if(qNorm){
    cols_out <- c(cols_out, paste0(feature_ids, '_ansc_qnorm'))
  }
  
  fwrite(
    cMatr[, cols_out, with = FALSE], 
    file = output_file, 
    append = FALSE, quote = FALSE, 
    sep = '\t', eol = '\n', na = 'NA', dec = '.', 
    row.names = FALSE, col.names = TRUE, scipen = 999999)
  
  message("Preprocessing complete: ", output_file)
  
}

