#' Preprocess ChIP-Seq counts per restriction fragment (RF)
#'
#' Applies Anscombe variance stabilising transformation, optional quantile normalisation,
#' and generates optional diagnostic plots.
#'
#' @param rmap_file Path to restriction digestion map (columns: chr, start, end)
#' @param fn_counts Path to counts per interval table (columns: chr, start, end, counts_target1, ..., counts_targetN)
#' @param fn_chr_info Path to chromosome info file (columns: chr, len, optional ploidy)
#' @param dir_out Output directory for processed data and plots
#' @param thresh_len_min Minimum RF length in bp (default: 100)
#' @param thresh_len_max Maximum RF length in bp (default: 50000)
#' @param perform_qnorm Logical; perform quantile normalisation between chromosomes
#' @param plot_transformations Logical; generate plots showing normalisation effects
#' @param fixed_dispersion Fixed dispersion parameter for Anscombe transformation; -1 to estimate automatically
#' @param feature_ids Comma-separated sample IDs
#' @param seed Random seed; -1 for no fixed seed
#'
#' @importFrom edgeR estimateDisp
#' @import data.table
#' @return Void
#' @export
preprocess_counts_per_rf <- function(
  rmap_file,
  fn_counts,
  fn_chr_info,
  dir_out,
  thresh_len_min = 100L,
  thresh_len_max = 50000L,
  perform_qnorm = FALSE,
  plot_transformations = FALSE,
  fixed_dispersion = -1,
  feature_ids = 'all',
  seed = -1
) {
  
  # Ensure output directory exists
  if (!dir.exists(dir_out)) {
    dir.create(dir_out, recursive = TRUE)
  }

  if (seed != -1) {
    set.seed(as.integer(seed))
  }

  message("Reading input files...")
  rmap <- fread(rmap_file, header = TRUE)
  rmap[, fid := paste0(chr, '_', start, '_', end)]
  dat_counts <- fread(fn_counts, header = TRUE, sep = "\t")
  
  # If all feature IDs were requested, get them from the dat_counts header:
  # the trailing columns should contain the feature counts
  if(feature_ids == 'all'){
    message('Using all feature IDs from counts matrix...')
    feature_ids <- colnames(dat_counts)[4:length(colnames(dat_counts))]
    message('Found feature IDs: ', paste0(feature_ids, collapse = ', '))
  } else {
    feature_ids <- unlist(strsplit(feature_ids, ","))
    if(!all(feature_ids %in% colnames(dat_counts))){
      stop("Feature IDs provided must match columns in counts matrix")
    }
    message('Using feature IDs: ', paste0(feature_ids, collapse = ', '))
  }
  
  dat_counts[, fid := paste0(chr, '_', start, '_', end)]
  dat_chr_info <- fread(fn_chr_info, header = TRUE, sep = "\t")
  
  # Select only chromosomes that are present in dat_chr_info
  dat_counts <- dat_counts[(chr %in% dat_chr_info$chr), ]
  
  message("Filtering restriction fragments by length...")
  rf_lengths <- rmap$end - rmap$start
  rmap <- rmap[(rf_lengths >= thresh_len_min & rf_lengths <= thresh_len_max), ]
  
  # Select only features in the rmap
  dat_counts <- dat_counts[(fid %in% rmap$fid), ]
  
  counts_only <- dat_counts[, feature_ids, with = FALSE]
  if (fixed_dispersion == -1) {
    message('Estimating dispersion using EdgeR common dispersion...')
    disp <- suppressMessages(edgeR::estimateDisp(counts_only)$common.dispersion)
    message('Estimated dispersion: ', disp)
  } else {
    disp <- fixed_dispersion
  }
  
  message("Performing Anscombe transformation...")
  counts_transformed <- vst(counts_only, method = "anscombe.nb", dispersion = disp)
  colnames(counts_transformed) <- paste0(colnames(counts_transformed), '_ansc')
  # Add transformed counts to dat_counts
  dat_counts <- cbind(dat_counts, counts_transformed)
  
  # Normalise quantiles between chromosomes
  # This is relevant when dealing with data with aneuploid origin
  if (perform_qnorm) {
    message("Performing between-chromosome quantile normalisation...")
    for(current_id in feature_ids){
      current_fdat <- dat_counts[, c('fid', paste0(current_id, '_ansc')), with = FALSE]
      counts_transformed <- normalizeQuantile.betweenChr(
        feature_dat = current_fdat,
        rmap = rmap
      )
      # Add the quantile-normalised data to dat_counts
      dat_counts <- merge.data.table(x = dat_counts, y = counts_transformed, by = 'fid')
    }
  }
  
  # Build plots that show the effect of anscombe and quantile transformations
  if (plot_transformations) {
    message('Plotting transformation effects')
    
    # The design determines the features to include
    design <- data.table(
      feature = feature_ids,
      anscombe = paste0(feature_ids, '_ansc')
    )
    if(perform_qnorm){
      design[, q_norm := paste0(anscombe, '_qnorm')]
    }
    
    # Build test plot
    plot_transform_effects(
      feature_dat = dat_counts,
      design = design,
      fn_out = file.path(dir_out, "'example_plot'"))
  }

  message("Saving processed counts...")
  
  output_file <- file.path(dir_out, "processed_counts.txt")
  # Remove ugly double fwd slash
  output_file <- gsub('\\/\\/', '/', output_file)
  
  cols_out <- c(
    'chr', 'start', 'end', 
    feature_ids, paste0(feature_ids, '_ansc'))
  if(perform_qnorm){
    cols_out <- c(cols_out, paste0(feature_ids, '_ansc_qnorm'))
  }
  
  fwrite(
    dat_counts[, cols_out, with = FALSE], 
    file = output_file, 
    append = FALSE, quote = FALSE, 
    sep = '\t', eol = '\n', na = 'NA', dec = '.', 
    row.names = FALSE, col.names = TRUE, scipen = 999999)
  
  message("Preprocessing complete: ", output_file)
  
}
