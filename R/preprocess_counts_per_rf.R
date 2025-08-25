#' Preprocess ChIP-Seq counts per restriction fragment (RF)
#'
#' Applies Anscombe variance stabilising transformation, optional quantile normalisation,
#' and generates optional diagnostic plots.
#'
#' @param fn_counts Path to counts per interval table (columns: chr, start, end, counts_target1, ..., counts_targetN)
#' @param dir_out Output directory for processed data and plots
#' @param fn_stub Prepended to output file names
#' @param thresh_len_min Minimum RF length in bp (default: 100)
#' @param thresh_len_max Maximum RF length in bp (default: 50000)
#' @param perform_qnorm_btchr Logical; perform quantile normalisation between chromosomes
#' @param plot_transformations Logical; generate plots showing normalisation effects
#' @param write_transf_to_file Logical; write intermediate transformed data to file
#' @param fixed_dispersion Fixed dispersion parameter for Anscombe transformation; -1 to estimate automatically
#' @param feature_ids Comma-separated sample IDs or 'all'
#' @param seed Random seed; -1 for no fixed seed
#'
#' @import fitdistrplus
#' @import data.table
#' @return list of 1)Data table with counts transformations 2) feature_ids
#' @export
preprocess_counts_per_rf <- function(
  fn_counts,
  dir_out,
  fn_stub, 
  thresh_len_min = 100L,
  thresh_len_max = 50000L,
  perform_qnorm_btchr = FALSE,
  plot_transformations = FALSE,
  write_transf_to_file = FALSE,
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
  
  message("Filtering restriction fragments by length...")
  rf_lengths <- dat_counts$end - dat_counts$start
  dat_counts <- dat_counts[(
    rf_lengths >= thresh_len_min & rf_lengths <= thresh_len_max), ]
  
  
  # Estimate dispersion parameter for each feature separately
  # NOTE:
  # consider ZIMB fit (pscl::zeroinfl) in stead of just negative binomial
  for(current_feature in feature_ids){
    
    counts_only <- dat_counts[, current_feature, with = FALSE]
    if (fixed_dispersion == -1) {
      current_fit <- fitdistrplus::fitdist(as.list(counts_only[, 1])[[1]], "nbinom")
      current_disp <- 1/current_fit$estimate['size']
      message(
        'Estimated dispersion on ', 
        current_feature, ': ', 
        round(current_disp, digits = 5))
    } else {
      current_disp <- fixed_dispersion
    }
    
    # message("Performing Anscombe transformation...")
    counts_transformed <- as.data.table(vst(
      counts_only, method = "anscombe.nb", dispersion = current_disp))
    
    setnames(
      counts_transformed, 
      old = current_feature, 
      new = paste0(current_feature, '_ansc'))
    
    # Add transformed counts to dat_counts
    dat_counts <- cbind(dat_counts, counts_transformed)
    
  }
  
  # Normalise quantiles between chromosomes
  # This is relevant when dealing with data with aneuploid origin
  if (perform_qnorm_btchr) {
    message("Performing between-chromosome quantile normalisation...")
    for(current_id in feature_ids){
      current_fdat <- dat_counts[, c(
        'chr', 'start', 'end', 'fid', 
        paste0(current_id, '_ansc')), with = FALSE]
      counts_transformed <- normaliseQuantile_betweenChr(
        feature_dat = current_fdat,
        value_col = paste0(current_id, '_ansc')
      )
      # Add the quantile-normalised data to dat_counts
      dat_counts <- merge.data.table(
        x = dat_counts, y = counts_transformed, by = 'fid')
    }
  }
  
  # Perform quantile normalisation between features
  dat_qnorm <- normaliseQuantile_betweenFeature(
    feature_dat = dat_counts,
    features = feature_ids,
    qnorm_bchr = perform_qnorm_btchr)
  # Add the quantile-normalised data to dat_counts
  dat_counts <- merge.data.table(
    x = dat_counts, y = dat_qnorm, by = 'fid', all.x = TRUE)
  
  # Build plots that show the effect of anscombe and quantile transformations
  if (plot_transformations) {
    message('Plotting transformation effects')
    
    # The design determines the features to include
    design <- data.table(
      feature = feature_ids,
      anscombe = paste0(feature_ids, '_ansc')
    )
    if(perform_qnorm_btchr){
      design[, q_norm_btchr := paste0(anscombe, '_qnorm-btchr')]
      design[, q_norm_btfeat := paste0(anscombe, '_qnorm-btchr-btfeat')]
    } else {
      design[, q_norm_btfeat := paste0(anscombe, '_qnorm-btfeat')]
    }
    
    # Build test plot
    plot_transform_effects(
      feature_dat = dat_counts,
      design = design,
      fn_stub = file.path(dir_out, fn_stub))
  }
  
  if (write_transf_to_file) {
    
    message("Saving processed counts...")
    fn_out <- paste0(fn_stub, '_processed_counts.tsv')
    
    output_file <- file.path(dir_out, fn_out)
    # Remove ugly double fwd slash
    output_file <- gsub('\\/\\/', '/', output_file)
    
    cols_out <- c(
      'chr', 'start', 'end', 
      feature_ids, paste0(feature_ids, '_ansc'))
    if(perform_qnorm_btchr){
      cols_out <- c(cols_out, paste0(feature_ids, '_ansc_qnorm-btchr'))
      cols_out <- c(cols_out, paste0(feature_ids, '_ansc_qnorm-btchr-btfeat'))
    } else {
      cols_out <- c(cols_out, paste0(feature_ids, '_ansc_qnorm-btfeat'))
    }
    
    fwrite_tsv(
      dat_counts[, cols_out, with = FALSE], 
      file = output_file
    )
  }
  
  message("Preprocessing complete")
  if(perform_qnorm_btchr){
    fids_out <- paste0(feature_ids, '_ansc_qnorm-btchr-btfeat')
  } else {
    fids_out <- paste0(feature_ids, '_ansc_qnorm-btfeat')
  }
  return(list(preproc = dat_counts, feature_ids = fids_out))
}
