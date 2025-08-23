#' Get the row from a quantile table to which a value belongs
#'
#' @param val numeric value
#' @param qTable data.frame or data.table with a column "end"
#' @return integer row number
#' @export
getQuantile <- function(val, qTable){
  min(which(qTable[, "end"] >= val))
}

#' Write data.table to file with preferred settings
#'
#' @param x a data.table object
#' @param file output file path
#' @export
fwrite_tsv <- function(x, file){
  data.table::fwrite(
    x, file, 
    append = FALSE, quote = FALSE, 
    sep = "\t", eol = "\n", na = "NA", dec = ".", 
    row.names = FALSE, col.names = TRUE, scipen = 999999)
}

#' Prompt user for a yes/no answer
#'
#' @return TRUE, FALSE, or NA
#' @export
promptUser_yn <- function(){
  if(interactive()){
    usrInput <- readline(prompt = "(y/n) ")
  } else {
    message("(y/n) ")
    usrInput <- readLines("stdin", n = 1)
  }

  if((grepl("^[Yy]$", usrInput)) | tolower(usrInput) == "yes"){
    TRUE
  } else if((grepl("^[Nn]$", usrInput)) | tolower(usrInput) == "no"){
    FALSE
  } else {
    NA
  }
}

#' Quantile normalize count data between chromosomes
#'
#' @param feature_dat data.table with columns fid and counts
#' @param chrs vector of chromosomes to normalise between
#' @param value_col name of the column in feature_dat that contains the values to normalize
#' @return data.table with normalized counts
#' @import data.table
#' @importFrom aroma.light normalizeQuantile
#' @export
normalizeQuantile.betweenChr <- function(feature_dat, 
                                         chrs = paste0('chr', c(1:22, "X", "Y")),
                                         value_col
                                         ){
  if(length(chrs) <= 1){
    stop("Cannot perform quantile normalization with fewer than two chromosomes")
  }
  if(!value_col %in% colnames(feature_dat)){
    stop("The value column ", value_col, " was not in the feature_dat columns")
  }
  
  # Get the desired chromosomes that are present in the feature data
  chrs_sel <- feature_dat[, unique(chr)]
  chrs_sel <- chrs_sel[chrs_sel %in% chrs]
  
  val_per_chr <- list()
  for(current_chr in chrs_sel){
    vals <- as.list(feature_dat[(chr == current_chr), 
                                value_col, with = FALSE])[[1]]
    val_per_chr[[current_chr]] <- vals
  }
  
  # Perform quantile normalisation
  val_per_chr.qnorm <- aroma.light::normalizeQuantile(val_per_chr)

  feature_dat[, value_qnorm := NA_real_]
  for(current_chr in chrs_sel){
    feature_dat[(chr == current_chr), 
             value_qnorm := val_per_chr.qnorm[[current_chr]]]
  }

  fdat_out <- feature_dat[, .(fid, value_qnorm)]
  data.table::setnames(fdat_out, 
                       old = 'value_qnorm', 
                       new = paste0(value_col, '_qnorm'))
  return(fdat_out)
}

#' Plot the transformed data
#' 
#' @param feature_dat data.table with data columns to plot
#' @param design table with design information
#' @param fn_stub prepended to output file name
#' @param chr_sel example chromosome to plot
#' @import data.table
#' @import ggplot2
#' @import ggrain
plot_transform_effects <- function(feature_dat,
                                   design,
                                   fn_stub,
                                   chr_sel = 'chr11'){
  
  design[, c(colnames(design)), with = FALSE]
  design <- as.list(design)
  names(design) <- NULL
  design <- unlist(design)
  
  if(!chr_sel %in% feature_dat$chr){
    message("Did not find ", chr_sel, " in feature_dat")
    chr_sel <- feature_dat[, unique(chr)][1]
    message("Selecting ", chr_sel, " for plotting")
  }
  
  # This warns about coercion. Suppress that
  feature_dat_m <- suppressWarnings(melt(
    data = feature_dat[(chr %in% chr_sel), ], 
    id.vars = c('fid', 'chr'), 
    measure.vars = design,
    variable.name = 'type', value.name = 'value'
    ))
  
  feature_dat_m[, type := factor(type, levels = design)]
  
  plot_out <- ggplot(data = feature_dat_m, aes(x = type, y = value)) +
    geom_rain()
  
  fn_out <- paste0(fn_stub, '_transf_effects_', chr_sel, '.pdf')
  
  pdf(file = fn_out, width = 10, height = 6)
  print(plot_out)
  dev.off()
  
}

#' Perform regression on score per restriction fragment
#'
#' @param score_per_frag data.table with preprocessed counts per fragment
#' @param rmap data.table with rfID, chr, len
#' @param regrType regression type: "ols", "pois", "nb", or "olslog"
#' @param isRef logical, if TRUE label as reference
#' @param plotResids logical, whether to plot residuals
#' @param outDir output directory
#' @param include FALSE or numeric vector of length 2: include scores in interval
#' @param sel numeric vector of row indices
#' @return list with residualsPerRF and lm.summary
#' @import data.table
#' @import ggplot2
#' @import MASS
#' @export
regressDistance <- function(score_per_frag,
                            feature_ids,
                            rmap,
                            regrType = "ols",
                            isRef = FALSE,
                            plotResids = FALSE,
                            outDir = getwd(),
                            include = FALSE,
                            sel = NA){
  stopifnot(regrType %in% c("ols", "pois", "nb", "olslog"))
  
  # data.table::setkey(score_per_frag, fid)
  # sampleID <- colnames(score_per_frag)[2]
  for(current_feature in feature_ids){
    
    data.table::setnames(score_per_frag, old = current_feature, new = 'score')
    
    if(is.numeric(include) && length(include) == 2){
      score_per_frag <- score_per_frag[(!(score < include[1] | score > include[2])), ]
    }
    
    browser()
    
    
    
    data.table::setnames(score_per_frag, old = 'score', new = current_feature)
  }
  
  
  
  

  if(is.na(sel)){
    sel <- sample(1:nrow(score_per_frag), size = min(30000, nrow(score_per_frag)))
  }

  if(regrType == "nb"){
    current.lm <- MASS::glm.nb(score ~ len, data = score_per_frag)
    resid.pearson <- residuals(current.lm, type = "pearson")
    resid.stud <- MASS:studres(current.lm)
  } else if(regrType == "pois"){
    current.lm <- glm(score ~ len, data = score_per_frag, family = poisson(link = "identity"))
    resid.pearson <- residuals(current.lm, type = "pearson")
    resid.stud <- MASS:studres(current.lm)
  } else if(regrType == "ols"){
    current.lm <- lm(score ~ len, data = score_per_frag)
    resid.pearson <- residuals(current.lm)
    resid.stud <- MASS:studres(current.lm)
  } else if(regrType == "olslog"){
    current.lm <- lm(score ~ log(len, base = 10), data = score_per_frag)
    resid.pearson <- residuals(current.lm)
    resid.stud <- MASS:studres(current.lm)
  }

  score_per_frag.resid <- cbind(score_per_frag, data.table(resid.pearson, resid.stud))

  if(plotResids){
    expType <- ifelse(isRef, "Ref", "IP")
    residsPlot <- ggplot2::ggplot(score_per_frag.resid[sel, ], ggplot2::aes(x = len, y = resid.pearson)) +
      ggplot2::geom_point(alpha = 1, size = 2, stroke = 0) +
      ggplot2::coord_cartesian(xlim = c(0, 50000)) +
      ggplot2::labs(title = paste0(sampleID, " ", toupper(regrType), " pearson residuals"),
                    x = "RFlength(bp)", y = "residual")

    ggplot2::ggsave(filename = paste0(outDir, expType, "_", sampleID, "_", regrType, "_pearson_resids.png"), plot = residsPlot, width = 8, height = 8, units = "in", dpi = 300)
    ggplot2::ggsave(filename = paste0(outDir, expType, "_", sampleID, "_", regrType, "_pearson_resids.pdf"), plot = residsPlot)
  }

  list(residualsPerRF = score_per_frag.resid, lm.summary = summary(current.lm))
}

