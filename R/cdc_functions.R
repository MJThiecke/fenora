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
#' @param x data.table
#' @param file output file path
#' @export
myfwrite <- function(x, file){
  data.table::fwrite(x, file,
                     quote = FALSE, sep = "\t", eol = "\n",
                     na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
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
#' @param rmap data.table with fid and chr
#' @param chrs vector of chromosomes to normalise between
#' @param exclude_values boolean or numeric vector length 2
#' @return data.table with normalized counts
#' @import data.table
#' @import aroma.light
#' @export
normalizeQuantile.betweenChr <- function(feature_dat, 
                                         rmap, 
                                         chrs = paste0('chr', c(1:22, "X", "Y")), 
                                         exclude_values = FALSE){
  
  # Modify the feature id, but keep it to place it back later
  orig_feature_id <- colnames(feature_dat)[2]
  colnames(feature_dat)[2] <- "value"
  
  data.table::setkey(feature_dat, fid)
  data.table::setkey(rmap, fid)
  fdat_mrg <- merge.data.table(
    x = feature_dat, y = rmap, by = 'fid', all.x = TRUE)
  
  # Get the desired chromosomes that are present in the feature data
  chrs_sel <- fdat_mrg[, unique(chr)]
  chrs_sel <- chrs_sel[chrs_sel %in% chrs]

  val_per_chr <- list()
  for(current_chr in chrs_sel){
    vals <- fdat_mrg[(chr == current_chr), value]
    if(is.logical(exclude_values)){
      if(exclude_values){
        minVal <- min(vals)
        vals[vals <= minVal] <- NA
      }
    } else if(is.numeric(exclude_values) & length(exclude_values) == 2){
      vals[vals >= exclude_values[1] & vals <= exclude_values[2]] <- NA
    } else {
      stop("exclude_values must be boolean or numeric vector of length 2 (lower- and upper bound)")
    }
    val_per_chr[[current_chr]] <- vals
  }

  val_per_chr.qnorm <- aroma.light::normalizeQuantile(val_per_chr)

  fdat_mrg[, value_qnorm := NA_real_]
  for(current_chr in chrs_sel){
    fdat_mrg[(chr == current_chr), 
             value_qnorm := val_per_chr.qnorm[[current_chr]]]
  }

  fdat_out <- fdat_mrg[, .(fid, value_qnorm)]
  data.table::setnames(fdat_out, old = 'value_qnorm', new = paste0(orig_feature_id, '_qnorm'))
  return(fdat_out)
}

#' Plot the transformed data
#' 
#' @param feature_dat data.table with data columns to plot
#' @param design table with design information
#' @param fn_out output file name
#' @import data.table
#' @import ggplot2
#' @import ggrain
plot_transform_effects <- function(feature_dat,
                                   design,
                                   fn_out,
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
  
  # This warns about coersion. Suppress that
  feature_dat_m <- suppressWarnings(melt(
    data = feature_dat[(chr %in% chr_sel), ], 
    id.vars = c('fid', 'chr'), 
    measure.vars = design,
    variable.name = 'type', value.name = 'value'
    ))
  
  feature_dat_m[, type := factor(type, levels = design)]
  
  plot_out <- ggplot(data = feature_dat_m, aes(x = type, y = value)) +
    geom_rain()
  
  fn_out <- paste0(fn_out, '_', chr_sel, '.pdf')
  pdf(file = fn_out, width = 10, height = 6)
  print(plot_out)
  dev.off()
  
}

#' Perform regression on score per restriction fragment
#'
#' @param scorePerRF data.table with rfID and counts
#' @param rmap data.table with rfID, chr, len
#' @param regrType regression type: "ols", "pois", "nb", or "olslog"
#' @param isRef logical, if TRUE label as reference
#' @param plotResids logical, whether to plot residuals
#' @param outDir output directory
#' @param exclude FALSE or numeric vector of length 2
#' @param sel numeric vector of row indices
#' @return list with residualsPerRF and lm.summary
#' @import data.table
#' @import ggplot2
#' @import MASS
#' @export
regressDistance <- function(scorePerRF,
                            rmap,
                            regrType = "ols",
                            isRef = FALSE,
                            plotResids = FALSE,
                            outDir = getwd(),
                            exclude = FALSE,
                            sel = NA){

  stopifnot(colnames(scorePerRF)[1] == "rfID")
  stopifnot(data.table::key(scorePerRF) == "rfID")
  stopifnot(data.table::key(rmap) == "rfID")
  stopifnot(regrType %in% c("ols", "pois", "nb", "olslog"))

  sampleID <- colnames(scorePerRF)[2]
  colnames(scorePerRF)[2] <- "score"

  data.table::setkey(scorePerRF, rfID)
  scorePerRF <- scorePerRF[(rmap), ][, .(chr, len, score), by = "rfID"]

  if(is.logical(exclude) && exclude){
    scorePerRF <- scorePerRF[score > min(score)]
  } else if(is.numeric(exclude) && length(exclude) == 2){
    scorePerRF <- scorePerRF[score < exclude[1] | score > exclude[2]]
  }

  if(is.na(sel)){
    sel <- sample(1:nrow(scorePerRF), size = min(30000, nrow(scorePerRF)))
  }

  if(regrType == "nb"){
    current.lm <- MASS::glm.nb(score ~ len, data = scorePerRF)
    resid.pearson <- residuals(current.lm, type = "pearson")
    resid.stud <- MASS:studres(current.lm)
  } else if(regrType == "pois"){
    current.lm <- glm(score ~ len, data = scorePerRF, family = poisson(link = "identity"))
    resid.pearson <- residuals(current.lm, type = "pearson")
    resid.stud <- MASS:studres(current.lm)
  } else if(regrType == "ols"){
    current.lm <- lm(score ~ len, data = scorePerRF)
    resid.pearson <- residuals(current.lm)
    resid.stud <- MASS:studres(current.lm)
  } else if(regrType == "olslog"){
    current.lm <- lm(score ~ log(len, base = 10), data = scorePerRF)
    resid.pearson <- residuals(current.lm)
    resid.stud <- MASS:studres(current.lm)
  }

  scorePerRF.resid <- cbind(scorePerRF, data.table(resid.pearson, resid.stud))

  if(plotResids){
    expType <- ifelse(isRef, "Ref", "IP")
    residsPlot <- ggplot2::ggplot(scorePerRF.resid[sel, ], ggplot2::aes(x = len, y = resid.pearson)) +
      ggplot2::geom_point(alpha = 1, size = 2, stroke = 0) +
      ggplot2::coord_cartesian(xlim = c(0, 50000)) +
      ggplot2::labs(title = paste0(sampleID, " ", toupper(regrType), " pearson residuals"),
                    x = "RFlength(bp)", y = "residual")

    ggplot2::ggsave(filename = paste0(outDir, expType, "_", sampleID, "_", regrType, "_pearson_resids.png"), plot = residsPlot, width = 8, height = 8, units = "in", dpi = 300)
    ggplot2::ggsave(filename = paste0(outDir, expType, "_", sampleID, "_", regrType, "_pearson_resids.pdf"), plot = residsPlot)
  }

  list(residualsPerRF = scorePerRF.resid, lm.summary = summary(current.lm))
}

