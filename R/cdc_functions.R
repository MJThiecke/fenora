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
#' @param countData data.table with columns rfID and counts
#' @param rmap data.table with rfID and chr
#' @param chrs vector of chromosomes
#' @param excludeVals boolean or numeric vector length 2
#' @return data.table with normalized counts
#' @import data.table
#' @import aroma.light
#' @export
normalizeQuantile.betweenChr <- function(countData, rmap, chrs = c(1:22, "X", "Y"), excludeVals = FALSE){
  sample.ID <- colnames(countData)[2]
  colnames(countData)[2] <- "value"

  data.table::setkey(countData, rfID)
  data.table::setkey(rmap, rfID)
  countData <- countData[(rmap), ][, .(chr, value), by = "rfID"]

  valPerChr <- list()
  for(current.chr in chrs){
    vals <- countData[(chr == current.chr), value]
    if(is.logical(excludeVals)){
      if(excludeVals){
        minVal <- min(vals)
        vals[vals <= minVal] <- NA
      }
    } else if(is.numeric(excludeVals) & length(excludeVals) == 2){
      vals[vals >= excludeVals[1] & vals <= excludeVals[2]] <- NA
    } else {
      stop("excludeVals must be boolean or numeric vector of length 2")
    }
    valPerChr[[current.chr]] <- vals
  }

  valPerChr.qnorm <- aroma.light::normalizeQuantile(valPerChr)

  countData[, value_qnorm := value]
  for(current.chr in chrs){
    idx <- which(countData$chr == current.chr)
    countData$value_qnorm[idx] <- valPerChr.qnorm[[current.chr]]
    countData$value_qnorm[is.na(countData$value_qnorm)] <- countData$value[is.na(countData$value_qnorm)]
  }

  data.table::setkey(countData, rfID)
  countData[, .(rfID, value, value_qnorm)]
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

