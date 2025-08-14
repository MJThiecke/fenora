#==============================================================================
#
# Perform linear regression on Anscombe-transformed ChIP-Seq data
# (score ~ rfLen * count_chr_mean)
# Build a plot for OLS on Anscombe transformed data
# Write studentized residuals to file
#
# Author: Michiel J. Thiecke
# Date: 01/11/2019
#
#==============================================================================

rm(list = ls())

library(data.table)
library(ggplot2)
library(gridExtra)
library(MASS)
library(quantreg)
library(argparser)

debug <- FALSE

if(debug == FALSE){
  aParse <- arg_parser(
    # TODO: check and update description
    description = "Calculates studentised residuals per RF, using OLS regression. Takes the output from preprocess_counts_per_rf.R. Regresses the preprocessed counts per RF vs the RF length. Outputs the studentised residuals.",
    name = "preprocess_counts_per_RF.R"
  )
  aParse <- add_argument(parser = aParse, arg = "-wd", default = ".", help = "Working directory")
  aParse <- add_argument(parser = aParse, arg = "-outDir", default = ".", help = "Output directory")
  aParse <- add_argument(parser = aParse, arg = "-rMap", default = "GRCh37_HindIII.rmap", help = "The restriction digestion map. Columns: chr, start, end, rfID")
  aParse <- add_argument(parser = aParse, arg = "-pcMatr", default = "counts_matrix_anscombe_qnorm.txt", help = "Preprocessed counts per RF matrix. Columns: rfID, counts_target1, ... counts_targetN")
  aParse <- add_argument(parser = aParse, arg = "-chrInfo", help = "A matrix with two required columns (chr and len for chromosome ID and chromosome length) and one optional column (ploidy). Note that the user should only include relevant chromosomes (e.g. HeLa is female so the Y-chromosome is excluded). Mitochondiral and nuclear DNA should not be analysed simultaneously while enabling between-chromosome quantile normalisation.")
  aParse <- add_argument(parser = aParse, arg = "-buildPlots", default = FALSE, help = "Set to TRUE to build figures that show the effect of between-chromosome quantile normalisation")
  aParse <- add_argument(parser = aParse, arg = "-chrLenQ", default = 5, help = "Set the number of chromosome length quantiles for plotting")
  aParse <- add_argument(parser = aParse, arg = "-exIDs", default = "3", help = "Examples for plotting. May be sample IDs or the number of samples to plot (randomly selected)")
  aParse <- add_argument(parser = aParse, arg = "-seed", default = -1, help = "If not -1, the value is passed to base::set.seed(). Must be positive integer")
  #aParse <- add_argument(parser = aParse, arg = "-", help = "")
  # Get the command line arguments from aParse
  args <- parse_args(aParse)
  setwd(args$wd)
} else {
  args <- list()
  args$wd <- "."
  args$outDir <- "/bi/group/sysgen/michiel/ChIPseq/cleanup/output/"
  args$rMap <- "GRCh37_HindIII.rmap"
  args$pcMatr <- "counts_matrix_anscombe_qnorm.txt"
  args$chrInfo <- "chromosome_info.txt" 
  args$buildPlots <- TRUE
  args$chrLenQ <- 5
  args$exIDs <- c("POLR2A,SCC1,YY1,H3M4me1,H3K27ac,H3K4me3,H3K9me3")
  #args$exIDs <- "2"
  args$seed <- 42
  setwd("/bi/group/sysgen/michiel/ChIPseq/cleanup/")
}

# Check types
stopifnot(is.character(args$wd) & file.exists(args$wd))
stopifnot(is.character(args$outDir) & file.exists(args$outDir))
stopifnot(is.character(args$rMap))
stopifnot(is.character(args$pcMatr))
stopifnot(is.character(args$chrInfo))
stopifnot(is.logical(args$buildPlots))
stopifnot(is.numeric(args$chrLenQ))
if(is.character(args$exIDs)){
  if(!is.na(suppressWarnings(as.integer(args$exIDs)))){#If args$exIDs can be cast to numeric
    args$exIDs <- as.integer(args$exIDs)
    stopifnot(args$exIDs > 0)
  }
} else{
  stop("exIDs must be a positive integer or a comma separated list of sample IDs")
}
if((args$seed != -1) & is.numeric(args$seed) & args$seed > 0){
  args$seed <- as.integer(args$seed)
} else if(args$seed == -1){
  args$seed <- NA
} else{
  stop("seed must be -1, or a positive integer")
}

source("functions/cdc_functions.R")

# Parse chromosome information
chr_info <- fread(paste0("input/", args$chrInfo))

counts_matrix_preproc <- fread(paste0("output/", args$pcMatr))

rmap <- fread(paste0("input/", args$rMap), col.names = c("chr", "start", "end", "rfID"), key = "rfID")
rmap[, len := (end - start)]

# Discard RFs from the rmap that are absent in the preprocessed counts matrix
rmap <- rmap[(rmap$rfID %in% counts_matrix_preproc$rfID), ]

# Warn the user if the RF lengths are not of recommended size
if(rmap[, min(len)] < 100){
  message("The minimum restriction fragment length is below 100bp. It is recommended to set the lower bound to 50bp or 100bp.")
}
if(rmap[, max(len)] > 50000){
  message("The maximum restriction fragment length is greater than 50000bp. It is recommended to set the upper bound to 50000bp.")
}

# Warn the user if the set of RFs have mitochondial as well as nuclear origin
if(any(grepl("MT|Mt|mt|M|m", rmap[, unique(chr)])) & any(grepl("X|x|\\d{1-2}", rmap[, (unique(chr))]))){
  message("The provided restriction fragments appear to have mitochondrial, as well as nuclear origin. It is not recommended to analyse these simultaneously.\nProceed only if you expressly want to analyse mitochondrial and nuclear data simultaneously!")
  if(debug == FALSE){
    message("Discard mitochondiral restriction fragments?")
    usrInput_dm <- promptUser_yn()
    while(is.na(usrInput_dm)){
      usrInput_dm <- promptUser_yn()
    }
    if(usrInput_dm == TRUE){
      rmap <- rmap[(!(grepl("MT|Mt|mt|M|m", rmap[, (chr)]))), ]
    } else if(usrInput_dm == FALSE){
      message("Discard nuclear restriction fragments?")
      usrInput_dn <- promptUser_yn()
      while(is.na(usrInput_dn)){
        usrInput_dn <- promptUser_yn()
      }
      if(usrInput_dn == TRUE){
        rmap <- rmap[((grepl("MT|Mt|mt|M|m", rmap[, (chr)]))), ]
      } else if(usrInput_dn == FALSE & usrInput_dm == FALSE){
        message("Proceeding with mitochondrial as well as nuclear restriction fragments...")
      }
    } 
  }
}

# If plotting is requested, set seed and prepare example figure sample IDs
if(args$buildPlots == TRUE){
  # Get a random set of rows for plotting
  if(!is.na(args$seed)){
    set.seed(args$seed)
  }
  sel <- sample(1:rmap[, .N], size = 30000, replace = FALSE)  
  # Get sample names for plotting
  if(is.character(args$exIDs)){
    current_ids_all <- colnames(counts_matrix_preproc)[2:ncol(counts_matrix_preproc)]
    # Split on comma
    current_elems <- strsplit(args$exIDs, ",")[[1]]
    if(length(current_elems) == 0){
      stop("Illegal argument provided for exIDs. This should be an integer or a comma separated string")
    }
    # Build grep pattern
    current_pattern <- paste0(current_elems, collapse = "|")
    # Overwrite exIDs with matched sample IDs
    args$exIDs <- current_ids_all[grep(current_pattern, current_ids_all)]
    if(length(args$exIDs) == 0){
      stop("None of the sample IDs provided in exIDs matched the samples (columns) in cMatr")
    }
    # Clean up
    rm(list = ls()[grep("current_", ls())])
  } else if(is.numeric(args$exIDs)){
    args$exIDs <- sample(colnames(counts_matrix_preproc)[2:ncol(counts_matrix_preproc)], args$exIDs)
  } else {
    stop("Illegal argument provided for exIDs. This should be an integer or a comma separated string")
  }
}

#------------------------------------------------------------------------------
# Build a chromosome length quantile table and get colours per chromosome based 
# on the chromosome length and ploidy

chrLenQ <- as.integer(args$chrLenQ)

# Build a chromosome-length-quantile table
lenQuantiles <- quantile(chr_info$len, seq(0, 1, by = 1/chrLenQ))
qTable <- matrix(nrow = chrLenQ, ncol = 2, data = NA)
qTable[1:chrLenQ, 1] <- lenQuantiles[1:chrLenQ]
qTable[1:chrLenQ, 2] <- lenQuantiles[2:(chrLenQ + 1)]
colnames(qTable) <- c("start", "end")

myPalette <- colorRampPalette(c("turquoise", "grey27", "goldenrod"))

# Colour by ploidy
ploidyColours <- myPalette(chr_info[, max(ploidy)])
chr_info[, colour_ploidy := ploidyColours[ploidy]]

# Colour by length quantile
chr_info[, lenQ := apply(chr_info[, .(len)], 1, getQuantile, qTable)]
lenQColours <- myPalette(chr_info[, max(lenQ)])
chr_info[, colour_lenQ := lenQColours[lenQ]]
chr_info$lenQ <- factor(chr_info$lenQ, levels = 1:chrLenQ)

rm(lenQuantiles, myPalette, ploidyColours, lenQColours)

# TODO: re-implement this to provide a comparison between OLS, nB and Pois regression
# #------------------------------------------------------------------------------
# # Build a scatter plot with counts vs RF length
# # Draw the lines for OLS, NB and Poisson
# # RF length cutoff hardcoded to 100 because too small RFs create zero inflation problems
# 
# # Perform NB, Pois and OLS regression on selected raw counts
# for(current_sampleID in args$exIDs){
#   message(current_sampleID)
#   counts_sel <- counts.ip[, c(1, grep(current_sampleID[1], colnames(counts.ip))), with = FALSE]
#   # set the minimum RF length to 100 (because negative binomial regression fails on too strongly zero-inflated data)
#   setkey(counts_sel, rfID)
#   counts_sel <- counts_sel[(rmap), ][(len >= 100), c(1, 2), with = FALSE]
# 
#   plotThreeRegressions(
#     scorePerRF = counts_sel, 
#     rmap = rmap, 
#     outDir = paste0(args$outDir, "scatterplots/"), 
#     sel = 1:rmap[, .N]
#   )
# }
# rm(list = ls()[grep("\\.sel", ls())])
# rm(list = ls()[grep("current\\.", ls())])
# 
# #------------------------------------------------------------------------------
# # Run OLS, Poisson and Negative binomial regression on raw count data and build 
# # scatter plots of the residuals
# # RF length cutoff hardcoded to 100 because too small RFs create zero inflation problems
# 
# # Perform NB, Pois and OLS regression on selected raw counts
# for(current_sampleID in args$exIDs){
#   message(current_sampleID)
#   counts_sel <- counts.ip[, c(1, grep(current_sampleID[1], colnames(counts.ip))), with = FALSE]
#   # set the minimum RF length to 100 (because negative binomial regression fails on too strongly zero-inflated data)
#   setkey(counts_sel, rfID)
#   counts_sel <- counts_sel[(rmap[(len >= 100), rfID]), ]
#   
#   for(current_regrType in c("nb", "ols", "pois")){
#     message(current_regrType)
#     dat <- regressDistance(
#       scorePerRF = counts_sel, 
#       rmap = rmap, 
#       regrType = current_regrType, 
#       isRef = FALSE, 
#       plotResids = TRUE, 
#       outDir = paste0(args$outDir, "scatterplots/"), 
#       exclude = FALSE, 
#       sel = 1:rmap[, .N]
#     )
#   }
# }
# rm(list = ls()[grep("\\.sel", ls())])
# rm(list = ls()[grep("current\\.", ls())])

#------------------------------------------------------------------------------
# Run OLS on pre-processed data (log-transform Rflen) and plot the scores and
# colour them acording to P-value cutoff (<0.05)

if(args$buildPlots == TRUE){
  message("\nBuilding scatter plots (note that these are built on a sub-selection of 30000 restriction fragments)...")
  for(current_sampleID in args$exIDs){
    message(current_sampleID)
    current_scores <- copy(counts_matrix_preproc[, (c(1, grep(current_sampleID, colnames(counts_matrix_preproc)))), with = FALSE])
    colnames(current_scores)[2] <- "score"
    setkey(current_scores, "rfID")
    setkey(rmap, "rfID")
    current_scores <- current_scores[(rmap), ]
    setkey(current_scores, chr)
    setkey(chr_info, chr)
    current_scores <- current_scores[(chr_info[, .(colour_ploidy, colour_lenQ), by = chr]), ]
    
    # Perform ordinary least squares regression on log-transformed data
    current_lm_ols <- lm(score ~ log(len, base = 10), data = current_scores)
    current_lm_ols_summary <- summary(current_lm_ols)
    
    # Get residuals for OLS, Studentize them, perform a one-sided t-test
    current_scores[, studres := studres(current_lm_ols)]
    #current_scores[, pT := pt(studres, df = (.N - 3), lower.tail = FALSE)]
    #current_scores[, pT_adj := p.adjust(pT, method = "fdr")]
    #current_scores[, sig_pT := ifelse(pT_adj < 0.05, "p<0.05", "p>=0.05")]
    #current_scores[, table(sig_pT)]
    
    current_plotTitle <- paste0(
      current_sampleID, 
      " OLS on log-processed data\nslope: ", 
      round(coefficients(current_lm_ols)[2], digits = 3), 
      " (orange), rsq: ",
      round(current_lm_ols_summary$r.squared, digits = 3)
    )
    
    current_plot <- ggplot(current_scores[(sel), .(score, len = log(len, base = 10))], aes(x = len, y = score)) +
      #geom_point(aes(color = current_scores[(sel), colour_ploidy]), alpha = 0.5)+
      geom_point(alpha = 0.5)+
      geom_abline(intercept = coefficients(current_lm_ols)[1], slope = coefficients(current_lm_ols)[2], color = "goldenrod", lwd = 0.8) +
      coord_cartesian(xlim = c(rmap[, log(min(len), base = 10)], rmap[, log(max(len), base = 10)])) +
      labs(
        title = current_plotTitle,
        x = "log10(RFlength(bp))",
        y = "processed score"
      ) + theme(legend.position = "none")
    
    current_plot.lin <- ggplot(current_scores[(sel), .(score, len)], aes(x = len, y = score)) +
      geom_point(alpha = 0.5)+
      coord_cartesian(xlim = c(rmap[, min(len)], rmap[, max(len)])) +
      labs(
        x = "RFlength(bp)",
        y = ""
      ) + theme(legend.position = "top")
    
    pdf(file = paste0(args$outDir, "regressionLine_anscombe_qnorm_logLen.pdf"))
    grid.arrange(current_plot, current_plot.lin, nrow = 1)
    dev.off()
  }
  rm(list = ls()[grep("_sel", ls())])
  rm(list = ls()[grep("current_", ls())])
  message("Done!")
}

#------------------------------------------------------------------------------
# Perform linear regression on all samples and write the studentised residuals 
# to file

message("\nPerforming OLS regression and retrieving studentised residuals...")
suppressWarnings(rm(resids_stud))
#TODO: add progress bar
for(current_sampleID in colnames(counts_matrix_preproc)[2:ncol(counts_matrix_preproc)]){
  
  current_data <- counts_matrix_preproc[, c("rfID", current_sampleID), with = FALSE]
  setkey(current_data, rfID)
  tmp <- regressDistance(
    scorePerRF = current_data, 
    rmap = rmap[, .(chr, start, end, rfID, len = log(len, base = 10))],#Note that the RF length is log-transformed here
    regrType = "ols",
    isRef = FALSE, 
    plotResids = FALSE,
    sel = 1:current_data[, .N]
  )
  
  if(exists("resids_stud") == FALSE){
    resids_stud <- copy(tmp$residualsPerRF[, .(rfID, resid.stud)])
    colnames(resids_stud)[2] <- current_sampleID
  } else {
    setkey(resids_stud, rfID)
    tmp.resids <- (tmp$residualsPerRF)
    setkey(tmp.resids, rfID)
    resids_stud <- resids_stud[(tmp.resids[, .(resid.stud), key = "rfID"]), ]
    colnames(resids_stud)[ncol(resids_stud)] <- current_sampleID
  }
  
  rm(tmp)
  
}
myfwrite(resids_stud, file = paste0(args$outDir, "counts_matrix_anscombe_qnorm_ols-studres.tsv"))
message("Done!")

