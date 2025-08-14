#' @importFrom grDevices axisTicks dev.off pdf png postscript
#' @importFrom stats as.dendrogram dist hclust is.leaf optimize order.dendrogram sd
#' @importFrom graphics par plot.new
#' @importFrom MASS Null
#' @importFrom edgeR calcNormFactors

# Internal helper function for asinh/log2-like transformation
#' @keywords internal
.log2ish_asinh <- function(a, b) {
    ab <- a * b
    log2(ab + sqrt(1 + ab * ab)) - log2(b) - 1
}

# List of supported VST methods
#' @keywords internal
vst_methods <- list(
    naive.poisson = list(
        description      = "Naive Poisson Variance Stabilizing Transformation.",
        units            = "sqrt",
        is.logish        = FALSE,
        needs.dispersion = FALSE,
        vst              = function(x) sqrt(x)
    ),
    naive.nb = list(
        description      = "Naive negative binomial Variance Stabilizing Transformation.",
        units            = "log2",
        is.logish        = TRUE,
        needs.dispersion = TRUE,
        vst              = function(x, dispersion)
            .log2ish_asinh(sqrt(x), sqrt(dispersion)) * 2.0
    ),
    anscombe.poisson = list(
        description      = "Anscombe's Variance Stabilizing Transformation for the Poisson distribution.",
        units            = "sqrt",
        is.logish        = FALSE,
        needs.dispersion = FALSE,
        vst              = function(x) sqrt(x + 0.375)
    ),
    anscombe.nb.simple = list(
        description      = "Anscombe's simplified Variance Stabilizing Transformation for the negative binomial distribution.",
        units            = "log2",
        is.logish        = TRUE,
        needs.dispersion = TRUE,
        vst              = function(x, dispersion) log2(x + 0.5 / dispersion)
    ),
    anscombe.nb = list(
        description      = "Anscombe's Variance Stabilizing Transformation for the negative binomial distribution.",
        units            = "log2",
        is.logish        = TRUE,
        needs.dispersion = TRUE,
        vst              = function(x, dispersion)
            .log2ish_asinh(sqrt(x + 0.375),
                            sqrt(1 / (1 / dispersion - 0.75))) * 2.0
    )
)

# Internal helper: calculate dispersion score
#' @keywords internal
.dispersion_score <- function(mat, design = NULL) {
    if (is.null(design)) design <- matrix(1, ncol = 1, nrow = ncol(mat))
    residuals <- mat %*% MASS::Null(design)
    rsd <- sqrt(rowMeans(residuals * residuals))
    sd(rsd) / mean(rsd)
}

# Internal helper: estimate optimal dispersion
#' @keywords internal
.optimal_dispersion <- function(x, method = "anscombe.nb", lib.size = NULL, design = NULL) {
    x <- x[rowMeans(x) >= 5,, drop = FALSE]
    if (nrow(x) == 0 || ncol(x) == 0) {
        warning("Insufficient data to estimate dispersion, default value used.")
        return(1.0)
    }
    optimize(
        function(d) {
            .dispersion_score(
                vst(x, method = method, lib.size = lib.size, dispersion = d),
                design = design
            )
        },
        lower = 1e-4,
        upper = 1.0
    )$minimum
}

#' Variance Stabilizing Transformation
#'
#' Perform a Variance Stabilizing Transformation (VST) of a matrix of count data.
#' @param x A matrix of counts. Rows are features, columns are samples.
#' @param method VST method; see details.
#' @param lib.size Optional library sizes.
#' @param cpm Return log2 CPM instead of simple log2.
#' @param dispersion Optional dispersion for negative binomial VST.
#' @param design Optional design matrix for dispersion estimation.
#' @return Transformed matrix with attributes for method and library sizes.
#' @export
vst <- function(x, method = "anscombe.nb", lib.size = NULL, cpm = FALSE, dispersion = NULL, design = NULL) {
    x <- as.matrix(x)
    method.info <- vst_methods[[method]]
    is.null(method.info) && stop("Unknown method")

    if (nrow(x) == 0 || ncol(x) == 0) return(x)

    true.lib.size <- colSums(x)
    if (is.null(lib.size)) {
        good <- true.lib.size > 0
        norm.factors <- rep(1, length(good))
        norm.factors[good] <- edgeR::calcNormFactors(x[, good, drop = FALSE])
        lib.size <- true.lib.size * norm.factors
        lib.size.method <- "TMM normalization"
    } else if (all(lib.size == true.lib.size)) {
        lib.size.method <- "no adjustment"
    } else {
        lib.size.method <- "unknown"
    }
    lib.size <- pmax(lib.size, 1)
    mean.size <- mean(lib.size)
    x.norm <- t(t(x) * (mean.size / lib.size))

    if (method.info$needs.dispersion) {
        if (is.null(dispersion)) {
            dispersion <- .optimal_dispersion(x, method = method, lib.size = lib.size, design = design)
            cat("Dispersion estimated as ", dispersion, "\n", sep = "")
        }
        result <- method.info$vst(x.norm, dispersion)
    } else {
        result <- method.info$vst(x.norm)
    }

    if (cpm) {
        method.info$is.logish || stop("Counts Per Million not meaningful with this transform")
        result <- result + log2(1e6 / mean.size)
    }

    if (!is.null(colnames(x))) colnames(result) <- colnames(x)
    if (!is.null(rownames(x))) rownames(result) <- rownames(x)

    attr(result, "lib.size") <- lib.size
    attr(result, "true.lib.size") <- true.lib.size
    attr(result, "lib.size.method") <- lib.size.method
    attr(result, "cpm") <- cpm
    attr(result, "method") <- method
    if (method.info$needs.dispersion) attr(result, "dispersion") <- dispersion
    result
}

#' Advise how VST will transform data
#' @param what Either output of vst() or method name
#' @param dispersion As per vst()
#' @param cpm As per vst()
#' @param lib.size As per vst()
#' @return Data frame showing transformed counts and twofold steps
#' @export
vst_advice <- function(what = "anscombe.nb", dispersion = NULL, cpm = FALSE, lib.size = NULL) {
    if (!is.character(what)) {
        (is.null(dispersion) && is.null(lib.size)) || stop("Extra parameters not needed")
        method <- attr(what, "method")
        dispersion <- attr(what, "dispersion")
        cpm <- attr(what, "cpm")
        lib.size <- mean(attr(what, "lib.size"))
    } else {
        method <- what
    }

    method.info <- vst_methods[[method]]
    is.null(method.info) && stop("Unknown method")
    method.info$needs.dispersion && is.null(dispersion) && stop("dispersion needed")
    cpm && is.null(lib.size) && stop("lib.size needed")

    count <- cbind(c(0, 2^(0:12)))
    y <- vst(count, method = method, dispersion = dispersion, cpm = cpm, lib.size = lib.size)
    step <- rep(NA, nrow(count))
    step[3:nrow(count)] <- y[3:nrow(count), 1] - y[2:(nrow(count)-1), 1]

    data.frame(
        count = count[, 1],
        transformed_count = y[, 1],
        twofold_step = step
    )
}

