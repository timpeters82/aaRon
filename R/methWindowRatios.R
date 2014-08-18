#' methWindowRatios
#'
#' Calculates methylation ratios from count data per each supplied window
#'
#' @param x \code{GRanges} of methylation count data
#' @param windows \code{GRanges} of windows to calculate methylation ratios across
#' @param samples \code{data.frame} describing the samples to calculate ratios for
#' @param minCov minimum amount of coverage to calculate a ratio
#' @param mc.cores Number of cores to use
#' @return \code{GRanges} of methylation ratios
#'
#' @export
#'
#' @importFrom parallel mclapply
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
methWindowRatios <- function(x, windows, samples, minCov=5, mc.cores=1) {
    tmp <- simplify2array(mclapply(1:nrow(samples),
        function(i) overlapRatios(x, windows, samples$C[i], samples$cov[i]),
        mc.cores = mc.cores))
    if (length(windows)==1) tmp <- matrix(tmp, nrow=1)
    values(windows) <- tmp
    names(values(windows)) <- samples$Sample
    windows
}
