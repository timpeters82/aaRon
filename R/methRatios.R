#' methRatios
#'
#' Calculates methylation ratios from count data
#'
#' @param x \code{GRanges} of methylation count data
#' @param samples \code{data.frame} describing the samples to calculate ratios for
#' @param minCov minimum amount of coverage to calculate a ratio
#' @return \code{GRanges} of methylation ratios
#'
#' @export
#'
#' @importFrom GenomicRanges values values<-
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
methRatios <- function(x, samples, minCov=5) {
    stopifnot(all(c("Sample", "C", "cov") %in% colnames(samples)))
	tmp <- as.matrix(values(x)[samples$cov])
	tmp[tmp<minCov] <- NA
	values(x) <- as.matrix(values(x)[samples$C])/tmp
	names(values(x)) <- samples$Sample
	x
}