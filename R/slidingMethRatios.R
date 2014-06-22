#' slidingMethRatios
#'
#' Calculates methylation ratios means across a sliding window of a number of CpG sites
#'
#' @param x \code{GRanges} of methylation count data
#' @param samples \code{data.frame} describing the samples to calculate ratios for
#' @param n Number of adjacent CpG sites to average
#' @param minCov minimum amount of coverage to calculate a ratio
#' @return \code{GRanges} of methylation ratio means across \code{n} CpG sites
#'
#' @export
#'
#' @importFrom GenomicRanges split seqnames values values<- GRangesList
#' @importFrom zoo rollapply
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
slidingMethRatios <- function(x, samples, n=5, minCov=5) {
    x <- methRatios(x, samples, minCov)
    x.split <- GenomicRanges::split(x, seqnames(x))
    x.result <- unlist(GRangesList(mclapply(names(x.split), function(i) {
        message("Processing ", i)
        tmp <- x.split[[i]]
        values(tmp) <- rollapply(as.matrix(values(tmp)), n, mean, na.rm=TRUE, fill=NA)
        message(i, " finished!")
        tmp
    })))
    x.result
}
