#' methMDS
#'
#' Plots multidimensional scaling plot of methylation count data
#'
#' @param x \code{GRanges} of methylation count data
#' @param samples \code{data.frame} describing the samples to calculate ratios for
#' @param smoothing Number of basepairs to smooth methylation data over, \code{1} indicates no smoothing
#' @param top \code{numeric vector} of percentage of top varying methlation ratios to include in the MDS, defaults to the top 1\%, 10\% and 100\%
#' @param minCov minimum amount of coverage to calculate a ratio
#' @param mc.cores Number of cores to use
#' @return Called for the side effect of plotting
#'
#' @export
#'
#' @importFrom GenomicRanges values
#' @importFrom IRanges as.matrix
#' @importFrom Repitools genomeBlocks
#' @importFrom minfi mdsPlot
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
methMDS <- function(x, samples, smoothing=1, top=c(1, 10, 100), minCov=5, mc.cores=1) {
    smoothing <- trunc(smoothing)
    stopifnot(smoothing>0)
    stopifnot (all(top>0) && all(top<=100))
    if (smoothing==1) ratios <- as.matrix(values(methRatios(x, samples, minCov)))
    else {
        tmp <- genomeBlocks(seqlengths(x), seqlevels(x), smoothing)
        tmp <- tmp[tmp %over% x]
        ratios <- as.matrix(values(methWindowRatios(x, tmp, samples, mc.cores=mc.cores)))
    }
    ratios <- ratios[rowMeans(is.na(ratios))<1, , drop=FALSE]
    if (nrow(ratios)<1000) stop("Fewer than 1000 sites have data above the minimum coverage, try inclreasing smoothing")    
    for (i in top) {
        npos <- trunc(nrow(ratios)*i/100)
        main <- paste0("MDS - top ", prettyNum(npos, big.mark=","), " sites")
        if (smoothing!=1) main <- paste0(main, " (", round(smoothing/1000,1), "kb smoothed)")
        mdsPlot(ratios, numPositions=npos, sampNames=samples$Sample, sampGroups=samples$Group, main=main)
    }
}
