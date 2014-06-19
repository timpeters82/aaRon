#' findDMVs
#'
#' Finds DNA methylation valleys using the procedure described in Hovestaft 2014
#'
#' @param x \code{GRanges} of methylation count data
#' @param samples \code{data.frame} describing the samples to discover DMVs for
#' @param minSize Minimum number of basepairs with methylation below \code{cutoff} for a region to be called as a DMV
#' @param wsize Window size to take weighted means within across the genome
#' @param step Distance to step windows across the genome by
#' @param cutoff Maximum weighted mean of a window to contribute to a DMV
#' @param minCov Minumum sequencing coverage for a CpG site to contribute to the weighted mean calculation
#' @return \code{GRangesList} of DMVs found in each sample
#' 
#' importFrom GenomicRanges values seqlengths seqnames GRangesList reduce width split
#' importFrom Repitools genomeBlocks
#' importFrom parallel multicore
#' 
#' @export
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
findDMVs <- function(x, samples, minSize=5000, wsize=1000, step=100, cutoff=0.15, minCov=5) {
	stopifnot("w" %in% names(values(x)))
	stopifnot(all(!is.na(seqlengths(x))))

	bins <- genomeBlocks(seqlengths(x), width=wsize, spacing=step)

    x.rat <- methRatios(x, samples, minCov)

    # split through chromosomes and mclapply
    bins <- split(bins, seqnames(bins))
    x.rat <- split(x.rat, seqnames(x.rat))
    w <- split(x$w, as.character(seqnames(x)))
    wm <- unlist(GRangesList(mclapply(names(bins), function(i) {
        message("Processing ", i)
        tmp <- bins[[i]]
        values(tmp) <- weightedMethylation(x.rat[[i]], bins[[i]], w[[i]])
        message(i, " Finished!")
        tmp
    })))
    rm(bins, x.rat, w)
    bins <- unvalue(wm)
    DMVs <- GRangesList(mclapply(1:ncol(values(wm)), function(i) {
        tmp <- reduce(bins[which(values(wm)[[i]]<cutoff)])
        tmp[width(tmp)>=minSize]
    }))
    names(DMVs) <- samples$Sample
    DMVs
}
