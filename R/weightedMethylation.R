#' cpgWeights
#'
#' Calculate weights for CpG density weighted means, as used in Berman BP 2011
#'
#' @param x A \code{GRanges} of CpG site positions to calculate weighting factors for
#' @return Integer vector of weights in the same order as the supplied \code{x}
#'
#' @export
#'
#' @importFrom GenomicRanges start seqnames
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
cpgWeights <- function(x) {
	n <- length(x)
	# Check that CpG sites within each chromosome are sorted
	stopifnot(all(unlist(sapply(split(start(x), as.character(seqnames(x))), diff))>0))
	x.d <- diff(start(x))
	x.d[x.d<=0] <- 0L
	w <- rep(0L, n)
	w[-c(1, n)] <- x.d[-(n-1)]+x.d[-1]
	w
}

#' weightedMethylation
#' 
#' Calculate a per sample weighted mean of CpG-site level methylation data
#'
#' @param x \code{GRanges} of CpG-site level methylation ratios
#' @param regions \code{GRanges} of genomic regions to be plotted
#' @param w \code{vector} of CpG-site level weights (optional)
#' @return \code{matrix} of weighted methylation values
#'
#' @export
#' 
#' @importFrom matrixStats colWeightedMeans
#' @importFrom GenomicRanges findOverlaps values
weightedMethylation <- function(x, regions, w) {
	if (missing(w)) w <- rep(1, length(x))
	stopifnot(length(x)==length(w))
	x.r <- as.matrix(values(x))
	ov <- as.matrix(findOverlaps(regions, x))
	res <- matrix(as.numeric(NA), nrow=length(regions), ncol=ncol(x.r))
	res[unique(ov[,1]),] <- t(sapply(split(ov[,2], ov[,1]), function(y) colWeightedMeans(x.r[y,,drop=FALSE], w[y], na.rm=TRUE)))
	res
}

