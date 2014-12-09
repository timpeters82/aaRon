#' grangesVenn
#'
#' Creates a table of overlapping \code{GRanges} from a \code{GRangesList}, vennCounts and optionally plots the first (up to) five
#'
#' @param x \code{GRangesList} of interest
#' @param plot Whether to plot after calculating the intersections
#' @param basepair Whether to calculate the intersections on a per-basepair level rather than per region
#' @param ... Additional parameters passed on to vennDiagram
#' @return A list containing the pooled and reduced ranges, a boolean matrix of whether each sample in \code{x} overlaps each \code{Ranges} and the vennCounts
#'
#' @export
#'
#' @importFrom GenomicRanges reduce unlist
#' @importFrom limma vennCounts vennDiagram
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
grangesVenn <- function(x, plot=TRUE, basepair=FALSE, ...) {
	stopifnot(inherits(x, "GRangesList"))
	x.pool <- reduce(unlist(x))
	if (basepair) {
		if (length(x)>5) {
			warning("Only using the first 5 elements of 'x'")
			x <- x[1:5]
		}
		# blatantly adapted from limma::vennCounts
		ncontrasts <- length(x)
		names <- names(x)
		if (is.null(names)) 
		    names <- paste("Group", 1:ncontrasts)
		noutcomes <- 2^ncontrasts
		outcomes <- matrix(0, noutcomes, ncontrasts)
		colnames(outcomes) <- names

		for (j in 1:ncontrasts) outcomes[, j] <- rep(0:1, times = 2^(j - 
		    1), each = 2^(ncontrasts - j))

		x.outcomes <- GRangesList(lapply(1:noutcomes, function(i) {
		    tmp <- x.pool
		    for (j in 1:ncontrasts)
		        if (outcomes[i,j]==1)
		            tmp <- intersect(tmp, x[[j]]) else tmp <- setdiff(tmp, x[[j]])
		    tmp
		}))

		x.counts <- structure(cbind(outcomes,
			Counts=paste0(prettyNum(round(sum(width(x.outcomes))/1e6,1), big.mark=","), "Mb")),
			class = "VennCounts")
		if (plot) vennDiagram(x.counts, ...)
		invisible(list("Ranges"=x.pool, "Counts"=x.counts))
	} else {
		x.ov <- sapply(x, function(y) x.pool %over% y)
		x.counts <- vennCounts(x.ov)
		if (plot) vennDiagram(x.ov[,1:min(5, ncol(x.ov))], ...)
		invisible(list("Ranges"=x.pool, "Overlaps"=x.ov, "Counts"=x.counts))
	}
}
