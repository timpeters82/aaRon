#' grangesVenn
#'
#' Creates a table of overlapping \code{GRanges} from a \code{GRangesList}, vennCounts and optionally plots the first (up to) five
#'
#' @param x \code{GRangesList} of interest
#' @param plot Whether to plot after calculating the 
#' @param ... Additional parameters passed on to vennDiagram
#' @return A list containing the pooled and reduced ranges, a boolean matrix of whether each sample in \code{x} overlaps each \code{Ranges} and the vennCounts
#'
#' @export
#'
#' @importFrom GenomicRanges reduce unlist
#' @importFrom limma vennCounts vennDiagram
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
grangesVenn <- function(x, plot=TRUE, ...) {
	x.pool <- reduce(unlist(x))
	x.ov <- sapply(x, function(y) x.pool %over% y)
	x.counts <- vennCounts(x.ov)
	if (plot) vennDiagram(x.ov[,1:min(5, ncol(x.ov))], ...)
	invisible(list("Ranges"=x.pool, "Overlaps"=x.ov, "Counts"=x.counts))
}
