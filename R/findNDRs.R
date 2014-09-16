#' findNDRs
#'
#' Find NDRs in NOMe-seq data using a sliding window chi-squared test versus whole genome background
#'
#' @param x \code{GRanges} of methylation count data at GCH sites
#' @param samples \code{data.frame} describing the samples to discover NDRs for
#' @param p.cutoff Minimum -log10 p.value for a window to be considered significant
#' @param windowWidth Size of each windoe
#' @param windowBy Distance between the starts of consecutive windows
#' @param minSize Minimum size of an NDR when overlapping significant windows are pooled
#' @return \code{GRangesList} of NDRs found in each sample with the mean pvalue for each region
#' 
#' importFrom compiler cmpfun
#' importFrom Repitools genomeBlocks
#' importFrom GenomicRanges seqlengths findOverlaps GRangesList values reduce width
#'
#' @export
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
findNDRs <- function(x, samples, p.cutoff=15, windowWidth=100, windowBy=20, minSize=140) {
	fast.chisq <- compiler::cmpfun(function(x, p) {
	    n <- rowSums(x)
	    E <- cbind(n * p[1], n * p[2])
	    STATISTIC <- rowSums((x - E)^2/E)
	    PARAMETER <- 1
	    PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
	    return(PVAL)
	})

	# x must have seqlengths defined
	stopifnot(all(!is.na(seqlengths(x))))

	# Define windows for pvalue calculation
	windows <- genomeBlocks(seqlengths(x), width=windowWidth, spacing=windowBy)

	# Only keep windows that overlap x
	windows <- windows[windows %over% x]

	# Find NDRs
	NDRs <- GRangesList(lapply(1:nrow(samples), function(i) {
		message(paste0("Processing ", samples$Sample[i]))
		# Counts of Cs and Ts in each window
		message(" - Counting Cs")
		windows.counts <- cbind("C"=overlapSums(x, windows, values(x)[[samples$C[i]]]))
		message(" - Counting Ts")
		windows.counts <- cbind(windows.counts, "T"=overlapSums(x, windows, values(x)[[samples$cov[i]]])-windows.counts[,"C"])

	    # Genomic background ratio to compare to
		Csum <- sum(values(x)[[samples$C[i]]])
		Tsum <- sum(values(x)[[samples$cov[i]]])-sum(values(x)[[samples$C[i]]])
	    C.T <- c(Csum,Tsum)/(Csum+Tsum)

	    # calculate chisq pvalue
		message(" - Calculating chisq pvalues")
	    windows.pvals <- -log10(fast.chisq(windows.counts, C.T))
	    windows.pvals[windows.pvals==Inf] <- max(windows.pvals[!windows.pvals==Inf])

	    # find regions that meet the cutoff on pvalue and size
		message(" - Finding significant regions")
        windows.cutoff <- reduce(windows[which(windows.pvals>p.cutoff)])
        windows.cutoff$p.mean <- overlapMeans(windows[!is.na(windows.pvals)], windows.cutoff, windows.pvals[!is.na(windows.pvals)])
#        windows.cutoff$p.mean <- overlapMeans(windows, windows.cutoff, windows.pvals, na.rm=TRUE)
        windows.cutoff <- windows.cutoff[width(windows.cutoff)>=minSize]
	}))

	names(NDRs) <- samples$Sample
	seqlevels(NDRs) <- seqlevels(x)
	seqlengths(NDRs) <- seqlengths(x)
	NDRs
}
