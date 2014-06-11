#' seqHeatmaps
#'
#' Plot a grid of heatmaps of ChIP-seq signals around regions of interest
#'
#' @param bamfiles Paths to BAM files with names
#' @param anno \code{GRangesList} or \code{GRanges} of regions to sample around
#' @param up How far upstream to sample
#' @param down How far downstream to sample
#' @param freq How many basepairs apart to sample
#' @param s.width How many basepairs to smooth the data by
#' @param nc Number of columns to plot per page
#' @param nr Number of rows to plot per page
#' @param cols colors to use in the heatmap
#' @param zlim limit to plot in the heatmap, defaults to 0-1 read per million
#' @param n.sample Number of regions from each anno to randomly sample for plotting, NULL disables the sampling
#' @return Plots a set of heatmaps and invisibly returns the plotted data
#'
#' @export
#' 
#' @importFrom Repitools featureScores tables
#' @importFrom GenomicRanges unlist
#' @importFrom IRanges elementLengths
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
seqHeatmaps <- function(bamfiles, anno, up=2000, down=2000, freq=20, s.width=200, nc=length(anno), nr=length(bamfiles), cols=colorRampPalette(c("white", "red"))(100), zlim=c(0, 1), n.sample=2000) {
	stopifnot(all(file.exists(bamfiles)))
	stopifnot(length(names(bamfiles))==length(bamfiles))
	if (class(anno)=="GRanges") {
		anno <- GRangesList("All"=anno)
	}
	stopifnot(length(names(anno))==length(anno))
	stopifnot(nc<=length(anno))
	stopifnot(nr<=length(bamfiles))
	stopifnot(all(elementLengths(anno)>=n.sample) | is.null(n.sample))

	# Get scores for all annos, scale up to rpm
	b.scores <- lapply(tables(featureScores(bamfiles, unlist(anno), up=up, down=down, freq=freq, s.width=s.width)), function(x) x*1e6)
	pos <- as.integer(colnames(b.scores[[1]]))
	
	# split per anno
	b.scores <- lapply(b.scores, function(x) lapply(split(as.data.frame(x), rep(names(anno), elementLengths(anno))), as.matrix))

	# sample (if enabled)
	if (is.null(n.sample)) b.sample <- b.scores else {
		b.sample <- lapply(elementLengths(anno), sample, n.sample)
		b.sample <- lapply(b.scores, function(x) lapply(1:length(anno), function(y) x[[y]][b.sample[[y]],]))
	}

	# plot
	par(mfrow=c(nr, nc))
	for (i in 1:length(bamfiles)) for (j in 1:length(anno))
		image(x=pos, z=t(b.sample[[i]][[j]]), zlim=zlim, col=cols, yaxt="n", xlab="Distance from annotation centre", main=if (i==1) names(anno)[j] else "", ylab=if (j==1) names(bamfiles)[i] else "")
    # Return the data
	invisible(list("Scores"=b.scores, "Sampled"=b.sample))
}
