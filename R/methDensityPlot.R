#' methDensityPlot
#'
#' Plots the distribution of methylatio ratios in multiple samples using ggplot
#'
#' @param x A \code{GRanges}, \code{matrix} or \code{data.frame} of methylation ratios
#' @param samples \code{data.frame} describing the samples to plot
#' @param title Title for the plot
#' @param cols Named vector of colours to use 
#' @return Called for the side effect of plotting
#'
#' @export
#'
#' @importFrom IRanges values
#' @importFrom ggplot2 ggplot aes ggtitle xlab ylab geom_line xlim scale_color_manual
#' 
#' @author Aaron Statham <a.statham@@garvan.org.au>
methDensityPlot <- function(x, samples, title, cols=NULL) {
	y <- Group <- Sample <- NULL #to shut up R CMD check
    if (is.null(samples$Group)) {
        samples$Group <- samples$Sample
    }
    if (!is.null(cols)) {
        stopifnot(all(names(cols) %in% samples$Group))
        stopifnot(all(samples$Group %in% names(cols)))
    }
	if (class(x)=="GRanges") x <- as.matrix(values(x))
    stopifnot(ncol(x)==nrow(samples))
    dists <- do.call(rbind, lapply(1:ncol(x), function(i) {
        x <- density(x[,i], na.rm=TRUE)
        data.frame(Sample=samples$Sample[i], Group=samples$Group[i], x=x$x, y=x$y, stringsAsFactors=FALSE)
    }))

    p <- ggplot(dists, aes(x=x, y=y, color=Group, group=Sample)) +
        ggtitle(title) + 
        xlab("Methylation ratio") + ylab("Relative Frequency") +
        geom_line() + xlim(-0.05, 1.05)
    if (!is.null(cols)) p <- p + scale_color_manual(values=cols)
    p
}
