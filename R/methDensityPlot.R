#' methDensityPlot
#'
#' Plots the distribution of methylatio ratios in multiple samples using ggplot
#'
#' @param x A \code{GRanges}, \code{matrix} or \code{data.frame} of methylation ratios
#' @param samples \code{data.frame} describing the samples to plot
#' @param title Title for the plot
#' @return Called for the side effect of plotting
#'
#' @export
#'
#' @importFrom IRanges values
#' @importFrom ggplot2 ggplot aes ggtitle xlab ylab geom_line xlim
#' 
#' @author Aaron Statham <a.statham@@garvan.org.au>
methDensityPlot <- function(x, samples, title) {
	y <- Group <- Sample <- NULL #to shut up R CMD check
	if (class(x)=="GRanges") x <- as.matrix(values(x))
    stopifnot(ncol(x)==nrow(samples))
    dists <- do.call(rbind, lapply(1:ncol(x), function(i) {
        x <- density(x[,i], na.rm=TRUE)
        data.frame(Sample=samples$Sample[i], Group=samples$Group[i], x=x$x, y=x$y, stringsAsFactors=FALSE)
    }))

    ggplot(dists, aes(x=x, y=y, color=Group, group=Sample)) +
        ggtitle(title) + 
        xlab("Methylation ratio") + ylab("Relative Frequency") +
        geom_line() + xlim(-0.05, 1.05)
}
