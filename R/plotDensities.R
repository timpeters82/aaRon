plotDensities <- function(s, col=NULL, xlim=NULL, ylim=NULL, na.rm=TRUE, main="", xlab="", ylab="Relative Frequency", ...) {
    if (is.matrix(s)) s <- data.frame(s)
    if (!is.list(s)) s <- list(s)
    if (is.null(col)) col <- rainbow(length(s))
    sD <- lapply(s, density, na.rm=na.rm)
    if (is.null(xlim)) xlim <- range(sapply(sD, function(x) range(x$x)))
    if (is.null(ylim)) ylim <- range(sapply(sD, function(x) range(x$y)))
    plot(0, type="n", xlim=xlim, ylim=ylim, main=main, xlab=xlab, ylab=ylab, ...)
    for (i in 1:length(sD)) lines(sD[[i]], col=col[i], ...)
}

