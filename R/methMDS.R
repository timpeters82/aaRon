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
#' @param cols Named vector of colours to use for each group
#' @return Called for the side effect of plotting
#'
#' @export
#'
#' @importFrom GenomicRanges values
#' @importFrom IRanges as.matrix
#' @importFrom Repitools genomeBlocks
#' @importFrom matrixStats rowVars
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
methMDS <- function(x, samples, smoothing=1, top=c(1, 10, 100), minCov=5, mc.cores=1, cols=NULL) {
    
    # mdsPlot directly stolen from minfi_1.12.0 as you cannot specify the colours
    mdsPlot <- function (dat, numPositions = 1000, sampNames = NULL, sampGroups = NULL, 
    xlim, ylim, pch = 1, cols=NULL, legendPos = "bottomleft", 
    legendNCol, main = NULL) {

        stopifnot(is(dat, "matrix"))
        b <- dat
        if (is.null(main)) 
            main <- sprintf("MDS\n%d most variable positions", 
                    numPositions)
        o <- order(-rowVars(b))[1:numPositions]
        d <- dist(t(b[o, ]))
        fit <- cmdscale(d)
        if (missing(xlim)) 
            xlim <- range(fit[, 1]) * 1.2
        if (missing(ylim)) 
            ylim <- range(fit[, 2]) * 1.2
        if (is.null(sampGroups)) 
            sampGroups <- rep(1, numPositions)

        groupNames <- unique(sampGroups)
        numGroups <- length(groupNames)
        if (is.null(cols)) {
            cols <- topo.colors(numGroups)
            names(cols) <- groupNames
        } else {
            stopifnot(all(names(cols) %in% samples$Group))
            stopifnot(all(samples$Group %in% names(cols)))
        }

        if (is.null(sampNames)) {
            plot(fit[, 1], fit[, 2], col = col, pch = pch, xlim = xlim, 
                ylim = ylim, xlab = "", ylab = "", main = main)
        } else {
            plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = "", 
                ylab = "", main = main)
            text(fit[, 1], fit[, 2], sampNames, col = cols[sampGroups])
        }
        
        if (missing(legendNCol)) 
            legendNCol <- numGroups
        if (numGroups > 1) {
            legend(legendPos, legend = names(cols), ncol = legendNCol, 
                text.col = cols)
        }
    }

    if (!is.null(cols)) {
        stopifnot(all(names(cols) %in% samples$Group))
        stopifnot(all(samples$Group %in% names(cols)))
    }
    smoothing <- trunc(smoothing)
    stopifnot(smoothing>0)
    stopifnot (all(top>0) && all(top<=100))
    stopifnot(all(c("Sample", "C", "cov") %in% colnames(samples)))
    if (smoothing==1) ratios <- as.matrix(values(methRatios(x, samples, minCov)))
    else {
        stopifnot(all(!is.na(seqlengths(x))))
        tmp <- genomeBlocks(seqlengths(x), seqlevels(x), smoothing)
        tmp <- tmp[tmp %over% x]
        ratios <- as.matrix(values(methWindowRatios(x, tmp, samples, minCov, mc.cores)))
    }
    ratios <- ratios[rowMeans(is.na(ratios))<1, , drop=FALSE]
    if (nrow(ratios)<1000) stop("Fewer than 1000 sites have data above the minimum coverage, try inclreasing smoothing")    
    for (i in top) {
        npos <- trunc(nrow(ratios)*i/100)
        main <- paste0("MDS - top ", prettyNum(npos, big.mark=","), " sites")
        if (smoothing!=1) main <- paste0(main, " (", round(smoothing/1000,1), "kb smoothed)")
        mdsPlot(ratios, numPositions=npos, sampNames=samples$Sample, sampGroups=samples$Group, main=main,
            cols=cols)
    }
}
