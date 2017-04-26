#' methylationPlotRegions
#'
#' Plot average methylation or occupancy ratio profiles around multiple sets of regions for a single sample
#'
#' @param meth \code{GRanges} of methylation count data (either CpG or GpC)
#' @param regions \code{GRangesList} to plot methylation/occupancy around, one line per element
#' @param samples \code{data.frame} describing the sample to calculate ratios for - must be a single row
#' @param GpC Whether the data is GpC data, and therefore plotting should be inverted
#' @param main Title for the plot
#' @param up Distance upstream of regions to plot
#' @param down Distance downstream of regions to plot
#' @param every Resolution to smooth at - distance between the start of adjacent smoothing windows
#' @param width Resolution to smooth at - width of each smoothing window
#' @param minCov Minimum sequencing coverage for a site to be used in calculations
#' @param addN Whether to add the number of regions to the legend
#' @return \code{GRangesList} of relative position of methylation sites to \code{regions} and corresponding methylation ratio 
#'
#' @export
#'
#' @importFrom GenomicRanges GRangesList values GRanges resize 
#' @importFrom IRanges IRanges
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_line xlab ylab ggtitle xlim ylim
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
methylationPlotRegions <- function(meth, regions, samples, GpC, main="", up=1000, down=1000, every=10, width=50, minCov=5, addN=TRUE) {
    position <- value <- Regions <- NULL # FUCK OFF R CMD CHECK
    stopifnot(nrow(samples)==1)
    stopifnot(all(c("Sample", "C", "cov") %in% colnames(samples)))
    if (class(regions)=="GRanges") {
        regions <- GRangesList(list("regions"=regions))
        cols <- 1
    }
    if (addN) names(regions) <- paste0(names(regions), " n=", prettyNum(elementNROWS(regions), big.mark=","))
    meth <- meth[values(meth)[[samples$cov]]>=minCov]

    meth.regions <- lapply(regions, function(x) {
        tmp <- overlapRegionFlanks(x, meth, up, down)
        tmp$ratio <- (values(meth)[[samples$C]]/values(meth)[[samples$cov]])[tmp$index]
        if (GpC) tmp$ratio <- 1-tmp$ratio
        GRanges(rep("chr", nrow(tmp)), IRanges(tmp$position, width=1), ratio=tmp$ratio)
    })

    at <- seq(-up, down, by=every)
    smoothing <- resize(GRanges(rep("chr", length(at)), IRanges(at, width=1)), width, fix="center")
    tmp <- data.frame(position=at, sapply(meth.regions, function(x) overlapMeans(x, smoothing, x$ratio)), check.names=FALSE)
    tmp <- melt(tmp, id.vars="position")
    names(tmp)[2] <- "Regions"
    p <- ggplot(tmp, aes(x=position, y=value, group=Regions, colour=Regions)) + 
        geom_line() + 
        xlab("Distance from feature (bp)") + 
        ylab(if(GpC) "1-GpC methylation" else "CpG methylation") +
        xlim(-up, down) + ylim(0, 1) +
        ggtitle(main)
    print(p)
    invisible(meth.regions)
}

#' methylationPlotSamples
#'
#' Plot average methylation or occupancy ratio profiles around a set of regions for multiple samples
#'
#' @param meth \code{GRanges} of methylation count data (either CpG or GpC)
#' @param regions \code{GRanges} to plot methylation/occupancy around, one line per element
#' @param samples \code{data.frame} describing the samples to calculate ratios for
#' @param GpC Whether the data is GpC data, and therefore plotting should be inverted
#' @param main Title for the plot
#' @param up Distance upstream of regions to plot
#' @param down Distance downstream of regions to plot
#' @param every Resolution to smooth at - distance between the start of adjacent smoothing windows
#' @param width Resolution to smooth at - width of each smoothing window
#' @param minCov Minimum sequencing coverage for a site to be used in calculations
#' @param addN Whether to add the number of regions to the title
#' @return \code{GRangesList} of per-sample relative position of methylation sites to \code{regions} and corresponding methylation ratio 
#'
#' @export
#'
#' @importFrom GenomicRanges GRangesList values GRanges resize 
#' @importFrom IRanges IRanges
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_line xlab ylab ggtitle xlim ylim
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
methylationPlotSamples <- function(meth, regions, samples, GpC, main="", up=1000, down=1000, every=10, width=50, minCov=5, addN=TRUE) {
    position <- value <- Sample <- NULL # FUCK OFF R CMD CHECK
    stopifnot(all(c("Sample", "C", "cov") %in% colnames(samples)))
    stopifnot(class(regions)=="GRanges")

    if (addN) main=paste0(main, " n=", prettyNum(length(regions), big.mark=","))

    meth.regions <- lapply(1:nrow(samples), function(i) {
        this.meth <- meth[values(meth)[[samples$cov[i]]]>=minCov]
        tmp <- overlapRegionFlanks(regions, this.meth, up, down)
        tmp$ratio <- (values(this.meth)[[samples$C[i]]]/values(this.meth)[[samples$cov[i]]])[tmp$index]
        if (GpC) tmp$ratio <- 1-tmp$ratio
        GRanges(rep("chr", nrow(tmp)), IRanges(tmp$position, width=1), ratio=tmp$ratio)
    })
    names(meth.regions) <- samples$Sample

    at <- seq(-up, down, by=every)
    smoothing <- resize(GRanges(rep("chr", length(at)), IRanges(at, width=1)), width, fix="center")
    tmp <- data.frame(position=at, sapply(meth.regions, function(x) overlapMeans(x, smoothing, x$ratio)), check.names=FALSE)
    tmp <- melt(tmp, id.vars="position")
    names(tmp)[2] <- "Sample"
    p <- ggplot(tmp, aes(x=position, y=value, group=Sample, colour=Sample)) + 
        geom_line() + 
        xlab("Distance from feature (bp)") + 
        ylab(if(GpC) "1-GpC methylation" else "CpG methylation") +
        xlim(-up, down) + ylim(0, 1) +
        ggtitle(main)
    print(p)
    invisible(meth.regions)
}

#' methylationBiPlot
#'
#' Plot average methylation and occupancy ratio profiles around a set of regions for a single sample
#'
#' @param methGpC \code{GRanges} of GpC methylation count data
#' @param methCpG \code{GRanges} of CpG methylation count data
#' @param regions \code{GRanges} to plot methylation and occupancy around
#' @param samples \code{data.frame} describing the sample to calculate ratios for - must be a single row
#' @param main Title for the plot
#' @param up Distance upstream of regions to plot
#' @param down Distance downstream of regions to plot
#' @param every Resolution to smooth at - distance between the start of adjacent smoothing windows
#' @param width Resolution to smooth at - width of each smoothing window
#' @param minCov Minimum sequencing coverage for a site to be used in calculations
#' @param addN Whether to add the number of regions to the title
#' @return \code{GRangesList} of relative position of GpC and CpG methylation sites to \code{regions} and corresponding methylation ratio
#'
#' @export
#'
#' @importFrom GenomicRanges values GRanges resize 
#' @importFrom IRanges IRanges
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_line xlab ylab ggtitle xlim ylim
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
methylationBiPlot <- function(methGpC, methCpG, regions, samples, main="", up=1000, down=1000, every=10, width=50, minCov=5, addN=TRUE) {
    position <- value <- Data <- NULL # FUCK OFF R CMD CHECK
    stopifnot(nrow(samples)==1)
    stopifnot(all(c("Sample", "C", "cov") %in% colnames(samples)))
    if (!class(regions)=="GRanges") stop("regions must be a GRanges")
    meth <- list("Occupancy"=methGpC[values(methGpC)[[samples$cov]]>=minCov],
               "Methylation"=methCpG[values(methCpG)[[samples$cov]]>=minCov])

    meth.regions <- lapply(meth, function(x) {
        tmp <- overlapRegionFlanks(regions, x, up, down)
        tmp$ratio <- (values(x)[[samples$C]]/values(x)[[samples$cov]])[tmp$index]
        GRanges(rep("chr", nrow(tmp)), IRanges(tmp$position, width=1), ratio=tmp$ratio)
    })
    meth.regions$Occupancy$ratio <- 1-meth.regions$Occupancy$ratio

    if (addN) main=paste0(main, " n=", prettyNum(length(regions), big.mark=","))

    at <- seq(-up, down, by=every)
    smoothing <- resize(GRanges(rep("chr", length(at)), IRanges(at, width=1)), width, fix="center")
    tmp <- data.frame(position=at, sapply(meth.regions, function(x) overlapMeans(x, smoothing, x$ratio)), check.names=FALSE)
    tmp <- melt(tmp, id.vars="position")
    names(tmp)[2] <- "Data"
    p <- ggplot(tmp, aes(x=position, y=value, group=Data, colour=Data)) + 
        geom_line() + 
        xlab("Distance from feature (bp)") + 
        ylab("1-GpC methylation/CpG methylation") +
        xlim(-up, down) + ylim(0, 1) +
        ggtitle(main)
    print(p)
    invisible(meth.regions)
}
