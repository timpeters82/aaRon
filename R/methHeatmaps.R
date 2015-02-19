#' methHeatmap
#'
#' Plots a region of the genome with per sample heatmaps of methylation
#'
#' @param x CpG-site level methylation data, including the \code{w} column of weights if \code{weighted} is true
#' @param region The genomic region to be plotted
#' @param samples \code{data.frame} describing the samples to plot
#' @param CpGislands \code{GRanges} of CpG island positions
#' @param tx \code{GRanges} of transcripts
#' @param extras \code{List} of extra \code{Gviz} tracks to plot above the heatmaps
#' @param smoothing If non-zero, the number of basepairs to summarise methyation over. If zero no smoothing is performed.
#' @param cols Named vector of colours to use for each group
#' @param minCov Minimum coverage for a CpG site to be included in the weighted mean
#' @param ... Additional parameters passed to \code{Gviz::plotTracks}
#' @return Create a heatmap of methylation
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom GenomicRanges subsetByOverlaps GRanges
#' @importFrom IRanges IRanges
#' @importFrom Gviz DataTrack GenomeAxisTrack AnnotationTrack plotTracks chromosome<-
#' @importFrom S4Vectors endoapply
#' @importFrom GenomeInfoDb seqnames
#'
#' @export
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
methHeatmap <- function(x, region, samples, CpGislands, tx, extras=list(), smoothing=0, cols=NULL, minCov=2, ...) {
    stopifnot(length(region)==1)

    groups <- unique(samples$Group)
    if (is.null(cols)) {
        cols <- brewer.pal(length(groups), "Dark2")[1:length(groups)]
        names(cols) <- groups
    } else {
        stopifnot(all(names(cols) %in% samples$Group))
        stopifnot(all(samples$Group %in% names(cols)))
    }

    # reduce the data we are working with
    x <- subsetByOverlaps(x, region)

    # calculate CpG ratio inside each bin PER GROUP
    if (smoothing==0) {
        dt.group <- lapply(groups, function(i)
            DataTrack(methRatios(x, samples[samples$Group==i,], minCov), name=i, background.title=cols[i],
                type="heatmap", showSampleNames=TRUE, ylim=c(0, 1),
                gradient=c("blue", "white", "red")))

        # now add type="b" plot of group means
        dt.group <- c(dt.group, list(DataTrack(methRatios(x, samples, minCov), groups=samples$Group, type="b",
            aggregateGroups=TRUE, col=cols, ylim=c(0, 1), name="Group means")))
    } else {
        bins <- GRanges(seqnames(region),
            IRanges(seq.int(start(region), end(region)-smoothing, smoothing), width=smoothing))
        dt.group <- lapply(groups, function(i)
            DataTrack(methWindowRatios(x, bins, samples[samples$Group==i,], minCov), name=i, background.title=cols[i],
                type="heatmap", showSampleNames=TRUE, ylim=c(0, 1),
                gradient=c("blue", "cyan", "green", "yellow", "red")))

        # now add type="b" plot of group means
        dt.group <- c(dt.group, list(DataTrack(methWindowRatios(x, bins, samples, minCov), groups=samples$Group, type="b",
            aggregateGroups=TRUE, col=cols, ylim=c(0, 1), name="Group means")))

    } 
    tx$group <- tx$id <- tx$gene_name
    extras <- endoapply(extras, function(x) {
        chromosome(x) <- as.character(seqnames(region))
        x
    })
    basetracks <- list(GenomeAxisTrack(),
        AnnotationTrack(subsetByOverlaps(tx, region), name="tx", showId=TRUE, fill="blue", col=NULL),
        AnnotationTrack(subsetByOverlaps(CpGislands, region), name="CpGi", fill="green", col=NULL),
        AnnotationTrack(unvalue(x), name="CpGs", fill="green", col=NULL, stacking="dense"))
    plotTracks(c(basetracks, extras, dt.group), from=start(region), to=end(region), ...)
}

