#' biasPlots
#'
#' Plot of sequencing coverage bias versus GC and CpG content in 100bp genomic bins
#'
#' @param x \code{GRanges} of methylation count data - is also used to calculate CpG density so should not be pre-filtered
#' @param samples \code{data.frame} describing the samples to calculate ratios for
#' @param build BSgenome object of the organism the data is from
#' @return Called for the side effect of plotting
#'
#' @export
#'
#' @importFrom IRanges values
#' @importFrom Repitools genomeBlocks
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom Biostrings letterFrequency
#' @importFrom GenomicRanges countOverlaps subsetByOverlaps
#' @importFrom matrixStats colMedians
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes ggtitle xlab ylab geom_line expand_limits
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
biasPlots <- function(x, samples, build) {
    y <- Sample <- NULL #to shut up R CMD check
    stopifnot(all(c("Sample", "C", "cov") %in% colnames(samples)))
    # minimum number of data points to have the median calculated
    minBlocks <- 100

    # scaling factor to normalise sequencing depth between samples
    cov.norm <- colSums(as.matrix(values(x)[,samples$cov]))
    cov.norm <- max(cov.norm)/cov.norm

    bl <- genomeBlocks(build, unique(as.character(seqlevels(x))), 100)
    bl$GC <- letterFrequency(getSeq(build, bl), "GC")
    bl$CpG <- countOverlaps(bl, x)

    # remove GC/CpG==0 regions
    bl <- bl[bl$GC>0 & bl$CpG>0]

    # find and plot median coverage for each GC and CpG bin
    bl.GC.counts <- table(bl$GC)
    bl.GC.use <- as.numeric(names(bl.GC.counts)[bl.GC.counts>=minBlocks])
    bl.GC <- bl[bl$GC %in% bl.GC.use]
    bl.GC.medians <- sapply(split(bl.GC, bl.GC$GC), function(y) 
        colMedians(as.matrix(values(subsetByOverlaps(x, y))[,samples$cov])))*cov.norm
    rownames(bl.GC.medians) <- samples$Sample
    bl.GC.melt <- melt(bl.GC.medians)
    colnames(bl.GC.melt) <- c("Sample", "x", "y")

    p <- ggplot(bl.GC.melt, aes(x=x, y=y, color=Sample)) + ggtitle("Sequencing depth bias vs GC content") +
        xlab("GC content % of 100bp windows") + ylab("Sequening coverage (Normalised median)") + geom_line() +
        expand_limits(y=0)
    print(p)

    bl.CpG.counts <- table(bl$CpG)
    bl.CpG.use <- as.numeric(names(bl.CpG.counts)[bl.CpG.counts>=minBlocks])
    bl.CpG <- bl[bl$CpG %in% bl.CpG.use]
    bl.CpG.medians <- sapply(split(bl.CpG, bl.CpG$CpG), function(y) 
        colMedians(as.matrix(values(subsetByOverlaps(x, y))[,samples$cov])))*cov.norm
    rownames(bl.CpG.medians) <- samples$Sample
    bl.CpG.melt <- melt(bl.CpG.medians)
    colnames(bl.CpG.melt) <- c("Sample", "x", "y")

    p <- ggplot(bl.CpG.melt, aes(x=x, y=y, color=Sample)) + ggtitle("Sequencing depth bias vs CpG content") +
        xlab("CpG content % of 100bp windows") + ylab("Sequening coverage (Normalised median)") + geom_line() +
        expand_limits(y=0)
    print(p)
}
