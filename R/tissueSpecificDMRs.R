#' tissueSpecificDMRs
#'
#' Find regions of tissue specific methylation differences using a sliding chi-squared test. Method adapted from "Adult tissue methylomes harbor epigenetic memory at embryonic enhancers", Hon et al, Nature Genetics 2013
#'
#' @param x CpG-site level methylation data
#' @param samples \code{data.frame} describing the samples to process
#' @param k Number of adjacent CpG sites to calculate chi-sq statistics upon
#' @param maxCov Maximum sequencing depth allowable at a single CpG site; essential for removing repeat regions
#' @param p.cutoff Chi-sq BY adjusted p-value for a CpG site to be considered significant
#' @param minCpGs Minimum number of adjacent CpG sites below \code{p.cutoff} for a region to be considered significant
#' @param minWidth Minimum size in basepairs for a tsDMR
#' @param minCpGdensity Minimum number of CpG sites per kilobase in a tsDMR
#' @param minSd Minimum standard deviation of methylaion ratios of a tsDMR across the tissues
#' @param mc.cores Number of cores to use
#' @return \code{GRanges} of significant regions, with mean-methyaltion per sample and chisq value in the values()
#'
#' @importFrom zoo rollsum
#' @importFrom GenomicRanges values GRangesList width seqnames GRanges countOverlaps
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom S4Vectors Rle runValue start end
#' @importFrom IRanges IRanges
#' @importFrom matrixStats rowSds
#'
#' @export
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
tissueSpecificDMRs <- function(x, samples, k=3, maxCov=200, p.cutoff=0.01, minCpGs=5, minWidth=100,
    minCpGdensity=5, minSd=0.1, mc.cores=1) {
    message("Masking out NAs and CpG sites of too high coverage")
    for (i in 1:nrow(samples)) {
        NAs <- is.na(values(x)[[samples$C[i]]]) | is.na(values(x)[[samples$cov[i]]])
        overCov <- values(x)[[samples$cov[i]]]>maxCov
        values(x)[[samples$C[i]]][NAs | overCov] <- values(x)[[samples$cov[i]]][NAs | overCov] <- 0
    }

    message(paste0("Performing sliding window sums across ", k," adjacent CpG sites using ", mc.cores, " cores"))
    unit.cov <- simplify2array(mclapply(samples$cov, function(i) rollsum(values(x)[[i]], k, na.pad=TRUE),
        mc.cores=mc.cores))
    message (" - completed coverage")
    unit.C <- simplify2array(mclapply(samples$C, function(i) rollsum(values(x)[[i]], k, na.pad=TRUE),
        mc.cores=mc.cores))
    message (" - completed methylation")
    unit.meth <- rowSums(unit.C)/rowSums(unit.cov)

    message("Calculating chi-squared statistic")
    unit.expected <- unit.cov*unit.meth
    unit.chisq <- rowSums(((unit.C-unit.expected)^2)/unit.expected)

    unit.chisq[is.na(unit.chisq)] <- 0
    unit.chisq[rowMeans(as.matrix(values(x))[, samples$cov])>maxCov] <- 0

    x$chisq <- unit.chisq
    rm(k, unit.cov, unit.C, unit.meth, unit.expected, unit.chisq)

    message("Calculating chi-squared p-value")
    x$pval <- pchisq(x$chisq, nrow(samples)-1, lower.tail=FALSE)
    x$adj.pval <- p.adjust(x$pval, "BY")

    message(paste0("Thresholding at chisq adj.pval<", p.cutoff))
    tsDMR <- unlist(GRangesList(lapply(seqlevels(x), function(chr) {
        message(paste("Processing", chr))
        tmp <- x[as.character(seqnames(x))==chr]
        signif <- Rle(tmp$adj.pval<p.cutoff)
        SOI <- runValue(signif) & width(signif)>=minCpGs
        out <- GRanges(rep(chr, sum(SOI)), IRanges(start(tmp)[start(signif)[SOI]], end(tmp)[end(signif)[SOI]]))
        # Filter for minimum width
        out <- out[width(out)>=minWidth]
        # Filter for minimum CpG density
        out <- out[countOverlaps(out, tmp)/width(out)*1000>=minCpGdensity]
        # Calculate methlation ratios
        out <- methWindowRatios(tmp, out, samples, mc.cores=mc.cores)
        # Filter for minimum standard deviation
        out <- out[which(rowSds(as.matrix(values(out)))>=minSd)]
        out$chisq <- overlapMeans(tmp, out, tmp$chisq)    
        message(paste("Completed", chr))
        out
    })))
    tsDMR
}
