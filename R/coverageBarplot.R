#' tabulateCoverage
#'
#' Tabulate how much of the genome is each annotation
#'
#' @param regions \code{GRangesList} of query regions, for example DMRs
#' @param background \code{GRanges} of interrogated genomic background
#' @param annotations \code{GRangesList} of subject regions, for example ChromHMM annotations
#' @return A \code{list} containing the elements:
#'  \item{coverage}{Bases intersected between ranges of \code{regions} and \code{annotations}}
#'  \item{ratios}{Percentage of genome covered by each intersection}
#'  \item{observed/expected}{Enrichment ratio of \code{regions} regions taking into account the \code{background}}
#'
#' @export
#'
#' @importFrom GenomicRanges GRangesList width intersect
#' @importFrom parallel mclapply
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
tabulateCoverage <- function(regions, background, annotations) {
    regions <- c(GRangesList("background"=background), regions)
    temp <- mclapply(1:length(regions), function(x) sapply(annotations, function(y) {
        sum(as.numeric(width(intersect(regions[[x]],y))))
    }))
    names(temp) <- names(regions)
    temp <- do.call(cbind, temp)
    temp <- cbind(temp, "genome"=sapply(annotations, function(x) sum(as.numeric(width(x)))))
    temp <- rbind(temp, "total"=c(sapply(regions, function(x) sum(as.numeric(width(x)))), NA))
    tempRat <- t(t(temp)/temp[nrow(temp),])*100
    tempRat2 <- tempRat[,-1]/tempRat[,"background"]
    list("coverage"=temp, "ratios"=tempRat, "observed/expected"=tempRat2[-nrow(tempRat2),-ncol(tempRat2), drop=FALSE])
}

#' coverageBarplot
#'
#' Create a barplot of enrichment/depletion of regions vs annotations, assess significance by permuting regions
#'
#' @param regions \code{GRangesList} of query regions, for example DMRs
#' @param background \code{GRanges} of interrogated genomic background
#' @param annotations \code{GRangesList} of subject regions, for example ChromHMM annotations
#' @param main Title of the plot
#' @param nperm Number of permutations to perform
#' @param cols Colours to use for the barplot
#' @param colBy Whether to colour by region or annotation
#' @param verbose Whether to print progress
#' @return Creates a barplot of observed/expected ratios of regions vs annotations with an asterisk to mark those which are significant, and returns a \code{list} of the plotted data containing the elements:
#'  \item{covTable}{Overlap of \code{regions} vs \code{annotations} as calculated by \code{tabulateCoverage}}
#'  \item{permTable}{Overlaps of \code{perms} vs \code{annotations}}
#'  \item{permQuant}{2.5\% and 97.5\% quantiles of observed/expected values obtained from permuted regions}
#'  \item{perms}{Permuted regions used to calculate p-values}
#'  \item{permPval}{P-values for enrichment and depletion}
#'
#' @export
#'
#' @importFrom GenomicRanges GRanges 
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
coverageBarplot <- function(regions, background, annotations, main, nperm=1000, cols=NULL, colBy=c("regions", "annotations"), verbose=FALSE) {
    colBy <- match.arg(colBy)
    if (colBy=="regions") 
        if (is.null(cols)) cols=rainbow(length(regions)) else stopifnot(length(cols)==length(regions))
    else
        if (is.null(cols)) cols=rainbow(length(annotations)) else stopifnot(length(cols)==length(annotations))
    #create sampled regions to assess significance
    nregions <- length(background)

    perms <- lapply(regions, function(x) {
        tmp <- sample(nregions, length(x)*nperm, replace=TRUE)
        tmp <- GRanges(background@seqnames[tmp], background@ranges[tmp])
        split(tmp, rep(1:nperm, each=length(x)))
    })
    if (verbose) message("Created permutations")

    covTable <- tabulateCoverage(regions, background, annotations)
    if (verbose) message("Tabulated coverage")

    permTable <- lapply(perms, tabulateCoverage, background, annotations)
    if (verbose) message("Tabulated coverage on permutations")
    permQuant <- lapply(permTable, function(x) apply(x$o, 1, quantile, c(0.025, 0.975)))
    permPval <- sapply(names(permTable), function(x)
        sapply(rownames(permTable[[x]]$o), function(y) min(
        (1-sum(covTable$o[y,x]>permTable[[x]]$o[y,])/nperm)*2,
        (1-sum(covTable$o[y,x]<permTable[[x]]$o[y,])/nperm)*2)))
    if (verbose) message("Calculated significance")

    #make barplot
    bp <- barplot(t(covTable$o), border=NA, beside=TRUE, col=cols, las=2, ylab="Observed/Expected", main=main, ylim=c(0, ceiling(max(covTable$o))))
    if (length(regions)==1) bp <- matrix(bp, nrow=1)
    for (i in 1:length(regions)) points(bp[i,], covTable$o[,i]+0.1, pch=ifelse(covTable$o[,i]<permQuant[[i]][1,] | covTable$o[,i]>permQuant[[i]][2,], "*", ""), cex=3)

    abline(h=c(0,1), col=c("black", "red"), lwd=c(1,3))
    if (colBy=="regions") {
        legend("topleft", fill=cols, legend=names(regions))
    }
    #return all the data
    list("covTable"=covTable, "permTable"=permTable, "permQuant"=permQuant, "perms"=perms, "permPval"=permPval)
}

