#' overlapMeans
#'
#' Calculate the mean of values associated with ranges which overlap query ranges
#'
#' @param x The \code{GRanges} object to find overlaps upon
#' @param y The \code{GRanges} object of regions of interest
#' @param z The \code{numeric} or \code{integer} vector associated with the GRanges \code{x}
#' @param na.rm Logical indicating whether or not to include missing values in the results
#' @return A \code{numeric} vector of the same length as \code{y} with per-range means of \code{z}
#'
#' @export
#' 
#' @importFrom GenomicRanges findOverlaps
#' @importFrom IRanges viewMeans Views ranges as.matrix
#' @importFrom S4Vectors Rle
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
overlapMeans <- function(x, y, z, na.rm=FALSE) {
    stopifnot(class(x)=="GRanges")
    stopifnot(class(y)=="GRanges")
    stopifnot(class(z)=="numeric" | class(z)=="integer")
    stopifnot(length(x)==length(z))
    ov <- as.matrix(findOverlaps(y, x))
    ovMeans <- rep(as.numeric(NA), length(y))
    ovMeans[unique(ov[,1])] <- viewMeans(Views(z[ov[,2]], ranges(Rle(ov[,1]))), na.rm=na.rm)
    ovMeans
}


#' overlapSums
#'
#' Calculate the sum of values associated with ranges which overlap query ranges
#'
#' @param x The \code{GRanges} object to find overlaps upon
#' @param y The \code{GRanges} object of regions of interest
#' @param z The \code{numeric} or \code{integer} vector associated with the GRanges \code{x}
#' @param na.rm Logical indicating whether or not to include missing values in the results
#' @return A \code{numeric} vector of the same length as \code{y} with per-range sums of \code{z}
#'
#' @export
#' 
#' @importFrom GenomicRanges findOverlaps
#' @importFrom IRanges viewSums Views ranges
#' @importFrom S4Vectors Rle
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
overlapSums <- function(x, y, z, na.rm=FALSE) {
    stopifnot(class(x)=="GRanges")
    stopifnot(class(y)=="GRanges")
    stopifnot(class(z)=="numeric" | class(z)=="integer")
    stopifnot(length(x)==length(z))
    ov <- as.matrix(findOverlaps(y, x))
    ovSums <- rep(as.numeric(NA), length(y))
    ovSums[unique(ov[,1])] <- viewSums(Views(z[ov[,2]], ranges(Rle(ov[,1]))), na.rm=na.rm)
    ovSums
}


#' overlapRatios
#'
#' Calculate the ratio of C/cov values associated with ranges which overlap query ranges
#'
#' @param x The \code{GRanges} object to find overlaps upon
#' @param y The \code{GRanges} object of regions of interest
#' @param C The numerator vector associated with the GRanges \code{x}
#' @param cov The denominator vector associated with the GRanges \code{x}
#' @param minCov Minimum per-range coverage required to return a ratio
#' @param na.rm Logical indicating whether or not to include missing values in the results
#' @return A \code{numeric} vector of the same length as \code{y} with per-range ratios of C/cov of \code{z}
#'
#' @export
#' 
#' @author Aaron Statham <a.statham@@garvan.org.au>
overlapRatios <- function(x, y, C, cov, minCov=5, na.rm=FALSE) {
    C.sum <- overlapSums(x, y, values(x)[[C]], na.rm=na.rm)
    cov.sum <- overlapSums(x, y, values(x)[[cov]], na.rm=na.rm)
    rat <- C.sum/cov.sum
    rat[cov.sum<minCov] <- NA
    rat
}


#' coverageRatio
#'
#' Calculate the amount of each \code{subject} ranges overlap with \code{query} ranges
#'
#' @param query query ranges
#' @param subject subject ranges
#' @param ratio If \code{TRUE} return the ratio, if \code{FALSE} return the number of basepairs overlap
#' @return A \code{vector} of the same length as \code{subject} containing the calculated overlap
#'
#' @export
#' 
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicRanges strand findOverlaps reduce order coverage width
#' @importFrom IRanges viewSums Views
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
coverageRatio <- function(query, subject, ratio=TRUE) {
    if (length(unique(strand(subject)))>1) {
        warning("Supplied 'subject' contains strand information - this is ignored by coverageRatio")
        strand(subject) <- "*"
    }
    subject <- reduce(subject)
    # Fucking Views use Fucking RangesLists which do not fucking preserve order
    oo <- GenomicRanges::order(query)
    chr.oo <- unique(as.character(seqnames(query[oo])))
    seqlevels(query, force=TRUE) <- chr.oo
    query.cov <- viewSums(Views(coverage(subject)[chr.oo], as(query[oo], "RangesList")))
    query$result[oo] <- if (ratio) unlist(query.cov)/width(query[oo]) else unlist(query.cov)
    as.numeric(query$result)
}

#' overlapStats
#'
#' Calculate statistics of the overlaps of regions of interest and transcription factor binding sites
#'
#' @param ROI \code{GRanges} of interest
#' @param TFs \code{GRangesList} of TF binding sites
#' @param TF.gx Index matching \code{TFs} to rows of \code{gx}
#' @param gx \code{GRanges} of gene expression
#' @param accessible Estimate of the number of basepairs of the genome that is interrogatable
#' @return A \code{data.frame} with with a row per transcription factor in \code{TFs} containing the columns:
#'  \item{TF}{Name of the transcription factor}
#'  \item{TF.n}{Number of TF sites}
#'  \item{Hits}{Number of \code{ROI} that overlap TF sites}
#'  \item{Hits.percent}{Percent of \code{ROI} that overlap TF sites}
#'  \item{pval}{P-value of significance of overlap from the hypergeometric distribution}
#'  \item{adj.pval}{Adjusted p-value using "fdr"}
#'  \item{fold}{Estimated fold enrichment of \code{ROI} overlapping this TF compared to random chance}
#'  \item{Expression}{Expression result from \code{gx}}
#'  \item{RPKM.SA}{Expression result from \code{gx}}
#'  \item{RPKM.VA}{Expression result from \code{gx}}
#'  \item{logFC}{Expression result from \code{gx}}
#'  \item{FDR}{Expression result from \code{gx}}
#'
#' @export
#'
#' @importFrom parallel mclapply
#' @importFrom GenomicRanges width
#' @importFrom IRanges %over% elementLengths
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
overlapStats <- function(ROI, TFs, TF.gx, gx, accessible=2.5e9) {
    nTest <- accessible/median(width(ROI))
    ROI.ov <- simplify2array(mclapply(TFs, function(x) ROI %over% x))
    ROI.cov <- simplify2array(mclapply(TFs, function(x) sum(coverageRatio(ROI, x, FALSE))))
    ROI.stats <- data.frame("TF"=names(TFs), 
                          "TF.n"=elementLengths(TFs), 
                          "Hits"=colSums(ROI.ov), 
                  "Hits.percent"=colMeans(ROI.ov)*100, 
                stringsAsFactors=FALSE)
    ROI.stats$pval <- sapply(names(TFs), function(s) 
        phyper(q=ROI.stats[s,]$Hits-1,
               m=ROI.stats[s,]$TF.n,
               n=nTest-ROI.stats[s,]$TF.n, 
               k=length(ROI), lower.tail=FALSE))
    ROI.stats$adj.pval <- p.adjust(ROI.stats$pval, "fdr")
    
    ROI.stats$fold <- (ROI.cov/sum(width(ROI)))/(sum(width(ROI))/accessible)*(sum(width(TFs))/accessible)
    ROI.stats$Expression <- gx$result[TF.gx]
    ROI.stats$RPKM.SA <- gx$SA[TF.gx]
    ROI.stats$RPKM.VA <- gx$VA[TF.gx]
    ROI.stats$logFC <- gx$logFC[TF.gx]
    ROI.stats$FDR <- gx$FDR[TF.gx]
    ROI.stats
}

