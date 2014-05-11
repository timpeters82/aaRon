#' shuffle
#'
#' Shuffle a set of \code{GRanges} across an available genome
#'
#' @param x Set of ranges to be shuffled
#' @return A shuffled \code{GRanges} of the same length as \code{x}
#'
#' @export
#'
#' @importFrom GenomicRanges width findOverlaps start
#' @importFrom IRanges as.matrix
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
shuffle <- function(x) {
    if (any(is.na(seqlengths(x)))) stop("seqlengths of x must be defined")
    x <- unvalue(x)
    available <- GRanges(seqlevels(x), IRanges(1, seqlengths(x)), seqlengths=seqlengths(x))
    x.random <- sample(sum(as.numeric(width(x))), length(x))
    x.pos <- projectOnto(x.random, available)
    suppressWarnings(width(x.pos) <- width(x))
    x.pos
}

projectOnto <- function(pos, avail) {
    offset <- c(0, cumsum(as.numeric(width(avail))))
    hits <- hits <- sapply(pos, function(x) match(TRUE, x<offset)-1)
    GRanges(seqnames(avail)[hits], 
            IRanges(start(avail)[hits]+(pos-offset[hits])-1, width=1),
            seqlengths=seqlengths(avail))
}
