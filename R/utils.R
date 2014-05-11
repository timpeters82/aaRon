#' strip
#'
#' Cleans up a \code{GRanges} object
#'
#' @param x A \code{GRanges} objext
#' @param genome \code{BSgenome} object
#' @param chrs Chromosomes to limit \code{x} to
#' @return Reduced \code{GRanges} of restricted chromosomes, with no strand or element metadata
#'
#' @export
#'
#' @importFrom GenomicRanges seqnames seqlevels reduce
#' @importFrom IRanges endoapply
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
strip <- function(x, genome, chrs=1:24) {
    if (class(x)=="GRangesList") return(endoapply(x, strip))
    seqlevels(x, force=TRUE) <- seqnames(genome)[chrs]
    reduce(unstrand(unvalue(x)))
}

#' unvalue
#'
#' Remove element metadata from a GRanges object
#'
#' @param x A \code{GRanges} objext
#' @return \code{x} with no element metadata
#'
#' @export
#'
#' @importFrom GenomicRanges values
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
unvalue <- function(x) {
    values(x) <- NULL
    x
}

#' unstrand
#'
#' Remove strand data from a GRanges object
#'
#' @param x A \code{GRanges} objext
#' @return \code{x} with strand set to '*'
#'
#' @export
#'
#' @importFrom GenomicRanges strand
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
unstrand <- function(x) {
    strand(x) <- "*"
    x
}

