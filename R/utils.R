#' strip
#'
#' Cleans up a \code{GRanges} object
#'
#' @param x A \code{GRanges} objext
#' @return Reduced \code{GRanges} with no strand or element metadata
#'
#' @export
#'
#' @importFrom GenomicRanges reduce
#' @importFrom IRanges endoapply
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
strip <- function(x) {
    if (class(x)=="GRangesList") return(endoapply(x, strip))
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
#' @importFrom GenomicRanges values<-
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
#' @importFrom GenomicRanges strand<-
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
unstrand <- function(x) {
    strand(x) <- "*"
    x
}

