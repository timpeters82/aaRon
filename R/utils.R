strip <- function(x) {
    if (class(x)=="GRangesList") return(endoapply(x, strip))
    x <- x[seqnames(x) %in% names(Hsapiens)[1:24]]
    seqlevels(x) <- names(Hsapiens)[1:25]
    values(x) <- NULL
    strand(x) <- "*"
    reduce(x)
}

unvalue <- function(x) {
    values(x) <- NULL
    x
}

unstrand <- function(x) {
    strand(x) <- "*"
    x
}

