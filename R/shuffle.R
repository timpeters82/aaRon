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

#' shuffleBed
#' 
#' A hack to use bedtools shuffleBed to shuffle ranges across a genome
#'
#' @param x Set of ranges to be shuffled
#' @return A shuffled \code{GRanges} of the same length as \code{x}
#'
#' @export
#'
#' @importFrom GenomicRanges seqlengths seqnames seqlevels start end width
#' @importFrom rtracklayer import.bed
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
shuffleBed <- function(x) {
    if (nchar(Sys.which("shuffleBed"))==0) stop("shuffleBed must be in the PATH")
    if (any(is.na(seqlengths(x)))) stop("seqlengths of x must be defined")
    genome_file <- tempfile()
    writeLines(paste(names(seqlengths(x)), seqlengths(x), sep="\t"), genome_file)
    
    bed_file <- tempfile()
    x <- x[order(width(x), decreasing=TRUE)] #makes it easier to find non-overlaps
    writeLines(paste(as.character(seqnames(x)), start(x), end(x), sep="\t"), bed_file)
    
    shuffle_file <- tempfile()
    system(paste("shuffleBed -maxTries 10000 -noOverlapping -i", bed_file, "-g", genome_file, ">", shuffle_file))
    shuffled <- import.bed(shuffle_file)
    seqlevels(shuffled, force=TRUE) <- seqlevels(x)
    seqlengths(shuffled) <- seqlengths(x)
    
    #clean up
    unlink(c(genome_file, bed_file, shuffle_file))
    shuffled
}
