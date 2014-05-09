#' overlapPositions
#'
#' Returns the positions where \code{subject} ranges overlap the \code{query}
#'
#' @param query \code{GRanges} of query regions
#' @param subject \code{GRanges} of subject regions (must be of width 1)
#' @return Positions where \code{subject} ranges overlap the \code{query} ranges
#'
#' @export
#'
#' @importFrom GenomicRanges width findOverlaps start
#' @importFrom IRanges as.matrix
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
overlapPositions <- function(query, subject) {
    stopifnot(all(width(subject)==1))
    mm <- as.matrix(findOverlaps(query, subject))
    start(subject)[mm[,2]]-start(query)[mm[,1]]
}

#' boundaryOverlaps
#'
#' Returns 
#'
#' @param query \code{GRanges} of query regions
#' @param subject \code{GRanges} of subject regions
#' @param boundary Whether to interrogate the start, end or both boundaries of the \code{query} range
#' @param distance Number of basepairs up and downstream of \code{subject} regions to interrogate
#' @return Positions where \code{subject} ranges overlap the boundaries of \code{query} ranges
#'
#' @export
#'
#' @importFrom GenomicRanges width strand
#' @importFrom IRanges as.matrix
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
boundaryOverlaps <- function(query, subject, boundary=c("both", "start", "end"), distance=10000) {
    stopifnot(all(width(subject)==1))
    boundary <- match.arg(boundary)
    query.boundary <- GRanges()

    # start boundaries
    if (boundary %in% c("both", "start")) {
        tmp <- resize(resize(query, 1, fix="start"), distance*2, fix="center")
        strand(tmp) <- ifelse(as.character(strand(tmp)) %in% c("*", "+"), "+", "-")
        query.boundary <- c(query.boundary, tmp)
    }

    # end boundaries
    if (boundary %in% c("both", "end")) {
        tmp <- resize(resize(query, 1, fix="end"), distance*2, fix="center")
        strand(tmp) <- ifelse(as.character(strand(tmp)) %in% c("*", "+"), "-", "+")
        query.boundary <- c(query.boundary, tmp)
    }

    # Remove ranges resized because of chromosome ends as coordinates may not be accurate
    query.boundary <- query.boundary[width(query.boundary)==distance*2]

    # Find overlap postions
    mm <- as.matrix(findOverlaps(query.boundary, subject))
    bp <- start(subject)[mm[,2]]-start(query.boundary)[mm[,1]]

    # Adjust to be relative to the boundary site
    bp <- ifelse(as.character(strand(query.boundary))[mm[,1]]=="+", bp-distance, distance-bp)
    bp
}

