#' overlapPositions
#'
#' Returns the positions where \code{subject} ranges overlap the \code{query}
#'
#' @param query \code{GRanges} of query regions
#' @param subject \code{GRanges} of subject regions (must be of width 1 and unstranded)
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
    stopifnot(all(strand(subject)=="*"))
    mm <- as.matrix(findOverlaps(query, subject))
    start(subject)[mm[,2]]-start(query)[mm[,1]]
}

#' overlapBoundaries
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
#' @importFrom GenomicRanges width strand strand<-
#' @importFrom IRanges as.matrix
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
overlapBoundaries <- function(query, subject, distance=10000, boundary=c("both", "start", "end")) {
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

#' overlapRegions
#'
#' Returns the relative position of \code{subject} sites within \code{query} regions
#'
#' @param query \code{GRanges} of query regions
#' @param subject \code{GRanges} of subject regions (must be of width 1 and unstranded)
#' @return \code{data.frame} of index and position
#'
#' @export
#'
#' @importFrom GenomicRanges width strand strand<- findOverlaps
#' @importFrom IRanges subjectHits queryHits
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
overlapRegions <- function(query, subject) {
    stopifnot(all(width(subject)==1))
    stopifnot(all(strand(subject)=="*"))
    if (!(all(strand(query)=="*") | all(strand(query)!="*")))
        stop("Query ranges must be either all stranded or all unstranded, not a mixture")
    strand(query)[strand(query)=="*"] <- "+"
    ov <- findOverlaps(query, subject)
    res <- data.frame(index=subjectHits(ov))
    # correct for strand
    res$position <- ifelse(as.character(strand(query))[queryHits(ov)]=="+", overlapPositions(query, subject),
                    width(query[queryHits(ov)])-overlapPositions(query, subject))
    if (!is.null(query$ref)) {
        res$position <- res$position-(query$ref[queryHits(ov)]-start(query)[queryHits(ov)])
    }
    res
}

#' overlapRegionFlanks
#'
#' Returns the relative position \code{subject} sites around \code{query} regions \code{start} positions if stranded, or around \code{query} region centres if unstranded
#'
#' @param query \code{GRanges} of query regions
#' @param subject \code{GRanges} of subject regions (must be of width 1 and unstranded)
#' @param up Distance upstream to find overlaps
#' @param down Distance downstream to find overlaps
#' @return \code{data.frame} of index and position
#'
#' @export
#'
#' @importFrom GenomicRanges strand resize start end promoters
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
overlapRegionFlanks <- function(query, subject, up=1000, down=1000) {
    stopifnot(up>0 & down>0)
    stopifnot(all(strand(subject)=="*"))
    if (!(all(strand(query)=="*") | all(strand(query)!="*")))
        stop("Query ranges must be either all stranded or all unstranded, not a mixture")

    # Create flanking regions
    if (all(strand(query)=="*")) { # Unstranded
        ROI <- resize(query, 1, fix="center")
        ROI$ref <- start(ROI)
        strand(ROI) <- "+"
        ROI <- promoters(ROI, up, down)
    } else {
        ROI <- promoters(query, up, down)
        ROI$ref <- ifelse(as.character(strand(query))=="+", start(query), end(query))
    }
    overlapRegions(ROI, subject)
}
