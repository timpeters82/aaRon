#' ampliconGRanges
#' 
#' 
#'
#' @param x \code{data.frame} of amplicons containing the columns "Amplicon", "Target", "FW", "RV", "Sequenom"
#' @param genome \code{BSgenome} to map the amplicons to
#' @param mc.cores Number of cores to use during mapping
#' @return \code{GRanges} of amplicon locations with metadata
#'
#' @export
#' 
#' @importFrom parallel mclapply
#' @importFrom Biostrings vmatchPattern
#' @importFrom IRanges elementLengths
#' @importFrom GenomicRanges unlist GRangesList values width start resize
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
ampliconGRanges <- function(x, genome, mc.cores=8) {
	x$FW_len <- ifelse(x$Sequenom, nchar(x$FW)-10, nchar(x$FW))
	x$RV_len <- ifelse(x$Sequenom, nchar(x$RV)-31, nchar(x$RV))

	#
	xHits <- mclapply(x$Target, vmatchPattern, genome, mc.cores=mc.cores)
	stopifnot(all(elementLengths(xHits)==1))
	#Keep x with a single hit (for now)
	x <- x[elementLengths(xHits)==1,]
	xHits <- xHits[elementLengths(xHits)==1]
	xHits <- unlist(GRangesList(xHits))
	values(xHits) <- x
	x <- xHits
	x$size <- width(x)-x$FW_len-x$RV_len
	x$internal <- substr(x$Target, x$FW_len+1, x$FW_len+x$size)
	rm(xHits)

	#annotate amplicon CG sites
	x$CGs <- mapply(function(cgs, offset) cgs+offset-1, gregexpr("CG", x$internal), start(x)+x$FW_len)

	#Add a GRangesList of primer positions
	x$primers <- GRangesList(lapply(1:length(x), function(i) {
	    tmp <- rep(x[i],2)
	    tmp <- resize(tmp, c(tmp$FW_len[1], tmp$RV_len[1]), fix=c("start", "end"))
	    values(tmp) <- NULL
	    tmp
	}))

	x
}
