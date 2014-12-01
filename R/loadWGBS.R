#' loadWGBS
#'
#' Loads a WGBS experiments bigTable and samples file and checks consistency
#'
#' @param path The path to the project.csv and bigTable.tsv
#' @param build BSgenome object of the organism the data is from
#' @return A list with the elements 'samples', 'CpGs' and if a NOMe experiment 'GCH'
#'
#' @export
#'
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges values
#' @importFrom GenomeInfoDb seqlevels seqlengths
#' @importFrom IRanges IRanges
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
loadWGBS <- function(path, build=NULL) {
	stopifnot(all(file.exists(paste0(path, c("/bigTable.tsv", "/project.csv")))))
	NOMe <- file.exists(paste0(path, "/bigTable.GCH.tsv"))

	# Load project.csv
	samples <- read.csv(paste0(path, "/project.csv"), stringsAsFactors=FALSE)
	stopifnot(all(!duplicated(samples$Sample)))
	samples$C <- paste0(samples$Sample, ".C")
	samples$cov <- paste0(samples$Sample, ".cov")

	# Load CpG bigTable
	tab <- fread(paste0(path, "/bigTable.tsv"))
	CpGs <- GRanges(tab$chr, IRanges(tab$position, width=1))
	values(CpGs) <- as.data.frame(tab)[,-c(1:2)]
	# fix for sample names that start with a number
	names(values(CpGs)) <- colnames(tab)[-c(1:2)]
	if (!is.null(build)) seqlengths(CpGs) <- seqlengths(build)[seqlevels(CpGs)]
	rm(tab)

	# Sanity checks
	stopifnot(all(samples$C %in% names(values(CpGs))))
	stopifnot(all(samples$cov %in% names(values(CpGs))))

	# Load NOMe GCH table
	if (NOMe) {
	    tab <- fread(paste0(path, "/bigTable.GCH.tsv"))
	    GCH <- GRanges(tab$chr, IRanges(tab$position, width=1))
	    values(GCH) <- as.data.frame(tab)[,-c(1:2)]
	    # fix for sample names that start with a number
	    names(values(GCH)) <- colnames(tab)[-c(1:2)]
	    if (!is.null(build)) seqlengths(GCH) <- seqlengths(build)[seqlevels(GCH)]
	    rm(tab)

	    # Sanity checks
	    stopifnot(all(samples$C %in% names(values(GCH))))
	    stopifnot(all(samples$cov %in% names(values(GCH))))
	}

	if (NOMe) list(samples=samples, CpGs=CpGs, GCH=GCH) else list(samples=samples, CpGs=CpGs)
}
