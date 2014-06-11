#' binarizeFromPeaks
#'
#' Converts bed files of peaks (eg from MACS or PeakRanger) to binary files for ChromHMM
#'
#' @param peaks \code{data.frame} with column names \code{sample}, \code{mark}, \code{bed}
#' @param outdir Path to directory to output binary files to
#' @param resolution Resolution in basepairs to bin the genome into
#' @param seqlen Named \code{character} vector of chromosome lengths
#' @return Called for side effect of creating binary files
#'
#' @export
#' 
#' @importFrom Repitools genomeBlocks
#' @importFrom rtracklayer import.bed
#' @importFrom IRanges %over%
#' @importFrom GenomicRanges seqnames
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
binarizeFromPeaks <- function(peaks, outdir, resolution=200, seqlen) {
	stopifnot(all(c("sample", "mark", "bed") %in% names(peaks))) # Required columns
	stopifnot(all(file.exists(peaks$bed))) # Check beds exist
	stopifnot(file.info(outdir)$isdir) # Check outdir exists/is a dir
	stopifnot(length(seqlen)>0 & min(seqlen>0)) # Check seqlengths are supplied and sane

	# Split peaks by sample
	peaks <- peaks[order(peaks$sample, peaks$mark),]
	peaks <- split(peaks, peaks$sample)

	# Ensure all samples have all marks
	peaks.marks <- sapply(peaks, "[[", "mark")
	stopifnot(all(apply(peaks.marks, 2, function(x) identical(peaks.marks[,1], x))))

	# Genomic tiles to be binarized
	tiles <- genomeBlocks(seqlen, width=resolution)
	for (i in 1:length(peaks)) {
		tiles.bin <- suppressWarnings(simplify2array(lapply(peaks[[i]]$bed, function(x) as.integer(tiles %over% import.bed(pipe(paste0("cut -f1,2,3 ", x)))))))
		colnames(tiles.bin) <- peaks[[i]]$mark
		tmp <- split(as.data.frame(tiles.bin), as.character(seqnames(tiles)))
		# Export per chromosome
		for (j in names(tmp)) {
			f <- file(paste0(outdir, "/", names(peaks)[i], "_", j, "_binary.txt"), open="wt")
			writeLines(paste(names(peaks)[i], j, sep="\t"), f)
			write.table(tmp[[j]], f, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
			close(f)
		}
	}

}
