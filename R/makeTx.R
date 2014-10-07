#' makeTx
#'
#' Create a tx object for annotating regions from an Ensembl GTF
#'
#' @param file Filename to uncompressed Ensembl GTF to process
#' @param genome \code{BSgenome} object the GTF will be associated with
#' @param type_attrib Name of the GTF attribute containing the type of transcript - usually "transcript_biotype" or "gene_biotype"
#' @return \code{GRanges} of transcripts with appropriate columns for \code{annotateRegions} namely \code{gene_name}, \code{gene_id}, \code{tx_id} and \code{tx_type}
#'
#' @export
#'
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges 
#' @importFrom GenomeInfoDb seqlevelsStyle<-
#' @importFrom IRanges IRanges
#' @importFrom stringr str_extract
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
makeTx <- function(file, genome, type_attrib="transcript_biotype") {
    stopifnot(length(file)==1)
    stopifnot(file.exists(file))

    gtf <- fread(file)
    gtf <- gtf[gtf$V3=="transcript", ]
    # UCSC style seqnames
    gtf$V1 <- paste0("chr", gtf$V1)
    gtf$V1[gtf$V1=="chrMT"] <- "chrM"
    tx <- GRanges(gtf$V1, IRanges(gtf$V4, gtf$V5), strand=gtf$V7)

    # collect attributes
    attribs <- c("gene_name"="gene_name", "gene_id"="gene_id", "tx_id"="transcript_id", "tx_type"=type_attrib)
    gtf.attribs <- data.frame(lapply(attribs, function(a) gsub('.*"', '', gsub('";', '', str_extract(gtf$V9, paste0(a, ".+?;"))))), stringsAsFactors=FALSE) 

    NAs <- is.na(gtf.attribs$gene_name)
    gtf.attribs$gene_name[NAs] <- gtf.attribs$gene_id[NAs]

    values(tx) <- gtf.attribs

    # Fix seqnames and seqlengths
    seqlevels(tx, force=TRUE) <- seqlevels(genome)[seqlevels(genome) %in% seqlevels(tx)]
    seqlengths(tx) <- seqlengths(genome)[seqlevels(tx)]

    tx
}
