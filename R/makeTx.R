#' makeTx
#'
#' Create a tx object for annotating regions from an Ensembl GTF
#'
#' @param file Filename to uncompressed Ensembl GTF to process
#' @param genome \code{BSgenome} object the GTF will be associated with
#' @param type_attrib Name of the GTF attribute containing the type of transcript - usually "transcript_biotype" but for Homo_sapiens "source"
#' @return \code{GRanges} of transcripts with appropriate columns for \code{annotateRegions} namely \code{gene_name}, \code{gene_id}, \code{tx_id} and \code{tx_type}
#'
#' @export
#'
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges 
#' @importFrom IRanges IRanges
#' @importFrom stringr str_extract
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
makeTx <- function(file, genome, type_attrib=c("transcript_biotype", "source")) {
    stopifnot(length(file)==1)
    stopifnot(file.exists(file))
    type_attrib <- match.arg(type_attrib)

    gtf <- fread(file)
    gtf <- gtf[gtf$V3=="transcript", ]
    # UCSC style seqnames
    gtf$V1 <- paste0("chr", gtf$V1)
    gtf$V1[gtf$V1=="chrMT"] <- "chrM"
    # Fix scaffolds
    dots <- grep("\\.1$", gtf$V1)
    gtf$V1[dots] <- gsub("\\.1$", "", gsub("^chr", "chrUn_", gtf$V1[dots]))
    tx <- GRanges(gtf$V1, IRanges(gtf$V4, gtf$V5), strand=gtf$V7)

    # collect attributes
    attribs <- c("gene_name"="gene_name", "gene_biotype", "gene_id"="gene_id", "tx_id"="transcript_id")
    if (type_attrib=="transcript_biotype") {
        if (!any(grepl(type_attrib, gtf$V9[1:1000]))) stop("'transcript_biotype' attribute not found in supplied GTF file!")
        attribs <- c(attribs, "transcript_biotype"="transcript_biotype")
    }
    gtf.attribs <- data.frame(lapply(attribs, function(a) gsub('.*"', '', gsub('";', '', str_extract(gtf$V9, paste0(a, ".+?;"))))), stringsAsFactors=FALSE) 
    if (type_attrib=="source") gtf.attribs$tx_type <- gtf$V2
    NAs <- is.na(gtf.attribs$gene_name)
    gtf.attribs$gene_name[NAs] <- gtf.attribs$gene_id[NAs]

    values(tx) <- gtf.attribs

    # Fix seqnames and seqlengths
    if (!all(seqlevels(tx) %in% seqlevels(genome)))
        warning("Some seqlevels in supplied gtf not present in genome and will be dropped.")
    n.tx <- length(tx)
    seqlevels(tx, force=TRUE) <- seqlevels(genome)
    if (n.tx!=length(tx)) warning(n.tx-length(tx), " rows dropped!")
    seqlengths(tx) <- seqlengths(genome)

    tx
}
