#' annotateRegionsByInteractions
#'
#' Annotate regions with respect to gene annotations using chromatin interaction data
#'
#' @param reg The \code{GRanges} to be annotated
#' @param tx \code{GRanges} of transcripts, annotated with \code{gene_name}, \code{gene_id}, \code{tx_id} and \code{tx_type}
#' @param interactions \code{GRanges} of interaction data
#' @param distance Distance interaction can be from a promoter to count as promoter linked
#' @return \code{reg} with additional per region metadata
#'
#' @export
#'
#' @importFrom GenomicRanges findOverlaps distanceToNearest resize setdiff flank
#' @importFrom IRanges IRanges
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
annotateRegionsByInteractions <- function(reg, tx, interactions, distance=2000) {
    if (!all(c("gene_name", "gene_id", "tx_id", "tx_type") %in% names(values(tx)))) stop("Supplied tx does not contain all required columns, was it created by makeTx?")

    if (!(all(as.character(seqnames(reg)) %in% seqlevels(tx)))) {
        message("Removing seqlevels/regions from 'reg' that are not represented present in supplied 'tx'")
        seqlevels(reg, force=TRUE) <- seqlevels(tx)
    }

    TSS <- resize(resize(tx, 1, fix="start"), distance*2, fix="center")

    # annotate interactions with promoters within distance
    message("Annotating interaction fragments with TSSs")
    tmp <- rep(list(integer(0)), length(interactions))
    int.TSS <- as.matrix(findOverlaps(interactions, TSS))
    int.TSS <- split(int.TSS[,2], int.TSS[,1])
    tmp[as.integer(names(int.TSS))] <- unname(int.TSS)
    interactions$TSSs <- IntegerList(tmp)
    rm(int.TSS, tmp)

    # Overlap regions and interactions
    reg.TSSs <- reg.frags <- reg.frag.TSSs <- reg.pairs <- reg.pair.TSSs  <- rep(list(integer(0)), length(reg))

    # TSSs directly overlapped by reg
    message("Finding regions directly overlapping a TSS")
    tmp <- as.matrix(findOverlaps(reg, TSS))
    tmp <- split(tmp[,2], tmp[,1])
    reg.TSSs[as.integer(names(tmp))] <- unname(tmp)
    reg$TSSs <- IntegerList(reg.TSSs)
    rm(reg.TSSs, tmp)

    # frags overlapped by reg
    message("Finding regions overlaping an interaction fragment")
    tmp <- as.matrix(findOverlaps(reg, interactions))
    tmp <- split(tmp[,2], tmp[,1])
    reg.frags[as.integer(names(tmp))] <- unname(tmp)
    reg$frags <- IntegerList(reg.frags)
    rm(reg.frags)

    # Extracts lists of subsets of IntegerLists, unlists and uniques each subset individuall
    # Could do with another rewrite for speed - full unlist then split per calculated subsets range then unique to avoid an unlist per subset
    extractUnique <- function(x, indList) {
        tmp <- data.frame("i"=rep(names(indList), elementLengths(indList)), "inds"=unname(unlist(indList)))
        lapply(split(x[tmp$inds], tmp$i), function(i) unique(unlist(i)))[names(indList)]
    }

    # TSS overlapping frags overlapped by reg
    message(" * Annotating whether the fragments overlap a TSS")
    tmp.TSSs <- extractUnique(interactions$TSSs, tmp)
    reg.frag.TSSs[as.integer(names(tmp.TSSs))] <- unname(tmp.TSSs)
    reg$frag.TSSs <- IntegerList(reg.frag.TSSs)
    rm(reg.frag.TSSs, tmp.TSSs)

    # Pairs of frags overlapped by reg
    message(" * Linking the fragments to their distal pairs")
    tmp.pairs <- extractUnique(interactions$pair, tmp)
    reg.pairs[as.integer(names(tmp.pairs))] <- unname(tmp.pairs)
    reg$pairs <- IntegerList(reg.pairs)
    rm(reg.pairs, tmp)

    # TSSs overlapping pairs of frags overlapped by reg
    message(" * Finding distal pairs that overlap a TSS")
    tmp.pair.TSSs <- extractUnique(interactions$TSSs, tmp.pairs)
    reg.pair.TSSs[as.integer(names(tmp.pair.TSSs))] <- unname(tmp.pair.TSSs)
    reg$pair.TSSs <- IntegerList(reg.pair.TSSs)
    rm(reg.pair.TSSs, tmp.pairs, tmp.pair.TSSs)

    # Try to summarize peak
    message("Summarising overlaps")
    reg$Summary <- ifelse(elementLengths(reg$TSSs)>0, "TSS", 
                        ifelse(elementLengths(reg$pair.TSSs)>0, "TSS_Interaction",
                            ifelse(elementLengths(reg$pairs)>0, "Non_TSS_Interaction", "Non_Interacting")))

    reg
}

