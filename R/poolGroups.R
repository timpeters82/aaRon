#' poolGroups
#'
#' Pools methylation count data of samples from the same group
#'
#' @param x \code{GRanges} of methylation count data
#' @param samples \code{data.frame} describing the samples and their grouping
#' @return A \code{list} containing the elements \code{CpGs} and \code{samples} with the groups pooled
#'
#' @export
#'
#' @importFrom IRanges as.matrix
#' @importFrom GenomicRanges values
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
poolGroups <- function(x, samples) {
    groups <- unique(samples$Group)
    tmp <- unvalue(x)
    x <- as.matrix(values(x))
    for (i in groups) {
        values(tmp)[[paste0(i, ".C")]] <- rowSums(x[, samples$C[samples$Group==i], drop=FALSE])
        values(tmp)[[paste0(i, ".cov")]] <- rowSums(x[, samples$cov[samples$Group==i], drop=FALSE])
    }
    samples <- data.frame(Sample=groups, Group=groups, C=paste0(groups, ".C"), cov=paste0(groups, ".cov"), stringsAsFactors=FALSE)
    list(CpGs=tmp, samples=samples)
}
