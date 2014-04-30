#tabulate how much of the genome is each annotation
tabulateCoverage <- function(reg1, background, reg2) {
    require(parallel)
    reg1 <- c(GRangesList("background"=background), reg1)
    temp <- mclapply(1:length(reg1), function(x) sapply(reg2, function(y) {
        sum(as.numeric(width(intersect(reg1[[x]],y))))
    }))
    names(temp) <- names(reg1)
    temp <- do.call(cbind, temp)
    temp <- cbind(temp, "genome"=sapply(reg2, function(x) sum(as.numeric(width(x)))))
    temp <- rbind(temp, "total"=c(sapply(reg1, function(x) sum(as.numeric(width(x)))), NA))
    tempRat <- t(t(temp)/temp[nrow(temp),])*100
    tempRat2 <- tempRat[,-1]/tempRat[,"background"]
    list("coverage"=temp, "ratios"=tempRat, "observed/expected"=tempRat2[-nrow(tempRat2),-ncol(tempRat2), drop=FALSE])
}

#create a barplot of enrichment/depletion of regions vs annotations, assess significance
#by permuting regions
coverageBarplot <- function(regions, background, annotations, main, nperm=1000, cols=NULL, colBy=c("regions", "annotations"), verbose=FALSE) {
    colBy <- match.arg(colBy)
    if (colBy=="regions") 
        if (is.null(cols)) cols=rainbow(length(regions)) else stopifnot(length(cols)==length(regions))
    else
        if (is.null(cols)) cols=rainbow(length(annotations)) else stopifnot(length(cols)==length(annotations))
    #create sampled regions to assess significance
    nregions <- length(background)

    perms <- lapply(regions, function(x) {
        tmp <- sample(nregions, length(x)*nperm, replace=TRUE)
        tmp <- GRanges(background@seqnames[tmp], background@ranges[tmp])
        split(tmp, rep(1:nperm, each=length(x)))
    })
    if (verbose) message("Created permutations")

    covTable <- tabulateCoverage(regions, background, annotations)
    if (verbose) message("Tabulated coverage")

    permTable <- lapply(perms, tabulateCoverage, background, annotations)
    if (verbose) message("Tabulated coverage on permutations")
    permQuant <- lapply(permTable, function(x) apply(x$o, 1, quantile, c(0.025, 0.975)))
    permPval <- sapply(names(permTable), function(x)
        sapply(rownames(permTable[[x]]$o), function(y) min(
        (1-sum(covTable$o[y,x]>permTable[[x]]$o[y,])/nperm)*2,
        (1-sum(covTable$o[y,x]<permTable[[x]]$o[y,])/nperm)*2)))
    if (verbose) message("Calculated significance")

    #make barplot
    bp <- barplot(t(covTable$o), border=NA, beside=TRUE, col=cols, las=2, ylab="Observed/Expected", main=main, ylim=c(0, ceiling(max(covTable$o))))
    if (length(regions)==1) bp <- matrix(bp, nrow=1)
    for (i in 1:length(regions)) points(bp[i,], covTable$o[,i]+0.1, pch=ifelse(covTable$o[,i]<permQuant[[i]][1,] | covTable$o[,i]>permQuant[[i]][2,], "*", ""), cex=3)

    abline(h=c(0,1), col=c("black", "red"), lwd=c(1,3))
    if (colBy=="regions") {
        legend("topleft", fill=cols, legend=names(regions))
    }
    #return all the data
    list("covTable"=covTable, "permTable"=permTable, "permQuant"=permQuant, "perms"=perms, "permPval"=permPval)
}

