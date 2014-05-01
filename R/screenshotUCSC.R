#' screenshotUCSC
#'
#' Load a UCSC track hub and save a screenshot of a given region as a pdf
#'
#' @param url UCSC URL to use, sets genome and screenshot resolution
#' @param hubfile URL of UCSC track hub to load
#' @param chr Chromosome to screenshot
#' @param start Start of the region to screenshot
#' @param end End of the region to screenshot
#' @param filename Filename of the output pdf
#' @param session URL of the UCSC session to load (optional)
#' @return Called for the side effect of downloading a pdf screenshot, returns invisible debugging output
#'
#' @export
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
screenshotUCSC <- function(url, hubfile, chr, start, end, filename, session=NULL) {
    oldpen <- options("scipen")
    options(scipen=100)
    geturl <- paste(url, "&hgt.psOutput=on&hubUrl=", hubfile, "&position=",chr,":",start,"-",end, sep="")
    if (!is.null(session)) geturl <- paste(geturl, "&hgS_loadUrlName=", session, "&hgS_doLoadUrl=submit", sep="")
    temp <- readLines(geturl)
    pdfurl <- paste("http://genome.ucsc.edu/trash", gsub(".*../trash","",gsub(".pdf.*", "", temp[grep("the current browser graphic in PDF", temp, fixed=TRUE)])), ".pdf", sep="")
    options(scipen=oldpen)
    download.file(pdfurl, filename, mode="wb", quiet=TRUE)	
    Sys.sleep(2)
    invisible(list(temp=temp, pdfurl=pdfurl))
}

#' plotManyUCSC
#'
#' Takes many UCSC screenshots and collates them together into a single pdf using the system ghostscript
#'
#' @param x \code{GRanges} of regions to screenshot
#' @param outFile Filename of the output pdf
#' @param mc.cores No of threads to use to take screenshots - too many will break UCSC/this function
#' @param hub URL of UCSC track hub to load
#' @param session URL of the UCSC session to load (optional)
#' @export
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
plotManyUCSC <- function(x, outFile, mc.cores=1, hub, session=NULL) {
    plotDir <- paste(tempdir(), "/plot/", sep="")
    unlink(plotDir, recursive=TRUE, force=TRUE)
    dir.create(plotDir)

    x <- as.data.frame(x)
    x.tmp <- mclapply(1:nrow(x), function(i)
        screenshotUCSC("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&pix=1800",
            hub,
            x$seqnames[i],
            x$start[i],
            x$end[i],
            paste(plotDir,  i, ".pdf", sep=""),
            session)
        , mc.cores=1)
    inFiles <- paste(plotDir, 1:nrow(x), ".pdf", sep="", collapse=" ")
    system(paste("gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=", outFile, " ", inFiles, sep=""))
    unlink(plotDir, recursive=TRUE, force=TRUE)
    invisible(x.tmp)
}

#' plotPromoterDMRs
#'
#' Takes UCSC screenshots which contain a region of interest, its closest TSS and (optionally) its closest CpG island
#'
#' @param x \code{GRanges} of regions to screenshot
#' @param tx \code{GRanges} of transcript expression
#' @param CpGislands \code{GRanges} of CpG island positions
#' @param outFile Filename of the output pdf
#' @param zoom The distance to pad either side of the regions of interest by
#' @param proteinOnly Whether to restrict closest TSS plotting to protein coding transcripts only, or to include non-coding transcripts
#' @param CpGi Whether to include the closest CpG island in the screenshot
#' @param mc.cores No of threads to use to take screenshots - too many will break UCSC/this function
#' @param hub URL of UCSC track hub to load
#' @param session URL of the UCSC session to load (optional)
#'
#' @export
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
plotPromoterDMRs <- function(x, tx, CpGislands, outFile, zoom=3000, proteinOnly=FALSE, CpGi=TRUE, mc.cores=1, hub, session=NULL) {
    x <- x[order(abs(x$meanDiff), decreasing=TRUE)]
    x.TSS <- resize(tx[match(if (proteinOnly) x$tx_name_prot else x$tx_name, tx$tx_name)], 1, fix="start")

    # come up with the window to screenshot which includes the TSS, the DMR 
    if (CpGi) { # and the island
        start(x) <- pmin(start(x)-zoom, start(x.TSS)-zoom, start(CpGislands)[x$CpGi]-zoom)
        end(x) <- pmin(end(x)+zoom, end(x.TSS)+zoom, end(CpGislands)[x$CpGi]+zoom)
    } else { # or not the island
        start(x) <- pmin(start(x)-zoom, start(x.TSS)-zoom)
        end(x) <- pmin(end(x)+zoom, end(x.TSS)+zoom)
    }
    plotManyUCSC(x, outFile=outFile, mc.cores=mc.cores, hub=hub, session=session)
}

