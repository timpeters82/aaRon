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

plotManyUCSC <- function(x, outFile, zoom=0, mc.cores=1, hub, session=NULL) {
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

plotPromoterDMRs <- function(x, outFile, zoom=3000, proteinOnly=FALSE, CpGi=TRUEhub, session=NULL) {
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
    plotManyUCSC(x, outFile=outFile, hub=hub, session=session)
}

