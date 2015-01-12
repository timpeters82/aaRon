#' ampliconGRanges
#' 
#' Align a \code{data.frame} describing amplicons against a supplied genome, returning an annotated amplicon \code{GRanges}
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
ampliconGRanges <- function(x, genome, mc.cores=1) {
    stopifnot(all(c("Amplicon", "Target", "FW", "RV", "Sequenom") %in% names(x)))
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
	x$CGs <- mapply(function(cgs, offset) cgs+offset-1, gregexpr("CG", x$internal),
        start(x)+x$FW_len, SIMPLIFY=FALSE)

	#Add a GRangesList of primer positions
	x$primers <- GRangesList(lapply(1:length(x), function(i) {
	    tmp <- rep(x[i],2)
	    tmp <- resize(tmp, c(tmp$FW_len[1], tmp$RV_len[1]), fix=c("start", "end"))
	    values(tmp) <- NULL
	    tmp
	}))

	x
}

#' ampliconAnalysis
#'
#' Start to finish analysis of MiSeq bisulfite amplicon sequencing data
#'
#' @param amplicon_file Filename of the \code{csv} to read in amplicon data from
#' @param bams A named \code{character} vector of bams to be analysed for the supplied amplicons
#' @param genome \code{BSgenome} to map the amplicons to
#' @param paired \code{boolean} of whether to only load properly paired reads from the bams
#' @param minCov Minumum sequencing coverage required to report methylation/conversion estimates
#' @param min.map.q Minumum mapping quality for reads to be included in the analysis, set to 0 for no filtering
#' @param mc.cores Number of cores to use during processing
#' @return List
#'
#' @export
#'
#' @importFrom GenomicRanges GRangesList GRanges seqnames start end strand values countOverlaps resize subsetByOverlaps
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom IRanges IRanges endoapply elementLengths %over% cbind as.data.frame
#' @importFrom S4Vectors DataFrame
#' @importFrom Biostrings getSeq complement letterFrequency
#' @importFrom parallel mclapply
#' @importFrom GenomicAlignments readGAlignments pileLettersAt cigar
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
ampliconAnalysis <- function(amplicon_file, bams, genome, paired=TRUE, minCov=50, min.map.q=40, mc.cores=1) {
    # Quick checks
    stopifnot(file.exists(amplicon_file))
    stopifnot(all(file.exists(bams)))
    stopifnot(length(names(bams))==length(bams))
    
    # Read in, validate and align amplicon file
    amplicons <- read.csv(amplicon_file, stringsAsFactors=FALSE)
    stopifnot(all(c("Amplicon", "Target", "FW", "RV", "Sequenom") %in% names(amplicons)))
    message("Mapping amplicons to the supplied genome")
    amplicons <- ampliconGRanges(amplicons, genome, mc.cores)

    # Handy for s/lapplying
    amps <- 1:length(amplicons)
    names(amps) <- amplicons$Amplicon
    
    # Create a GRanges of all the bases in the amplicons we want to get a pileup at
    amplicon.bases <- sort(unique(unlist(GRangesList(lapply(1:length(amplicons), function(i) 
        GRanges(seqnames(amplicons)[i], 
            IRanges(start(amplicons)[i]:end(amplicons)[i], width=1), 
            strand=strand(amplicons)[i]))))))
    # Annotate with amplicon & base
    tmp <- sapply(1:length(amplicons), function(i) amplicon.bases %over% amplicons[i])
    amplicon.bases$amplicons <- apply(tmp, 1, function(i) paste(amplicons$Amplicon[i], collapse="+"))
    amplicon.bases$base <- getSeq(genome, amplicon.bases)
    # Remove primers from analysis
    amplicon.bases <- amplicon.bases[!amplicon.bases %over% unlist(amplicons$primers)]
    rm(tmp)
    
    # Read in libraries
    message("Reading in aligned sequencing libraries")
    if (paired) flags <- scanBamFlag(isPaired=TRUE, isProperPair=TRUE) else flags <- scanBamFlag()
    libs <- mclapply(bams, readGAlignments, param=ScanBamParam(flag=flags, what=c("seq", "mapq")), mc.cores=mc.cores)
    if (min.map.q>0) libs <- endoapply(libs, function(x) x[values(x)$mapq>=min.map.q])
    seqlevels(amplicons, force=TRUE) <- seqlevels(amplicon.bases, force=TRUE) <- seqlevels(libs[[1]]) # for pileup
    
    # Summary of how many reads hit amplicons
    amplicon_summary <- data.frame("No_Reads"=elementLengths(libs),
                            "On_Target_Reads"=sapply(libs, function(x) sum(unstrand(x) %over% amplicons)))
    amplicon_summary$On_Target_Percentage <- amplicon_summary$On_Target_Reads/amplicon_summary$No_Reads*100
    amplicon_counts <- sapply(libs, function(x) countOverlaps(amplicons, unstrand(x)))
    amplicon_counts <- if (length(amplicons)==1) as.matrix(amplicon_counts) else t(amplicon_counts)
    colnames(amplicon_counts) <- paste0(names(amps), "_Reads")
    amplicon_summary <- data.frame(amplicon_summary, amplicon_counts)
    rm(amplicon_counts)
    
    # Pileup of each library at each non-primer base
    message("Creating pileup at each base of each amplicon")
    nucl_piles <- GRangesList(mclapply(libs, function(x) {
        tmp <- amplicon.bases
        tmp.piles <- pileLettersAt(values(x)$seq, seqnames(x), start(x), cigar(x), amplicon.bases)
        minus <- as.character(strand(tmp))=="-"
        if (any(minus)) tmp.piles[minus] <- complement(tmp.piles[minus])
        values(tmp) <- cbind(values(tmp), DataFrame(letterFrequency(tmp.piles, c("A", "C", "G", "T"))))
        tmp$cov <- tmp$A+tmp$C+tmp$G+tmp$T
        tmp$A.rat <- tmp$A/tmp$cov
        tmp$C.rat <- tmp$C/tmp$cov
        tmp$G.rat <- tmp$G/tmp$cov
        tmp$T.rat <- tmp$T/tmp$cov
        tmp
    }, mc.cores=mc.cores))

    # Subset into C's only, change $base to dinucleotide
    message("Calculating methylation and conversion ratios")
    result <- amplicon.bases[amplicon.bases$base=="C"]
    result$base <- getSeq(genome, resize(result, 2, fix="start"), as.character=TRUE)
    tmp.C <- sapply(nucl_piles, function(x) x$C.rat[x$base=="C"])
    tmp.cov <- sapply(nucl_piles, function(x) x$cov[x$base=="C"])
    values(result) <- cbind(values(result), C=DataFrame(tmp.C), cov=DataFrame(tmp.cov))
    rm(tmp.C)
    rm(tmp.cov)
    
    # Subset to CpGs only
    result.CpGs <- result[result$base=="CG"]
    
    # Add to summary
    conversion_summary <- sapply(amps, function(i) {
        tmp <- subsetByOverlaps(result, amplicons[i])
        tmp <- as.data.frame(tmp[tmp$base!="CG"])
        tmp.C <- as.matrix(tmp[,grep("^C\\.", names(tmp))])
        tmp.cov <- as.matrix(tmp[,grep("^cov\\.", names(tmp))])
        # Remove low coverage data
        tmp.C[tmp.cov<minCov] <- NA
        # Conversion == 1-methylation
        1-colMeans(tmp.C, na.rm=TRUE)
    })
    colnames(conversion_summary) <- paste0(colnames(conversion_summary), "_Conversion")
    
    methylation_summary <- sapply(amps, function(i) {
        tmp <- subsetByOverlaps(result, amplicons[i])
        tmp <- as.data.frame(tmp[tmp$base=="CG"])
        tmp.C <- as.matrix(tmp[,grep("^C\\.", names(tmp))])
        tmp.cov <- as.matrix(tmp[,grep("^cov\\.", names(tmp))])
        # Remove low coverage data
        tmp.C[tmp.cov<minCov] <- NA
        colMeans(tmp.C, na.rm=TRUE)
    })
    colnames(methylation_summary) <- paste0(colnames(methylation_summary), "_Methylation")
    
    amplicon_summary <- data.frame(amplicon_summary, conversion_summary, methylation_summary)

    # Return all results
    list(summary=amplicon_summary, CpGs=result.CpGs, Cs=result, all_bases=nucl_piles, amplicons=amplicons)
}

#' regionalAmpliconAnalysis
#'
#' Condenses the analysis of individual amplicons to those grouped in regions
#'
#' @param x Result returned by \link{ampliconAnalysis}
#' @param grouping \code{character} vector of region names by which amplicons are grouped together
#' @param minCov Minumum sequencing coverage required to report methylation/conversion estimates
#' @return List
#'
#' @export
#'
#' @importFrom GenomicRanges GRangesList seqnames countOverlaps split findOverlaps values
#' @importFrom IRanges subjectHits
#' @importMethodsFrom GenomicRanges range
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
regionalAmpliconAnalysis <- function(x, grouping, minCov=50) {
    stopifnot(length(grouping)==length(x$amplicons))

    # Create GRanges of boundaries of regions
    regions <- unlist(GRangesList(lapply(split(x$amplicons, grouping), function(x) {
        stopifnot(length(unique(seqnames(x)))==1)
        range(x, ignore.strand=TRUE)
    })))
    stopifnot(all(countOverlaps(regions)==1))

    # Annotate the regions with which amplicons they comprise of
    x.ov <- as.matrix(findOverlaps(regions, x$amplicons))
    stopifnot(!any(duplicated(x.ov[,2])))
    x.names <- tapply(x$amplicons$Amplicon[x.ov[,2]], x.ov[,1], paste, collapse="+")
    regions$Amplicons[as.integer(names(x.names))] <- unname(x.names)


    # Reannotate Cs & CpGs as regions not amplicons
    stopifnot(names(values(x$Cs))[[1]]=="amplicons", names(values(x$CpGs))[[1]]=="amplicons")
    regions.Cs <- x$Cs
    names(values(regions.Cs))[1] <- "region"
    regions.Cs$region <- names(regions)[subjectHits(findOverlaps(regions.Cs, regions))]
    regions.CpGs <- x$CpGs
    names(values(regions.CpGs))[1] <- "region"
    regions.CpGs$region <- names(regions)[subjectHits(findOverlaps(regions.CpGs, regions))]

    # Summarise per region
    regions_summary <- lapply(split(regions.Cs, regions.Cs$region), function(y) {
        y <- as.data.frame(values(y))
        # Remove methylation data below minCov
        toNA <- y[,grep("^cov\\.", names(y))] < minCov
        y[,grep("^C\\.", names(y))][toNA] <- NA
        CG <- y$base=="CG"
        res <- cbind("coverage"=colMeans(y[,grep("^cov\\.", names(y))], na.rm=TRUE),
                  "methylation"=colMeans(y[CG,grep("^C\\.", names(y))], na.rm=TRUE),
                   "conversion"=1-colMeans(y[!CG,grep("^C\\.", names(y))], na.rm=TRUE))
        rownames(res) <- gsub("^cov.", "", rownames(res))
        res
    })
    c(x, list("regions"=regions, "regions_summary"=regions_summary, "regions_Cs"=regions.Cs, "regions_CpGs"=regions.CpGs))
}
