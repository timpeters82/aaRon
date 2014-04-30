annotateRegions <- function(reg) {

    # Want to do most analysis for "just" protein coding as well as all genes
    tx2 <- tx[tx$tx_type=="protein_coding"]

    # Number of different genes TSSs within the DMR
    reg.TSS <- as.matrix(findOverlaps(reg, resize(tx, 1, fix="start")))
    reg.TSS <- tapply(tx$gene_name[reg.TSS[,2]], reg.TSS[,1], function(x) length(unique(x)))
    reg$nGeneTSS <- 0
    reg$nGeneTSS[as.integer(names(reg.TSS))] <- unname(reg.TSS)
    #just for protein coding genes
    reg.TSS <- as.matrix(findOverlaps(reg, resize(tx2, 1, fix="start")))
    reg.TSS <- tapply(tx2$gene_name[reg.TSS[,2]], reg.TSS[,1], function(x) length(unique(x)))
    reg$nProtGeneTSS <- 0
    reg$nProtGeneTSS[as.integer(names(reg.TSS))] <- unname(reg.TSS)
    

    # Distance to closest TSS
    reg.dist <- as.data.frame(distanceToNearest(reg, resize(tx, 1, fix="start")))
    reg$distanceTSS <- reg.dist$distance
    reg$TSS <- reg.dist$subjectHits
    reg$tx_name <- tx$tx_name[reg$TSS]
    reg$tx_type <- tx$tx_type[reg$TSS]
    reg$gene_id <- tx$gene_name[reg$TSS]
    reg$gene_name <- tx$symbol[reg$TSS]
    reg$result <- gx[tx$gene_name[reg$TSS]]$result

    # Distance to closest protein-coding TSS
    reg.dist <- as.data.frame(distanceToNearest(reg, resize(tx2, 1, fix="start")))
    reg$distanceTSS_prot <- reg.dist$distance
    reg$TSS_prot <- reg.dist$subjectHits
    reg$tx_name_prot <- tx2$tx_name[reg$TSS_prot]
    reg$gene_id_prot <- tx2$gene_name[reg$TSS_prot]
    reg$gene_name_prot <- tx2$symbol[reg$TSS_prot]
    reg$result_prot <- gx[tx2$gene_name[reg$TSS_prot]]$result

    # Distance to closest CpG island
    reg.dist <- as.data.frame(distanceToNearest(reg, CpGislands))
    reg$distanceCpGi <- reg.dist$distance
    reg$CpGi <- reg.dist$subjectHits

    # % promoter/genebody/intergenic (2kb up and down)
    promotersGR <- reduce(strip(resize(resize(tx, 1, fix="start"), 4000, fix="center")))
    genebodyGR <- setdiff(reduce(strip(tx)), promotersGR)
    genomeGR <- GRanges(seqlevels(tx), IRanges(1, seqlengths(tx)))
    intergenicGR <- setdiff(setdiff(genomeGR, genebodyGR), promotersGR)

    reg$promoter <- coverageRatio(reg, promotersGR)
    reg$genebody <- coverageRatio(reg, genebodyGR)
    reg$intergenic <- coverageRatio(reg, intergenicGR)
    
    # % CpG island/CpG shore/nonCpG
    reg$CpGisland <- coverageRatio(reg, CpGislands)
    reg$CpGshores <- coverageRatio(reg, CpGshores)
    reg$nonCpG <- coverageRatio(reg, setdiff(setdiff(genomeGR, CpGislands), CpGshores))

    reg
}

