

# Function Name: 	disambiguateMultihits
# Description: 		among multiple alignments of the same read (i.e. multihit),
#					select the alignment corresponding to the bin with maximum
#					posterior for the enriched hidden state.
# Input: 			GAlignments with 1-read-(>=1)-alignments
#					posterior decoding from HMM (GRanges) 
# Output:			the same GAlignments with 1-read-1-alignment 
#
# Author: Yue Li
###############################################################################

disambiguateMultihits2 <- function(alignGal, 
                                   nbhGRList, 
                                   postprobCutoff=0){
  require(dplyr)
  # unlist GRangesList to search multihits over all chromosomes
  nbhGR <- unlist(nbhGRList)
  names(nbhGR) <- NULL
  
  # find indices of multihits labeled as FALSE for column "uniqueHit"
  # labelling is expected to have been performed by getAlignGal
  idx <- which(values(alignGal)$uniqueHit == FALSE)
  
  # subset alignGal multihit reads
  mhitReads <- alignGal[idx]
  
  # for each multihit, get indices of the corresponding genomic regions from nbhGR
  # mhitAllOverlaps <- findOverlaps(mhitReads, nbhGR)
  tmp = as.data.frame(nbhGR)
  sq = c(seq(1, nrow(tmp), 100000),nrow(tmp))
  tmplist <- lapply(2:length(sq), function(i){
    a = sq[i-1]
    b = sq[i]-1
    cat(a, b, "\n")
    sub = tmp[a:b,]
    sub <- as(sub, "GRanges")
    sublist <- GRangesList(sub)
    
    ### disambiguateMultihits
    nbhGRList <- sublist
    nbhGR <- unlist(nbhGRList)
    names(nbhGR) <- NULL
    idx <- which(values(alignGal)$uniqueHit == FALSE)
    mhitReads <- alignGal[idx]
    over = findOverlaps(mhitReads, nbhGR)
    over@to = as.integer(over@to + (sq[i-1]-sq[1]))
    over
  })
  x = tmplist[[1]]
  x@from = lapply(tmplist, function(x) x@from) %>% unlist
  x@to = lapply(tmplist, function(x) x@to) %>% unlist
  x@nRnode = nrow(tmp)
  mhitAllOverlaps = x
  
  # index multireads in the order as nbhGR interesected regions occur
  mhitQuery <- mhitReads[queryHits(mhitAllOverlaps)]
  
  # then assign enriched state posterior to each alignment				
  values(mhitQuery) <- values(nbhGR[subjectHits(mhitAllOverlaps)])$state_2_postprob
  names(values(mhitQuery)) <- "state_2_postprob"
  mhitMatrix <- cbind(names(mhitQuery), 1:length(mhitQuery), 							
                      values(nbhGR[subjectHits(mhitAllOverlaps)])$state_2_postprob)
  mhitList <- split(mhitMatrix, mhitMatrix[,1])
  
  
  # central step: split alignment by read names into a list
  # for each list element (containing > 1 alignments for the same read)
  # select the alignment corresponding to the region with max postprob
  # return only the unique index of that alignment in order to unlist the results
  desiredIdx <- lapply(mhitList,
                       function(x) {
                         y <- matrix(x, ncol=3)
                         maxIdx <- which.max(y[,3])
                         if(max(y[maxIdx,3]) > postprobCutoff) as.numeric(y[maxIdx, 2])
                       }
  )
  
  desiredIdx <- unlist(desiredIdx)
  selectGal <- mhitQuery[desiredIdx, ]
  values(selectGal) <- FALSE
  names(values(selectGal)) <- "uniqueHit"
  alignGalfiltered <- c(alignGal[values(alignGal)$uniqueHit == TRUE], selectGal)
  message(sprintf("\n%d/%d multihit reads corresponding to %d ambiguous alignments\nhave been assigned to %d unique regions with maximum posterior for the enriched state\n", 
                  length(selectGal), length(alignGalfiltered), 
                  length(idx), length(selectGal)))
  
  return(alignGalfiltered)			
}

annotateRIP2 <- function(sigGRanges, 
                         biomaRt_dataset, 
                         featureType="TSS", 
                         goAnno,
                         strandSpecific=FALSE, 
                         exportFormat = "txt",
                         hasGOdb=!missing(goAnno), 
                         goPval=0.1, 
                         outDir, 
                         ...){
  
  stopifnot(!missing(sigGRanges))
  stopifnot(!missing(biomaRt_dataset))
  stopifnot(require(biomaRt))
  stopifnot(require(ChIPpeakAnno))
  
  ##### hack useMart to ignore unused arguments
  formals(useMart) <- c(formals(useMart), alist(... = ))	
  # mart <- useMart(dataset=biomaRt_dataset, ...)	
  # anno <- getAnnotation(mart, featureType=featureType)
  message(sprintf("\n\n*** Annotating %d genomic ranges with %s.\n",
                  length(sigGRanges), biomaRt_dataset))
  
  names(sigGRanges) <- NULL
  
  ###### hack annotatePeakInBatch to ignore unused arguments
  formals(annotatePeakInBatch) <- c(formals(annotatePeakInBatch), alist(... = ))
  annotatedPeak <- annotatePeakInBatch(GRanges(sigGRanges), AnnotationData = anno, output="both", ...)
  if(strandSpecific) subsetByOverlaps(annotatedPeak, sigGRanges, ignore.strand=FALSE)
  
  ##### more useful information based on ensembl gene ID
  geneInfo <- getBM(mart=mart, attributes=c("ensembl_gene_id", "external_gene_name", "description"),
                    filters="ensembl_gene_id", values = annotatedPeak$feature)
  geneInfo  <- geneInfo[match(annotatedPeak$feature, geneInfo$ensembl_gene_id),]
  sigGRangesIdx <- as.numeric(annotatedPeak$peak)
  sigGRangesAnnotated <- sigGRanges[sigGRangesIdx]
  combinedInfo <- cbind(geneInfo, as.data.frame(values(sigGRangesAnnotated)),
                        as.data.frame(values(annotatedPeak)) )
  
  ##### remove useless information: "space" (chr)
  combinedInfo <- combinedInfo[, grep("space", colnames(combinedInfo), invert=TRUE)]
  colnames(combinedInfo)[colnames(combinedInfo) == "strand"] <- "feature_strand"
  values(sigGRangesAnnotated) <- combinedInfo
  sigGRangesAnnotated <- sigGRangesAnnotated[order(values(sigGRangesAnnotated)$peak)]
  if(!missing(goAnno)) {
    hasGOdb <- require(package=goAnno, character.only=TRUE)
    if(hasGOdb) {
  
      message(sprintf("\n\n*** GO analysis for %d associated features with %s.\n",
                      length(annotatedPeak$feature), goAnno))
      
      # hack functions to ignore unused arguments
      formals(getEnrichedGO) <- c(formals(getEnrichedGO), alist(... = ))
      enrichedGO <- getEnrichedGO(annotatedPeak, orgAnn=goAnno, maxP = goPval,
                                  multiAdj = TRUE, multiAdjMethod="BH", ...)
      # remove redundant sets
      enrichedGO$bp <- unique( enrichedGO$bp[,-11] )
      enrichedGO$mf <- unique( enrichedGO$mf[,-11] )
      enrichedGO$cc <- unique( enrichedGO$cc[,-11] )
      enrichedGO <- rbind(enrichedGO$bp[order(enrichedGO$bp$pvalue), ],
                          enrichedGO$mf[order(enrichedGO$mf$pvalue), ],
                          enrichedGO$cc[order(enrichedGO$cc$pvalue), ])
    } else {
      warning(sprintf("%s is not found!", goAnno))
    }
  }
  

  ################ save and export results to outDir ################		
  if(!missing(outDir)) {
    
    # remove backslash to avoid double backslash in the following path names
    outDir <- sub("/$", "", outDir)
    outfile <- paste(outDir, "/RIPGRanges_annotated.RData", sep="")
    message(sprintf("\n\n*** Saving RData to %s\n", outfile))
 
    ################ save results in RData ################
    
    save(sigGRangesAnnotated, file=outfile)

    ################ export RIP regions ################
    outfile <- paste(outDir, "/RIPregions_annotated.", exportFormat, sep="")
    message(sprintf("\n\n*** Exporting %s\n", outfile))
    exportGRanges(gRanges=sigGRangesAnnotated, outfile=outfile, 
                  exportFormat=exportFormat)

    ################ export enriched GO ################
    if(hasGOdb) {
      outfile <- paste(outDir, "/RIPregions_enrichedGO.txt", sep="")
      
      message(sprintf("\n\n*** Exporting %s\n", outfile))
      
      write.table(enrichedGO, file=outfile, row.names=F, quote=F, sep="\t")						
    }
  }
  
  if(hasGOdb) return(list(sigGRangesAnnotated=sigGRangesAnnotated, enrichedGO=enrichedGO))
  if(!hasGOdb) return(sigGRangesAnnotated)
  
}

