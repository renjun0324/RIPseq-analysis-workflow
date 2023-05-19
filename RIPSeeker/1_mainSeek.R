
args<-commandArgs(T) 
i <- as.numeric(args[1])
choose <- as.character(args[2])
binsize <- as.numeric(args[3])
root_path <- as.character(args[4])

#-----------------------------------------------------------------
#                                                                
#                         路径及文件读取准备
#                                                                 
#-----------------------------------------------------------------

library(RIPSeeker)
library(ChIPpeakAnno)
library(ChIPseeker)
library(clusterProfiler)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(biomaRt)
library(dplyr)

source(paste0(root_path,"/7_RIPSeeker/0_RIPSeeker.R"))

## 1. 建立group-binsize-chr文件夹
chr = paste0("chr", c(1:22,"M","X","Y"))
save_path <- paste0(root_path, "/7_RIPSeeker/", choose, "/binsize", binsize, "/", chr[i])
system(paste0("mkdir -p ", save_path))

## 2. 获取原始bam文件的完整路径
bam_path = paste0(root_path, "/4_bam/", choose, "/split")
base_name = system(paste0("tmp=",list.files(path=bam_path)[1]," && echo ${tmp%%_*}"), intern = TRUE)
bam_path = paste0(bam_path, "/", base_name, "_", chr[i], ".bam")

# ## 非染色体部分
# ip_other <- setdiff(list.files(paste0(bampath,"/ip/split/")),
#                     paste0("inputAligned.sortedByCoord.out.REF_",chr,".bam"))
# input_other <- setdiff(list.files(paste0(bampath,"/input/split/")),
#                        paste0("inputAligned.sortedByCoord.out.REF_",chr,".bam"))
# ip_other_n <- stringr::str_split_fixed(ip_other,"_",2)[,2]
# input_other_n <- stringr::str_split_fixed(input_other,"_",2)[,2]
# intersect_other <- intersect(input_other_n,ip_other_n)
# intersect_other <- stringr::str_split_fixed(intersect_other,"\\.bam",2)[,1]
# ip_otherFiles <- paste0(bampath,"/ip/split/inputAligned.sortedByCoord.out.REF_",intersect_other,".bam")
# input_otherFiles <- paste0(bampath,"/input/split/inputAligned.sortedByCoord.out.REF_",intersect_other,".bam")
# 
# ## 所有的文件
# chr <- c(chr, intersect_other)
# ipFiles <- c(ipFiles,ip_otherFiles)
# inputFiles <- c(inputFiles,input_otherFiles)

#-----------------------------------------------------------------
#                                                                
#                            parameters
#                                                                 
#-----------------------------------------------------------------

## getAlignGal
paired = TRUE
format = "bam"
genomeBuild = "x"
deleteGeneratedBAM = FALSE
reverseComplement=FALSE
returnDuplicate = FALSE
flagMultiHits = TRUE
returnAllResults = TRUE
returnOnlyUniqueHits=FALSE
rerunWithDisambiguatedMultihits = TRUE

## mainSeekSingleChrom
K = 2
binSize = binsize
minBinSize = 200
maxBinSize = 1200
minReadCount = 10
backupNumBins = 10
increment = 5
pathToSavePlotsOfBinSizesVersusCosts = save_path
verbose = TRUE
allowSecondAttempt = TRUE

#-----------------------------------------------------------------
#                                                                
#                             running
#                                                                 
#-----------------------------------------------------------------

alignGal <- getAlignGal(bam_path[1],
                        format = format, 
                        genomeBuild = genomeBuild, 
                        deleteGeneratedBAM = deleteGeneratedBAM, 
                        reverseComplement = reverseComplement, 
                        returnDuplicate = returnDuplicate, 
                        flagMultiHits = flagMultiHits, 
                        returnOnlyUniqueHits = returnOnlyUniqueHits, 
                        paired = paired)
alignGR <- as(alignGal, "GRanges")
alignGR <- addPseudoAlignment(alignGR)
alignGRList <- GRangesList(as.list(split(alignGR, seqnames(alignGR))))
# runViterbi <- all(values(alignGal)$uniqueHit) || !rerunWithDisambiguatedMultihits
mainseek = mainSeekSingleChrom(alignGRList[[chr[i]]],
                               K = K,
                               runViterbi = TRUE,
                               # binSize = median(width(get(ip)$nbhGRList)),
                               binSize = binSize,
                               minReadCount = minReadCount,
                               backupNumBins = backupNumBins,
                               minBinSize = minBinSize,
                               maxBinSize = maxBinSize,
                               increment = increment,
                               pathToSavePlotsOfBinSizesVersusCosts = pathToSavePlotsOfBinSizesVersusCosts,
                               verbose = verbose,
                               allowSecondAttempt = allowSecondAttempt)
alignGRList <- GRangesList(mainseek)
alignGalFiltered <- disambiguateMultihits2(alignGal, alignGRList)
alignGR <- as(alignGalFiltered, "GRanges")
alignGR <- addPseudoAlignment(alignGR)
alignGRList <- GRangesList(as.list(split(alignGR, seqnames(alignGR))))
mainseek2 = mainSeekSingleChrom(alignGRList[[chr[i]]],
                                K = K,
                                runViterbi = TRUE,
                                # binSize = median(width(get(ip)$nbhGRList)),
                                binSize = binSize,
                                minReadCount = minReadCount,
                                backupNumBins = backupNumBins,
                                minBinSize = minBinSize,
                                maxBinSize = maxBinSize,
                                increment = increment,
                                pathToSavePlotsOfBinSizesVersusCosts = pathToSavePlotsOfBinSizesVersusCosts,
                                verbose = verbose,
                                allowSecondAttempt = allowSecondAttempt)
x = list(nbhGRList=mainseek2, alignGal=alignGal, alignGalFiltered=alignGalFiltered)
  
#-----------------------------------------------------------------
#                                                                
#                              save
#                                                                 
#-----------------------------------------------------------------

name <- paste0(choose,"_", chr[i])
assign(name, x)
save(list = name, file = paste0(save_path, "/", name,".rda"))
