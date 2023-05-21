

##################################################################
#                                                                
#            分割染色体之后按照提交任务的方式运行
#                                                                 
##################################################################

## 1. 按照染色体分割bam文件     
bam_paths=("/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/3_sequence_bowtie2/HuR"  \
           "/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/3_sequence_bowtie2/input")
for p in ${bam_paths[*]};do
mkdir $p"split"
cd $p"split"
echo $($p"split")
bam_file=`find $p -name "*.redup.bam"`
bamtools split -in $bam_file -reference
mv $p*_*.bam .
done

## 2. mainSeek: running 0_run.sh
bash /share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/0_run.sh input 3 20 20

##################################################################
#                                                                
#               peak calling METHOD3 - Ripseeker                
#                                                                 
##################################################################

library(RIPSeeker)
bamFiles=c("/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/3_sequence_bowtie2/HuR/HuR.sorted.redup.bam",
           "/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/3_sequence_bowtie2/input/input.sorted.redup.bam")

## 1. mainseek
## 优化出一个mainsize，后面计算seekRIP的时候会因为binsize不同而到时片段数不一样无法计算差异
## 所以如果有对照组的情况下，最好手动设置binsize
mainseek_hur = mainSeek(bamFiles[1],
                        binSize=200,
                        reverseComplement = FALSE,
                        genomeBuild = "mm10",
                        uniqueHit = TRUE,
                        assignMultihits = TRUE,
                        strandType = NULL,
                        paired=TRUE,
                        rerunWithDisambiguatedMultihits = TRUE,
                        silentMain = FALSE,
                        multicore = TRUE,
                        returnAllResults = TRUE)
save(mainseek_hur, file = "mainseek_hur.rda")

mainseek_input = mainSeek(bamFiles[2],
                          binSize=200,
                          reverseComplement = FALSE,
                          genomeBuild = "mm10",
                          uniqueHit = TRUE,
                          assignMultihits = TRUE,
                          strandType = NULL,
                          paired=TRUE,
                          rerunWithDisambiguatedMultihits = TRUE,
                          silentMain = FALSE,
                          multicore = TRUE,
                          returnAllResults = TRUE)
save(mainseek_input, file = "mainseek_input.rda")

## 3. ripseek
load("mainseek_input.rda")
load("mainseek_hur.rda")
seekrip_list = list()
for(i in 1:length(mainseek_hur$nbhGRList)){
  cat(i, "\n")
  x = names(mainseek_hur$nbhGRList)[i]
  # if(as.character(mainseek_hur$nbhGRList[[x]])[1]=="none" | as.character(mainseek_input$nbhGRList[[x]])[1]=="none") next
  if(nrow(data.frame(mainseek_hur$nbhGRList[[x]])) != nrow(data.frame(mainseek_input$nbhGRList[[x]]))) next
  seekrip = seekRIP(nbhGRRIP = mainseek_hur$nbhGRList[[x]],
                    nbhGRCTL = mainseek_input$nbhGRList[[x]])
  seekrip_list = append(seekrip_list, list(seekrip))
}
rip <- map_dfr(seekrip_list, function(x) data.frame(x))
save(seekrip_list, file = "seekrip_list.rda")
save(rip, file = "rip.rda")

##################################################################
#                                                                
#                     annotation - ripseeker
#                                                                 
##################################################################

biomaRt_dataset = "mmusculus_gene_ensembl"
featureType = "TSS"
goAnno = "org.Mm.eg.db"
strandSpecific = FALSE
exportFormat = "csv"
goPval = 1

### 导入reference
# mart = biomaRt::useMart(biomart = "ensembl", dataset = biomaRt_dataset)
# anno = ChIPpeakAnno::getAnnotation(mart, featureType=featureType)
# load("/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/7_RIPSeeker/ref/anno.rda")
# load("/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/7_RIPSeeker/ref/mart.rda")
load("/share/data6/tmp/renjun/Reference/ripseeker_ref/org.Mm.eg.db_tss_mmusculus_gene_ensembl/anno.rda")
load("/share/data6/tmp/renjun/Reference/ripseeker_ref/org.Mm.eg.db_tss_mmusculus_gene_ensembl/mart.rda")

### annotation
x <- as(rip, "GRanges")
rip_anno = annotateRIP2(sigGRanges = x,
                        biomaRt_dataset = biomaRt_dataset,
                        featureType = featureType,
                        goAnno = NULL,
                        strandSpecific = FALSE,
                        exportFormat = exportFormat,
                        hasGOdb=!missing(goAnno),
                        goPval = goPval,
                        outDir = NULL)
rip_anno = as.data.frame(rip_anno)
write.csv(rip_anno, file = "RIPregions_annotated.csv")

##################################################################
#                                                                
#                   annotation - ChIPpeakAnno
#                                                                 
##################################################################

library(ChIPseeker)
library(GenomicFeatures)
txdb <- makeTxDbFromGFF(file = "/share/data6/tmp/renjun/Reference/refdata-gex-mm10-2020-A/genes/genes.gtf",
                        format = "gtf", organism = "Mus musculus")
load("rip.rda")
peaks = as(rip, "GRanges")
peakAnno <- annotatePeak(peaks, 
                         tssRegion = c(-5000, 5000), 
                         TxDb = txdb, 
                         addFlankGeneInfo = TRUE, 
                         flankDistance = 5000,
                         annoDb = "org.Mm.eg.db")
df <- as.data.frame(as.GRanges(peakAnno))

write.csv(df, file = paste0("chipseeker_annotatedPeak.csv"))
save(df, file = paste0("chipseeker_annotatedPeak.rda"))

## 过滤
rip = df
rip %>% 
  filter(eFDR<0.05, annotation=="Promoter (<=1kb)") %>% 
  dplyr::select(SYMBOL) %>% 
  unlist %>% 
  unique

rip %>% 
  filter(pvalAdj<0.05, annotation=="Promoter (<=1kb)") %>% 
  dplyr::select(SYMBOL) %>% 
  unlist %>% 
  unique


##################################################################
#                                                                
#                   annotation - ChIPpeakAnno
#                                                                 
##################################################################

### 1. 导入reference
# mart = useMart(biomart = "ensembl", dataset = biomaRt_dataset)
# anno = getAnnotation(mart, featureType=featureType)
# load("/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/7_RIPSeeker/ref/anno.rda")
# load("/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/7_RIPSeeker/ref/mart.rda")

### 2. 注释
# annotatedPeak <- annotatePeakInBatch(GRanges(sigGRanges), AnnotationData = anno, output="both", ...)
peakAnno <- annotatePeakInBatch(x, AnnotationData = anno, output="both", bindingRegion = c(-5000, 5000))

### 3. id转换
geneInfo <- getBM(mart=mart,
                  attributes=c("ensembl_gene_id", "external_gene_name", "description"),
                  filters="ensembl_gene_id", 
                  values = peakAnno$feature)
geneInfo  <- geneInfo[match(peakAnno$feature, geneInfo$ensembl_gene_id),]

### 4. 输出
df = data.frame(as.data.frame(peakAnno), geneInfo)
write.csv(df, file = "chippeakanno_annotatedPeak.csv")

##################################################################
#                                                                
#                     mainSeekSingleChrom
#                                                                 
##################################################################

# ## 2. single mainseek chrom - HuR
# alignGal <- getAlignGal(alignFilePath = bamFiles[1], format = "BAM", genomeBuild = "mm10", paired = TRUE)
# alignGR <- as(alignGal, "GRanges")
# alignGRList <- GRangesList(as.list(split(alignGR, seqnames(alignGR))))
# mainseek_list = list()
# for(i in 1:length(alignGRList)){
#   cat(i, "\n")
#   if(nrow(data.frame(alignGRList[[i]]))<10){
#     mainseek_list = append(mainseek_list, list("none"))
#   }else{
#     nbhGR <- mainSeekSingleChrom(alignGR=alignGRList[[i]],  binSize= 200)
#     mainseek_list = append(mainseek_list, list(nbhGR))
#   }
# }
# mainseek_hur = mainseek_list
# names(mainseek_hur) = names(alignGRList)
# save(mainseek_hur, file = "mainseek_hur.rda")
# 
# ## 2. single mainseek chrom - input
# alignGal <- getAlignGal(alignFilePath = bamFiles[2], format = "BAM", genomeBuild = "mm10", paired = TRUE)
# alignGR <- as(alignGal, "GRanges")
# alignGRList <- GRangesList(as.list(split(alignGR, seqnames(alignGR))))
# mainseek_list = list()
# for(i in 1:length(alignGRList)){
#   cat(i, "\n")
#   if(nrow(data.frame(alignGRList[[i]]))<10){
#     mainseek_list = append(mainseek_list, list("none"))
#   }else{
#     nbhGR <- mainSeekSingleChrom(alignGR=alignGRList[[i]],  binSize= 200)
#     mainseek_list = append(mainseek_list, list(nbhGR))
#   }
# }
# mainseek_input = mainseek_list
# names(mainseek_input) = names(alignGRList)
# save(mainseek_input, file = "mainseek_input.rda")

# ## 3. ripseek
# load("mainseek_input.rda")
# load("mainseek_hur.rda")
# seekrip_list = list()
# for(i in 1:length(alignGRList)){
#   cat(i, "\n")
#   x = names(mainseek_hur)[i]
#   if(as.character(mainseek_hur[[x]])[1]=="none" | as.character(mainseek_input[[x]])[1]=="none") next
#   if(nrow(data.frame(mainseek_hur[[x]])) != nrow(data.frame(mainseek_input[[x]]))) next
#   seekrip = seekRIP(nbhGRRIP = mainseek_hur[[x]],
#                     nbhGRCTL = mainseek_input[[x]])
#   seekrip_list = append(seekrip_list, list(seekrip))
# }
# rip <- map_dfr(seekrip_list, function(x) data.frame(x))    
# save(seekrip_list, rip, file = "seekrip_list.rda")
