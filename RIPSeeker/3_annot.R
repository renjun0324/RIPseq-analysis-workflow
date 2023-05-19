
#-----------------------------------------------------------------
#                                                                
#                       parameter setting
#                                                                 
#-----------------------------------------------------------------

biomaRt_dataset = "hsapiens_gene_ensembl"
featureType = "TSS"
goAnno = "org.Hs.eg.db"
strandSpecific = FALSE
exportFormat = "csv"
goPval = 1

#-----------------------------------------------------------------
#                                                                
#                     annotation - ripseeker
#                                                                 
#-----------------------------------------------------------------

### 导入reference
# mart = useMart(biomart = "ensembl", dataset = biomaRt_dataset)
# anno = getAnnotation(mart, featureType=featureType)
load("/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/7_RIPSeeker/ref/anno.rda")
load("/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/7_RIPSeeker/ref/mart.rda")

### annotation
rip_anno = annotateRIP2(sigGRanges = ripGR,
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

#-----------------------------------------------------------------
#                                                                
#             annotation - ChIPpeakAnno (=annotateRIP2)
#                                                                 
#-----------------------------------------------------------------

### 1. 导入reference
# mart = useMart(biomart = "ensembl", dataset = biomaRt_dataset)
# anno = getAnnotation(mart, featureType=featureType)
load("/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/7_RIPSeeker/ref/anno.rda")
load("/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/7_RIPSeeker/ref/mart.rda")

### 2. 注释
# annotatedPeak <- annotatePeakInBatch(GRanges(sigGRanges), AnnotationData = anno, output="both", ...)
peakAnno <- annotatePeakInBatch(ripGR, AnnotationData = anno, output="both", bindingRegion = c(-5000, 5000))

### 3. id转换
geneInfo <- getBM(mart=mart,
                  attributes=c("ensembl_gene_id", "external_gene_name", "description"),
                  filters="ensembl_gene_id", 
                  values = peakAnno$feature)
geneInfo  <- geneInfo[match(peakAnno$feature, geneInfo$ensembl_gene_id),]

### 4. 输出
df = data.frame(as.data.frame(peakAnno), geneInfo)
write.csv(df, file = "chippeakanno_annotatedPeak.csv")

#-----------------------------------------------------------------
#                                                                
#                     annotation - chipseeker
#                                                                 
#-----------------------------------------------------------------

### 1. 准备reference
txdb <- makeTxDbFromGFF(file = "/share/data6/tmp/renjun/Reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
                        format = "gtf", organism = "Homo sapiens")

### 2. 注释peaks
peakAnno <- annotatePeak(ripGR, 
                         tssRegion=c(-5000, 5000), 
                         TxDb=txdb, 
                         addFlankGeneInfo=TRUE, 
                         flankDistance=5000,
                         annoDb = "org.Hs.eg.db")

### 3. 输出
df <- as.data.frame(as.GRanges(peakAnno))
# df$SYMBOL[which(df$pvalAdj<0.05)] %>% unique
write.csv(df, file = paste0("chipseeker_annotatedPeak.csv"))
save(df, file = paste0("chipseeker_annotatedPeak.rda"))