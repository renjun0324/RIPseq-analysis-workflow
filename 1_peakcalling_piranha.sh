
##################################################################
#                                                                
#               peak calling METHOD 1 - Piranha 
#                                                                 
##################################################################

### 1. 将bam文件转为bed文件
bam_paths=("/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/3_sequence_bowtie2/HuR"  \
           "/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/3_sequence_bowtie2/input")
for p in ${bam_paths[*]};do
echo $p
cd $p
bamfile=`find . -type f |grep .redup.bam$`
# bamfile=`find . -name *.out.bam`
outbedfile=$(basename -s .bam $bamfile).bed
echo $bamfile
echo $outbedfile
bamToBed -i $bamfile>$outbedfile
done      

### 2. 保留chr1-chr22, chrM, X, Y, 生成新的bed文件
# cat IgG-IPAligned.sortedByCoord.out.redup.bed | cut -f 1 | uniq -c
# cat inputAligned.sortedByCoord.out.redup.chr.bed | cut -f 1 | uniq -c
for p in ${bam_paths[*]};do
echo $p
cd $p
bed_file=`find . -name *.redup.bed`
# bed_file=`find . -name *.out.bed`
outbed_file=$(basename -s .bed $bed_file).chr.bed
cat $bed_file | grep chr > $outbed_file
done   

### 3. Piranha
cd 4_piranh
control="../3_sequence_bowtie2/input/input.sorted.redup.chr.bed"
treat="../3_sequence_bowtie2/HuR/HuR.sorted.redup.chr.bed"

Piranha $treat $control -s -z 100 -p 0.05 -o result1.bed
# Piranha $treat $control -s -b 100 -a 100 -p 0.05 -o result2.bed

### 4. debug
# 提示协变量必须要和treat的基因位点对应，确实不对应
# 只能用下面的方法操作然后注释
# Piranha $treat -s -z 1 -p 0.05 -o HA-IP.bed

#################################################################
#                                                                
#                            ChIPseeker                     
#                                                                 
##################################################################

setwd("4_piranha/");peaks <- readPeakFile("result1.bed")

## 注释
library(ChIPseeker)
library(GenomicFeatures)
txdb <- makeTxDbFromGFF(file = "/share/data6/tmp/renjun/Reference/refdata-gex-mm10-2020-A/genes/genes.gtf",
                        format = "gtf", organism = "Mus musculus")
peakAnno <- annotatePeak(peaks, 
                         tssRegion = c(-5000, 5000), 
                         TxDb = txdb, 
                         addFlankGeneInfo = TRUE, 
                         flankDistance = 5000,
                         annoDb = "org.Mm.eg.db")
df <- as.data.frame(as.GRanges(peakAnno))
colnames(df)[6:10] = c("x","peak_score","peak_strand", "pval", "p.adjust")
write.csv(df, file = paste0("chipseeker_annotatedPeak.csv"))
save(df, file = paste0("chipseeker_annotatedPeak.rda"))    

## 过滤
pir = df
pir %>% 
  filter(p.adjust<0.01, annotation=="Promoter (<=1kb)") %>% 
  dplyr::select(SYMBOL) %>% 
  unlist %>% 
  unique