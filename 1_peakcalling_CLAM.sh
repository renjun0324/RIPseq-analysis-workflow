##################################################################
#                                                                
#                peak calling METHOD 2 - CLAM                      
#                                                                 
##################################################################

bam_paths=("/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/3_sequence_bowtie2/HuR"  \
           "/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/3_sequence_bowtie2/input")

### 1. 从bam文件中提取部分染色体的数据并建立索引
for p in ${bam_paths[*]};do
echo $p
cd $p
bamfile=`find . -type f |grep .redup.bam$`
outbamfile=$(basename $p).chr.bam
echo $bamfile
echo $outbamfile
samtools view -b $bamfile \
chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
chr21 chr22 chrM chrX chrY > $outbamfile

samtools index $outbamfile
done

### 2. realigner
cd 6_CLAM
for p in ${bam_paths[*]};do
dir_name=`basename $p`
mkdir $dir_name
cd $dir_name 

# bam_file=`find $p -name "*.redup.bam"`
bam_file=`find $p -name "*.chr.bam"`
echo $bam_file
# CLAM preprocessor -i $bam_file -o . --read-tagger-method start
CLAM realigner -i $bam_file -o . --read-tagger-method median
cd ..
done

### 3. peakcaller
s=/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/5_CLAM
mkdir -p $s/result_HuR_input
/public/home/renjun/tool/CLAM-1.2.0/bin/CLAM peakcaller \
-i $s/HuR/unique.sorted.bam $s/HuR/realigned.sorted.bam \
-c $s/input/unique.sorted.bam $s/input/realigned.sorted.bam \
-o $s/result_HuR_input \
--gtf /share/data6/tmp/renjun/Reference/refdata-gex-mm10-2020-A/genes/genes.gtf \
-p 3 \
--binsize 50 \
--qval-cutoff 0.1 

# /public/home/renjun/tool/CLAM-1.2.0/bin/CLAM peakcaller \
# 	-i $s/HA-IP/unique.sorted.bam $s/HA-IP/realigned.sorted.bam \
# 	-c $s/IgG-IP/unique.sorted.bam $s/IgG-IP/realigned.sorted.bam \
# 	-o $s/result_HA-IP_IgG-IP \
# 	--gtf /share/data6/tmp/renjun/Reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
# 	-p 3 \
# 	--binsize 50 \
# 	--qval-cutoff 0.1 

##################################################################
#                                                                
#                   annotation - original                     
#                                                                 
##################################################################

## 1. annotation using CLAM
s=/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/5_CLAM/
/public/home/renjun/tool/CLAM-1.2.0/bin/CLAM peak_annotator \
	-i $s/result_HuR_input/narrow_peak.combined.bed \
	-g mm10 \
	-o $s/result_HuR_input/annotate.txt

## 2. annotation using EnsDb.Hsapiens.v79
library(AnnotationDbi)
# library(EnsDb.Hsapiens.v79)
library(EnsDb.Mmusculus.v79)
peaks = read.table("annotate.txt")
colnames(peaks) = c("chr", "peak_start", "peak_end", "peak_name", "peak_score", "peak_strand",
					          "peak_signal_value","peak_pvalue", "peak_qvalue",
					          "Point-source called for this peak", "Genomic region chromosome", "Genomic region start", "Genomic region end",
					          "gene_id", "quality_score", "Genomic region strand", "Genomic region type")

sapply(1:nrow(peaks), function(i){
  cat(i, "\n")
  AnnotationDbi::select(EnsDb.Mmusculus.v79, keys = peaks$gene_id[i], columns = "SYMBOL", keytype = "GENEID")                 
}) -> tmpg                   
peaks$symbol = tmpg[2,]
write.csv(as.matrix(peaks), file = "v79_annotatedPeak.csv", row.names=FALSE, quote=FALSE)

##################################################################
#                                                                
#                   annotation - original                     
#                                                                 
##################################################################

peaks <- readPeakFile("narrow_peak.combined.bed")

peakAnno <- annotatePeak(peaks, 
                         tssRegion = c(-5000, 5000), 
                         TxDb = txdb, 
                         addFlankGeneInfo = TRUE, 
                         flankDistance = 5000,
                         annoDb = "org.Mm.eg.db")
df <- as.data.frame(as.GRanges(peakAnno))
colnames(df)[10:11] = c("pval","p.adjust")
write.csv(df, file = paste0("chipseeker_annotatedPeak.csv"))
save(df, file = paste0("chipseeker_annotatedPeak.rda"))

clam = df
clam %>% 
  filter(p.adjust<0.05, annotation=="Promoter (<=1kb)") %>% 
  dplyr::select(SYMBOL) %>% 
  unlist %>% 
  unique