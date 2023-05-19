
# https://link.springer.com/protocol/10.1007/978-1-4939-2291-8_18
# 环境，59.77.41.69，conda cR4.1

#-----------------------------------------------------------------
#                                                                
#                          create index                     
#                                                                 
#-----------------------------------------------------------------

# --runThreadN：线程数
# --runMode genomeGenerate：构建基因组索引
# --genomeDir：索引目录。（index_dir一定要是存在的文件夹，需提前建好）
# --genomeFastaFiles：基因组文
# --sjdbGTFfile：基因组注释文件
# --sjdbOverhang：reads长度减1

STAR  \
--runMode genomeGenerate \
--genomeDir index \
--runThreadN 20 \
--genomeFastaFiles /share/data6/tmp/renjun/Reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
--sjdbGTFfile /share/data6/tmp/renjun/Reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
--sjdbOverhang 149

# 查看序列长度
zcat HuR_1.fq.gz | awk '{if(NR%4==2) print length($1)}' | head

#-----------------------------------------------------------------
#                                                                
#                         fastqc & multiqc               
#                                                                 
#-----------------------------------------------------------------

len=149
fq_paths=("/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/0_data/input/soapnuke/clean/input/"  \
		  "/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/0_data/ip/soapnuke/clean/HuR/" )

for p in ${fq_paths[*]};do
	echo `ls $p`
	fq1=`ls $p*1.fq.gz`
	fq2=`ls $p*2.fq.gz`
	outfile=$(echo $(basename $fq1) | cut -d \_ -f 1)
	mkdir $outfile && cd $outfile

	fastqc -o ./ -t 20 $fq1 $fq2
	cd ..
done

# 主要是包括前面的各种选项和最后面的可以加入N个文件
# -o --outdir FastQC生成的报告文件的储存路径，生成的报告的文件名是根据输入来定的
# --extract 生成的报告默认会打包成1个压缩文件，使用这个参数是让程序不打包
# -t --threads 选择程序运行的线程数，每个线程会占用250MB内存，越多越快咯
# -c --contaminants 污染物选项，输入的是一个文件，格式是Name [Tab] Sequence，里面是可能的污染序列，如果有这个选项，FastQC会在计算时候评估污染的情况，并在统计的时候进行分析，一般用不到
# -a --adapters 也是输入一个文件，文件的格式Name [Tab] Sequence，储存的是测序的adpater序列信息，如果不输入，目前版本的FastQC就按照通用引物来评估序列时候有adapter的残留
# -q --quiet 安静运行模式，一般不选这个选项的时候，程序会实时报告运行的状况。

multiqc `find . -type f |grep .*zip` -o multiqc

#-----------------------------------------------------------------
#                                                                
#                            sequence                     
#                                                                 
#-----------------------------------------------------------------

len=149
genome_dir="/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/2_genome"
fq_paths=("/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/0_data/input/soapnuke/clean/input/"  \
		  "/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/0_data/ip/soapnuke/clean/HuR/" )

for p in ${fq_paths[*]};do
	echo `ls $p`
	fq1=`ls $p*1.fq.gz`
	fq2=`ls $p*2.fq.gz`
	outfile=$(echo $(basename $fq1) | cut -d \_ -f 1)
	mkdir $outfile && cd $outfile

	STAR \
	--twopassMode Basic \
	--genomeDir $genome_dir \
	--readFilesIn $fq1 $fq2 \
	--readFilesCommand  zcat \
	--sjdbOverhang  $len \
	--outFileNamePrefix $outfile \
	--runThreadN 20 \
	--outSAMtype BAM SortedByCoordinate

	cd ..
done

# Log.out（一些进程日志，方便之后debug）
# Log.progress.out（进程中的统计工作，每一分钟记录一次，包括当前比对率等信息）
# Log.final.out（最终的比对结果统计文件，对之后的质控很有用！）
# SAM文件（Aligned.out.sam）
# Unsorted或sorted-by-coordinate的BAM文件
# 通过调整--outSAMtype 可以得到不同形式的bam文件Aligned.out.bam或Aligned.sortedByCoord.out.bam两类文件
# splice junction，SJ.out.tab文件中包含了可信度很高的剪切位点信息，具体每一列的意义参考帮助文档

#-----------------------------------------------------------------
#                                                                
#               picard duplicate + samtools index                 
#                                                                 
#-----------------------------------------------------------------

cd 4_bam
for p in $(ls);do
	echo $p
	cd $p
	bamfile=`find . -type f |grep .out.bam`
	outbamfile=$(basename -s .bam $bamfile).redup.bam

	# duplicated
	picard MarkDuplicates \
		REMOVE_DUPLICATES=true \
	    I=$bamfile \
	    O=$outbamfile \
	    M=picard.txt

	# make index
	# 方便用IGV查看
	samtools index $outbamfile

	cd ..
done

#-----------------------------------------------------------------
#                                                                
#                        bam文件进一步操做                     
#                                                                 
#-----------------------------------------------------------------

### 检查bam文件的完整性
samtools quickcheck inputAligned.sortedByCoord.out.bam && echo "ok" || echo $i error

## 统计染色体数量
samtools view x.bam | awk '{print $3}' | uniq -c

#-----------------------------------------------------------------
#                                                                
#                    Piranha - peak calling                     
#                                                                 
#-----------------------------------------------------------------

### 1. 将bam文件转为bed文件
bam_paths=("/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/4_bam/input/"  \
		   "/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/4_bam/HA-IP/" \
		   "/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/4_bam/IgG-IP/")
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
cd 5_piranha
control_1="../4_bam/input/inputAligned.sortedByCoord.out.redup.chr.bed"
control_2="../4_bam/IgG-IP/IgG-IPAligned.sortedByCoord.out.redup.chr.bed"
treat="../4_bam/HA-IP/HA-IPAligned.sortedByCoord.out.redup.chr.bed"

Piranha $treat $control_1 -s -z 100 -p 0.05 -o HA-IP_input.bed
Piranha $treat $control_2 -s -b 100 -a 100 -p 0.05 -o HA-IP_IgG-IP.bed

### 4. debug
# 提示协变量必须要和treat的基因位点对应，确实不对应
# 只能用下面的方法操作然后注释
Piranha $treat -s -z 1 -p 0.05 -o HA-IP.bed

#-----------------------------------------------------------------
#                                                                
#                      CLAM - peak calling                     
#                                                                 
#-----------------------------------------------------------------

### https://github.com/Xinglab/CLAM
bam_paths=("/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/4_bam/input/"  \
		   "/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/4_bam/HA-IP/" \
		   "/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/4_bam/IgG-IP/")

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
s=/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/6_CLAM
/public/home/renjun/tool/CLAM-1.2.0/bin/CLAM peakcaller \
	-i $s/HA-IP/unique.sorted.bam $s/HA-IP/realigned.sorted.bam \
	-c $s/input/unique.sorted.bam $s/input/realigned.sorted.bam \
	-o $s/result_HA-IP_input \
	--gtf /share/data6/tmp/renjun/Reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
	-p 3 \
	--binsize 50 \
	--qval-cutoff 0.1 

/public/home/renjun/tool/CLAM-1.2.0/bin/CLAM peakcaller \
	-i $s/HA-IP/unique.sorted.bam $s/HA-IP/realigned.sorted.bam \
	-c $s/IgG-IP/unique.sorted.bam $s/IgG-IP/realigned.sorted.bam \
	-o $s/result_HA-IP_IgG-IP \
	--gtf /share/data6/tmp/renjun/Reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
	-p 3 \
	--binsize 50 \
	--qval-cutoff 0.1 

### 4. annotation + EnsDb.Hsapiens.v79
# s=/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/6_CLAM
# /public/home/renjun/tool/CLAM-1.2.0/bin/CLAM peak_annotator \
# 	-i $s/result_HA-IP_input/narrow_peak.combined.bed \
# 	-g hg19 \
# 	-o $s/result_HA-IP_input/annotate.txt

# /public/home/renjun/tool/CLAM-1.2.0/bin/CLAM peak_annotator \
# 	-i $s/result_HA-IP_IgG-IP/narrow_peak.combined.bed \
# 	-g hg19 \
# 	-o $s/result_HA-IP_IgG-IP/annotate.txt

# library(AnnotationDbi)
# library(EnsDb.Hsapiens.v79)
# peaks = read.table("annotate.txt")
# colnames(peaks) = c("chr", "peak_start", "peak_end", "peak_name", "peak_score", "peak_strand", 
# 					"peak_signal_value","peak_pvalue", "peak_qvalue",
# 					"Point-source called for this peak", "Genomic region chromosome", "Genomic region start", "Genomic region end",
# 					"gene_id", "quality_score", "Genomic region strand", "Genomic region type")
# geneDF = data.frame(row.names = 1:nrow(peaks),
#                     index = 1:nrow(peaks),
#                     geneid = peaks$gene_id,
#                     stringsAsFactors = FALSE)
# tmp = AnnotationDbi::select(EnsDb.Hsapiens.v79, keys = geneDF$geneid, columns = "SYMBOL", keytype = "GENEID")
# ind = which(!is.na(tmp$SYMBOL))
# df = merge(geneDF, tmp[ind,], 
#            by.x = "geneid", by.y = "GENEID", 
#            all.x = TRUE, all.y = TRUE, all = TRUE, sort = FALSE)
# df = df[order(df$index),]
# df = merge(peaks,df, by.x="gene_id",by.y="geneid")
# write.csv(df, file = "v79_annotatedPeak.csv", row.names=FALSE, quote=FALSE)

#-----------------------------------------------------------------
#                                                                
#                      Ripseeker - peak calling                     
#                                                                 
#-----------------------------------------------------------------
                                                               
## 1. 按照染色体分割bam文件     
bam_paths=("/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/4_bam/input/"  \
		   "/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/4_bam/HA-IP/" \
		   "/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq/4_bam/IgG-IP/")
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

## 3. seekRIP
## 4. annot
## 5. enrichment

#-----------------------------------------------------------------
#                                                                
#                            ChIPseeker                     
#                                                                 
#-----------------------------------------------------------------
  
# peaks <- readPeakFile("HA-IP.bed")
# peaks <- readPeakFile("narrow_peak.combined.bed")

library(ChIPseeker)
library(GenomicFeatures)
txdb <- makeTxDbFromGFF(file = "/share/data6/tmp/renjun/Reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
                        format = "gtf", organism = "Homo sapiens")
peakAnno <- annotatePeak(peaks, 
                         tssRegion = c(-5000, 5000), 
                         TxDb = txdb, 
                         addFlankGeneInfo = TRUE, 
                         flankDistance = 5000,
                         annoDb = "org.Hs.eg.db")
df <- as.data.frame(as.GRanges(peakAnno))
colnames(df)[6:12] = c("peak_name", "peak_score", "peak_strand", "peak_signal_value",
					   "peak_pvalue", "peak_qvalue",
					   "Point-source called for this peak")

write.csv(df, file = paste0("chipseeker_annotatedPeak.csv"))
save(df, file = paste0("chipseeker_annotatedPeak.rda"))