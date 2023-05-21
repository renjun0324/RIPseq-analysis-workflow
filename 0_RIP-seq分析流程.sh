
# https://link.springer.com/protocol/10.1007/978-1-4939-2291-8_18
# 环境，59.77.41.69，conda cR4.1

##################################################################
#                                                                
#                      fastqc & multiqc               
#                                                                 
##################################################################

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

##################################################################
#                                                                
#                       create index                     
#                                                                 
##################################################################

## 1. 查看序列长度
zcat HuR_1.fq.gz | awk '{if(NR%4==2) print length($1)}' | head
zcat input_1.fq.gz | awk '{if(NR%4==2) print length($1)}' | head

## 2. STAR构建索引
# --runThreadN：线程数
# --runMode genomeGenerate：构建基因组索引
# --genomeDir：索引目录。（index_dir一定要是存在的文件夹，需提前建好）
# --genomeFastaFiles：基因组文
# --sjdbGTFfile：基因组注释文件
# --sjdbOverhang：reads长度减1

# STAR  \
# --runMode genomeGenerate \
# --genomeDir index \
# --runThreadN 20 \
# --genomeFastaFiles /share/data6/tmp/renjun/Reference/refdata-gex-mm10-2020-A/fasta/genome.fa \
# --sjdbGTFfile /share/data6/tmp/renjun/Reference/refdata-gex-mm10-2020-A/genes/genes.gtf \
# --sjdbOverhang 149

## 3. bowtie2构建索引
bowtie2-build --threads 20 /share/data6/tmp/renjun/Reference/refdata-gex-mm10-2020-A/fasta/genome.fa genome

##################################################################
#                                                                
#                           sequence                     
#                                                                 
##################################################################

## 1. star比对（可直接对输出结果进行排序）
# len=149
# genome_dir="/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/2_genome_star/"
# fq_paths=("/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/0_data/input/soapnuke/clean/input/"  \
# 		  "/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/0_data/ip/soapnuke/clean/HuR/" )
# 
# for p in ${fq_paths[*]};do
# 	echo `ls $p`
# 	fq1=`ls $p*1.fq.gz`
# 	fq2=`ls $p*2.fq.gz`
# 	outfile=$(echo $(basename $fq1) | cut -d \_ -f 1)
# 	mkdir $outfile && cd $outfile
# 
# 	STAR \
# 	--twopassMode Basic \
# 	--genomeDir $genome_dir \
# 	--readFilesIn $fq1 $fq2 \
# 	--readFilesCommand  zcat \
# 	--sjdbOverhang  $len \
# 	--outFileNamePrefix $outfile \
# 	--runThreadN 20 \
# 	--outSAMtype BAM SortedByCoordinate
# 
# 	cd ..
# done

# Log.out（一些进程日志，方便之后debug）
# Log.progress.out（进程中的统计工作，每一分钟记录一次，包括当前比对率等信息）
# Log.final.out（最终的比对结果统计文件，对之后的质控很有用！）
# SAM文件（Aligned.out.sam）
# Unsorted或sorted-by-coordinate的BAM文件
# 通过调整--outSAMtype 可以得到不同形式的bam文件Aligned.out.bam或Aligned.sortedByCoord.out.bam两类文件
# splice junction，SJ.out.tab文件中包含了可信度很高的剪切位点信息，具体每一列的意义参考帮助文档

## 2. bowtie2比对
genome_dir="/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/2_genome_bowtie2/mm10"
fq_paths=("/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/0_data/input/soapnuke/clean/input/"  \
		  "/share/data6/tmp/renjun/ExtraJobs/202305_RIPseq_analysis/0_data/ip/soapnuke/clean/HuR/" )

for p in ${fq_paths[*]};do
	echo `ls $p`
	fq1=`ls $p*1.fq.gz`
	fq2=`ls $p*2.fq.gz`
	outfile=$(echo $(basename $fq1) | cut -d \_ -f 1)
	mkdir $outfile && cd $outfile

  bowtie2 --end-to-end -p 20 -x $genome_dir -1 $fq1 -2 $fq2 2> ${outfile}.bowtie2 | samtools view -bS - >${outfile}.bam  
  
	cd ..
done
  
## 备用参数
# -I 10 
# -X 700 
# --very-sensitive \
# --no-mixed \
# --no-discordant \
# --phred33 \

## 3. bowtie2之后排序 - 必须排序，否则去重没输出
samtools sort HuR.bam -o HuR.sorted.bam
samtools sort input.bam -o input.sorted.bam

##################################################################
#                                                                
#               picard duplicate + samtools index                 
#                                                                 
##################################################################

cd 3_sequence_bowtie2
for p in $(ls);do
	echo $p
	cd $p
	bamfile=`find . -type f |grep sorted`
	outbamfile=$(basename -s .bam $bamfile).redup.bam

	# duplicated
  picard MarkDuplicates REMOVE_DUPLICATES=true I=$bamfile O=$outbamfile M=picard.txt

	# make index
	# 方便用IGV查看
	samtools index $outbamfile

	cd ..
done

# 通过对两个组别进行验证，我发现如果没有sort，就不会去重成功

##################################################################
#                                                                
#                     检查bam文件的完整性                     
#                                                                 
##################################################################

### 检查bam文件的完整性
samtools quickcheck inputAligned.sortedByCoord.out.bam && echo "ok" || echo $i error
samtools quickcheck input.bam && echo "ok" || echo $i error

## 统计染色体数量
samtools view input.bam | awk '{print $3}' | uniq -c

##################################################################
#                                                                
#                         peak calling                    
#                                                                 
##################################################################

# 1. piranha
# 2. CLAM
# 3. ripseeker

#################################################################
#                                                                
#                      查看方法之间的交集                    
#                                                                 
##################################################################

# Mus musculus
# Homo sapiens

clam = read.csv("5_CLAM/result_HuR_input/chipseeker_annotatedPeak.csv")
rip = read.csv("7_RIPSeeker/chipseeker_annotatedPeak.csv")
pir = read.csv("4_piranha/chipseeker_annotatedPeak.csv")

clam %>% 
  filter(p.adjust<0.05, annotation=="Promoter (<=1kb)") %>% 
  dplyr::select(SYMBOL) %>% 
  unlist %>% 
  unique -> c

pir %>% 
  filter(p.adjust<0.01, annotation=="Promoter (<=1kb)") %>% 
  dplyr::select(SYMBOL) %>% 
  unlist %>% 
  unique -> b
  
rip %>% 
  filter(pvalAdj<0.05, annotation=="Promoter (<=1kb)") %>% 
  dplyr::select(SYMBOL) %>% 
  unlist %>% 
  unique -> a
  
Reduce(intersect, list(a,b,c))