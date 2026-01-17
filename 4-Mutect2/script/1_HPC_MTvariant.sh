#!/usr/bin/env bash

#############################################
# 本脚本主要用于线粒体基因组的变异检测流程。
# 说明：在不改变原始脚本功能和内容的前提下，对其进行了中文注释以使流程更清晰。
#############################################

file=$1

# 设置输出目录、参考基因组及其它路径
DIR=/public/home/tangting/MTvariant/20231113/01.output
REF=/public/share/wchirdzhq2022/zhaoxingkai/RefDIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
RefDIR=/public/share/wchirdzhq2022/zhaoxingkai/RefDIR
thread_num=32

# 加载所需软件及其依赖
module load apps/GATK/4.1.2.0
module load apps/singularity/3.8.3

gatk="singularity exec -B /public:/public /public/share/wchirdzhq2022/tangting/biosoft/varflow-mt.sif gatk"
picard="singularity exec -B /public:/public /public/share/wchirdzhq2022/tangting/biosoft/varflow-mt.sif java -jar /opt/picard.jar"
bwa="singularity exec -B /public:/public /public/share/wchirdzhq2022/tangting/biosoft/varflow-mt.sif bwa"
samtools="singularity exec -B /public:/public /public/share/wchirdzhq2022/tangting/biosoft/varflow-mt.sif samtools"
bgzip="singularity exec -B /public:/public /public/share/wchirdzhq2022/tangting/biosoft/varflow-mt.sif bgzip"
gzip="singularity exec -B /public:/public /public/share/wchirdzhq2022/tangting/biosoft/varflow-mt.sif gzip"
bcftools="singularity exec -B /public:/public /public/share/wchirdzhq2022/tangting/biosoft/varflow-mt.sif bcftools"
vep="singularity exec -B /public:/public /public/share/wchirdzhq2022/tangting/biosoft/varflow-mt.sif vep"
vcfanno="singularity exec -B /public:/public /public/share/wchirdzhq2022/tangting/biosoft/varflow-mt.sif vcfanno"

###########################################
# 下方为整理线粒体基因组的参考文件的过程 (示例性命令，已被注释)
# 其中将chrM部分提取并索引，供后续调用。
# grep -n ">chrM" $RefDIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
# grep -n ">chr1_KI270706v1_random" $RefDIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
# sed -n '1,333p' $RefDIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta > $RefDIR/ucsc.hg38.chrM.fa
# time $bwa index $RefDIR/ucsc.hg38.chrM.fa
# time $samtools dict $RefDIR/ucsc.hg38.chrM.fa > $RefDIR/ucsc.hg38.chrM.fa.dict
# time $samtools faidx $RefDIR/ucsc.hg38.chrM.fa
# time $bwa index -a bwtsw $RefDIR/ucsc.hg38.chrM.fa
###########################################

cat $file | while read SampleID Bam; do

  # 为每个Sample创建输出目录
  mkdir -p ${DIR}/${SampleID}
  OutputDIR=${DIR}/${SampleID}
  cd $OutputDIR

  ########################################################
  # 第4步：提取chrM上的reads
  # 使用GATK的PrintReads命令，主要参数：
  #   -I：输入BAM文件
  #   -L：指定只处理chrM染色体
  #   --read-filter：过滤mate在同一个contig或未比对的read等。
  # 输出：SampleID.chrM.bam
  ########################################################
  time $gatk --java-options "-Xmx20g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  PrintReads \
  -I $Bam \
  -L chrM \
  --read-filter MateOnSameContigOrNoMappedMateReadFilter \
  --read-filter MateUnmappedAndUnmappedReadFilter \
  -O $OutputDIR/${SampleID}.chrM.bam

  ########################################################
  # 对提取出的chrM BAM文件进行排序
  ########################################################
  time $samtools sort \
    -@ $thread_num \
    -O bam $OutputDIR/${SampleID}.chrM.bam \
    > $OutputDIR/${SampleID}.chrM.sorted.bam

  ########################################################
  # 使用Picard构建BAM索引
  ########################################################
  time $picard \
  BuildBamIndex \
  INPUT=${OutputDIR}/${SampleID}.chrM.sorted.bam \
  OUTPUT=${OutputDIR}/${SampleID}.chrM.sorted.bai

  ########################################################
  # 第5步：计算线粒体基因组上的平均覆盖深度
  # 说明：
  #   1）参考基因组必须与产生BAM的基因组一致
  #   2）interval文件头部必须与BAM头部一致，否则结果可能为0覆盖
  ########################################################
  time $picard \
  CollectWgsMetrics \
  I=$OutputDIR/$SampleID.chrM.sorted.bam \
  O=$OutputDIR/$SampleID.chrM.CollectWgsMetrics.txt \
  R=$REF \
  INTERVALS=${RefDIR}/hg38.chrM.interval

  ########################################################
  # 第6步：使用Mutect2进行线粒体非控制区变异检测 (non-control region)
  #   --mitochondria-mode
  #   --annotation StrandBiasBySample
  #   --max-reads-per-alignment-start 75
  #   --max-mnp-distance 0
  ########################################################
  time $gatk --java-options "-Xmx20g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  Mutect2 \
  -R $RefDIR/ucsc.hg38.chrM.fa \
  -L chrM \
  --mitochondria-mode \
  --annotation StrandBiasBySample \
  --max-reads-per-alignment-start 75 \
  --max-mnp-distance 0 \
  -I $OutputDIR/$SampleID.chrM.sorted.bam \
  -O $OutputDIR/$SampleID.chrM.NonControlRegion.vcf

  ########################################################
  # 清理VCF头部信息，去除不需要的contig信息
  # 并压缩索引
  ########################################################
  sed -e '/^##contig\(.*\)random\(.*\)/d' \
      -e '/^##contig\(.*\)Un\(.*\)/d' \
      -e '/^##contig\(.*\)_alt\(.*\)/d' \
      -e '/^##contig\(.*\)HLA\(.*\)/d' \
      -e '/^##contig\(.*\)EBV\(.*\)/d' \
      -e '/^##contig\(.*\)chr1\(.*\)/d' \
      -e '/^##contig\(.*\)chr2\(.*\)/d' \
      -e '/^##contig\(.*\)chr3\(.*\)/d' \
      -e '/^##contig\(.*\)chr4\(.*\)/d' \
      -e '/^##contig\(.*\)chr5\(.*\)/d' \
      -e '/^##contig\(.*\)chr6\(.*\)/d' \
      -e '/^##contig\(.*\)chr7\(.*\)/d' \
      -e '/^##contig\(.*\)chr8\(.*\)/d' \
      -e '/^##contig\(.*\)chr9\(.*\)/d' \
      -e '/^##contig\(.*\)chrX\(.*\)/d' \
      -e '/^##contig\(.*\)chrY\(.*\)/d' \
      $OutputDIR/$SampleID.chrM.NonControlRegion.vcf > $OutputDIR/$SampleID.chrM.NonControlRegion.reform.vcf

  $bgzip -f -@ $thread_num $OutputDIR/$SampleID.chrM.NonControlRegion.reform.vcf
  $bcftools index --threads $thread_num -f -t $OutputDIR/$SampleID.chrM.NonControlRegion.reform.vcf.gz

  ########################################################
  # 第7步：将chrM的reads比对到“shifted”后的chrM参考基因组上
  # GATK RevertSam 用于还原BAM文件到未比对状态
  # Picard SamToFastq 将BAM转换成Fastq用于重新比对
  ########################################################
  time $gatk RevertSam \
  -I $OutputDIR/$SampleID.chrM.sorted.bam \
  -O $OutputDIR/$SampleID.chrM.sorted.reverted.bam

  time $picard SamToFastq \
  I=$OutputDIR/$SampleID.chrM.sorted.reverted.bam \
  FASTQ=$OutputDIR/${SampleID}_R1.chrM.reverted.fastq.gz \
  SECOND_END_FASTQ=$OutputDIR/${SampleID}_R2.chrM.reverted.fastq.gz

  RG="@RG\\tID:group1\\tSM:$SampleID\\tPL:illumina\\tLB:CasCADE\\tPU:unit1"

  time $bwa mem \
  -t $thread_num \
  -M \
  -R $RG $RefDIR/ucsc.hg38.chrM.shifted.fa \
  $OutputDIR/${SampleID}_R1.chrM.reverted.fastq.gz \
  $OutputDIR/${SampleID}_R2.chrM.reverted.fastq.gz \
  > $OutputDIR/$SampleID.chrM.shifted.sorting

  time $samtools sort --threads $thread_num -O bam $OutputDIR/$SampleID.chrM.shifted.sorting > $OutputDIR/$SampleID.chrM.shifted.sorted.bam

  ########################################################
  # 为BAM文件添加或替换Read Group信息
  ########################################################
  time $picard AddOrReplaceReadGroups \
  I=$OutputDIR/$SampleID.chrM.shifted.sorted.bam \
  O=$OutputDIR/$SampleID.chrM.shifted.sorted.RG.bam \
  RGID=4 \
  RGLB=lib1 \
  RGPL=illumina \
  RGPU=unit1 \
  RGSM=$SampleID

  ########################################################
  # 去重复并建立索引
  ########################################################
  time $picard \
  MarkDuplicates \
  INPUT=$OutputDIR/$SampleID.chrM.shifted.sorted.RG.bam \
  OUTPUT=$OutputDIR/$SampleID.chrM.shifted.dedup.bam \
  METRICS_FILE=$OutputDIR/$SampleID.chrM.shifted.dedup.metrics

  time $picard \
  BuildBamIndex \
  INPUT=$OutputDIR/$SampleID.chrM.shifted.dedup.bam \
  OUTPUT=$OutputDIR/$SampleID.chrM.shifted.dedup.bai

  ########################################################
  # 第8步：计算转换后线粒体基因组的平均深度
  ########################################################
  time $picard \
  CollectWgsMetrics \
  I=$OutputDIR/$SampleID.chrM.shifted.dedup.bam \
  O=$OutputDIR/$SampleID.chrM.shifted.CollectWgsMetrics.txt \
  R=$RefDIR/ucsc.hg38.chrM.shifted.fa

  ########################################################
  # 第9步：使用Mutect2在“shifted”参考基因组上进行控制区变异检测
  ########################################################
  time $gatk \
  Mutect2 \
  -R $RefDIR/ucsc.hg38.chrM.shifted.fa \
  -L chrM \
  --mitochondria-mode \
  --annotation StrandBiasBySample \
  --max-reads-per-alignment-start 75 \
  --max-mnp-distance 0 \
  -I $OutputDIR/$SampleID.chrM.shifted.dedup.bam \
  -O $OutputDIR/$SampleID.chrM.ControlRegion.vcf.gz

  ########################################################
  # 第10步：使用Picard的LiftoverVcf将变异位点重新映射回原chrM坐标
  ########################################################
  time $picard LiftoverVcf \
  I=$OutputDIR/$SampleID.chrM.ControlRegion.vcf.gz \
  O=$OutputDIR/$SampleID.chrM.ControlRegion.liftover.vcf.gz \
  CHAIN=$RefDIR/ucsc.hg38.chrM.shifted.chain \
  REJECT=$OutputDIR/$SampleID.chrM.ControlRegion.liftover.rejected_variants.vcf.gz \
  R=$RefDIR/ucsc.hg38.chrM.fa

  ########################################################
  # 第11步：合并非控制区和控制区的变异检测结果
  ########################################################
  time $picard MergeVcfs \
  I=$OutputDIR/$SampleID.chrM.NonControlRegion.reform.vcf.gz \
  I=$OutputDIR/$SampleID.chrM.ControlRegion.liftover.vcf.gz \
  D=/public/home/tangting/hg38/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.dict \
  O=$OutputDIR/$SampleID.chrM.raw.vcf.gz

  ########################################################
  # 第12步：对合并后的变异进行过滤
  #   --max-alt-allele-count 4  限制备用等位基因数量
  #   --mitochondria-mode      线粒体模式
  #   --min-allele-fraction 0.01  等位基因频率阈值
  ########################################################
  time $gatk FilterMutectCalls \
  -R $RefDIR/ucsc.hg38.chrM.fa \
  -V $OutputDIR/$SampleID.chrM.raw.vcf.gz \
  --stats $OutputDIR/$SampleID.chrM.ControlRegion.vcf.gz.stats \
  --max-alt-allele-count 4 \
  --mitochondria-mode \
  --min-allele-fraction 0.01 \
  -O $OutputDIR/$SampleID.chrM.filtered.vcf.gz

done
