# De Novo Genome Assembly Pipeline
基于 PacBio HiFi + Illumina 双端测序数据的基因组组装流程（Shell 脚本），包含组装、去冗余、碱基校正三步。

## 流程概述
用canu对三代测序数据进行初步组装，得到contigs；再用Purge_Dups基于比对覆盖度识别并移除冗余序列；最后用二代测序数据进行polish。

## 依赖工具
canu、minimap2、purge_dups、bwa、samtools、freebayes、bcftools、bgzip、tabix

## 可通过conda安装：
```
conda install -c bioconda -c conda-forge canu minimap2 purge_dups bwa samtools freebayes bcftools
```

## 用法：
```
bash de_novo_genome_assembly.sh <hifi_reads> <illumina_R1> <illumina_R2> <SAMPLE> <output_dir> <genomeSize>
```
说明：hifi_reads：PacBio HiFi reads 文件路径（FASTQ 或 FASTQ.GZ）；
     illumina_R1：Illumina 双端测序 R1 文件路径；
     illumina_R2：Illumina 双端测序 R2 文件路径；
     SAMPLE：样本名称（用于输出文件命名）；
     output_dir：输出目录（自动创建）；
     genomeSize：物种基因组大小，如 50m。

## 输出结果
{SAMPLE}.final_genome.fa：最终基因组；
{SAMPLE}.assembly.log：全流程运行日志；

## 注意事项
genomeSize 需根据物种实际情况填写。
