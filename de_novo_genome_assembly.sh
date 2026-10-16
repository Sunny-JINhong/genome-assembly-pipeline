#!/bin/bash
# This shell script is for de novo genome assembly
# Usage:./de_novo_genome_assembly.sh <hifi_reads> <illumina_R1> <illumina_R2> <SAMPLE> <output_dir> <genomeSize>

if [ "$#" -ne 6 ]; then
	echo "Usage: $0 <hifi_reads> <illumina_R1> <illumina_R2> <SAMPLE-name> <output_directory> <genomeSize>"
	echo ""
	echo "Arguments:"
	echo " <hifi_reads> Path to PacBio HiFi reads (FASTQ or FASTQ.GZ)"
	echo " <illumina_R1> Path to Illumina paired-end R1 reads"
	echo " <illumina_R2> Path to Illumina paired-end R2 reads"
	echo " <SAMPLE-name> the name of your SAMPLE"
	echo " <output_directory> Directory for all output files"
    echo " <genomeSize> the genome size of your species"
	exit 1
fi

HiFi_reads="$1"
Illumina_R1="$2"
Illumina_R2="$3"
SAMPLE="$4"
outdir="$5"
genomeSize="$6"

THREADS=${THREADS:-16}

mkdir -p ${outdir}
cd ${outdir} || exit 1

LOG="${SAMPLE}.assembly.log"

echo "=======================" | tee -a "${LOG}"
echo "样本：${SAMPLE}" | tee -a "${LOG}"
echo "HiFi：${HiFi_reads}" | tee -a "${LOG}"
echo "Illumina R1：${Illumina_R1}" | tee -a "${LOG}"
echo "Illumina R2：${Illumina_R2}" | tee -a "${LOG}"
echo "输出目录：${outdir}" | tee -a "${LOG}"
echo "=======================" | tee -a "${LOG}"

echo "[$(date)] Step 1: Canu 组装" | tee -a "${LOG}"
canu -p ${SAMPLE} -d ${SAMPLE}_canu_output genomeSize=${genomeSize} useGrid=false maxThreads=$THREADS -pacbio-hifi ${HiFi_reads} >> "${LOG}" 2>&1

CANU_contigs="${SAMPLE}_canu_output/${SAMPLE}.contigs.fasta"

if [ ! -f "${CANU_contigs}" ];then
	echo "错误：Canu 组装失败，未找到 ${CANU_contigs}" | tee -a "${LOG}"
	exit 1
fi

echo "[$(date)] Canu 完成：${CANU_contigs}" | tee -a "${LOG}"

echo "[$(date)] Step 2: Purge_Dups 去冗余" | tee -a "${LOG}"

minimap2 -cx map-pb -t $THREADS ${CANU_contigs} ${HiFi_reads} > ${SAMPLE}_hifi_to_contigs.paf
pbcstat ${SAMPLE}_hifi_to_contigs.paf
split_fa ${CANU_contigs} > ${SAMPLE}.split.fa
minimap2 -x asm5 -D -t $THREADS ${SAMPLE}.split.fa ${SAMPLE}.split.fa > ${SAMPLE}.self_align.paf
calcuts PB.stat > cutoffs.txt 2>calcults.log
purge_dups -2 -T cutoffs.txt -c PB.base.cov ${SAMPLE}.self_align.paf > ${SAMPLE}.dups.bed 2>${SAMPLE}.purge_dups.log
get_seqs -e ${SAMPLE}.dups.bed ${CANU_contigs} purged
mv purged.fa ${SAMPLE}.purged_contigs.fa

PURGED_contigs="${SAMPLE}.purged_contigs.fa"

if [ ! -f "${PURGED_contigs}" ];then
	echo "错误：Purge_Dups 失败" | tee -a "${LOG}"
	exit 1
fi

echo "[$(date)] Purge_Dups 完成：${PURGED_contigs}" | tee -a "${LOG}"

echo "[$(date)] Step 3: BWA + FreeBayes 碱基校正" | tee -a "${LOG}"

bwa index ${PURGED_contigs}
bwa mem -t $THREADS ${PURGED_contigs} ${Illumina_R1} ${Illumina_R2} | samtools sort -@ 4 -m 2G -o ${SAMPLE}_illumina_sorted.bam
samtools index ${SAMPLE}_illumina_sorted.bam
freebayes -f ${PURGED_contigs} -b ${SAMPLE}_illumina_sorted.bam --ploidy 2 --min-coverage 10 -v ${SAMPLE}.variants.vcf
bgzip ${SAMPLE}.variants.vcf
tabix -p vcf ${SAMPLE}.variants.vcf.gz 
bcftools consensus -f ${PURGED_contigs} ${SAMPLE}.variants.vcf.gz > ${SAMPLE}.final_genome.fa

FINAL_genome="${SAMPLE}.final_genome.fa"

if [ ! -f "${FINAL_genome}" ];then
	echo "错误：BWA + FreeBayes 碱基校正 失败" | tee -a "${LOG}"
	exit 1
fi

echo "[$(date)] 碱基校正完成：${FINAL_genome}" | tee -a "${LOG}"

echo "========================================" | tee -a "${LOG}"
echo "[$(date)] 所有步骤完成！" | tee -a "${LOG}"
echo "最终基因组: ${FINAL_genome}" | tee -a "${LOG}"
echo "========================================" | tee -a "${LOG}"
