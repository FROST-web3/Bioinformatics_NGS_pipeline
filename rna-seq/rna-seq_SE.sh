#!/bin/bash

##定位在原始fa文件的文件夹中

###先单独质控
conda activate rna-seq
mkdir qc
mkdir hisat2
fastqc  *.fq -o qc  -t 28  #改路径和文件名
multiqc -d qc/ -o multiqc #生成整合化的质控报告

###过滤
ls -1 *.fq|while read i; do echo "fastp -i $i -o hisat2/${i%.fq.*}_clean.fq -w 4 -q 20 -u 30 -n 10 --detect_adapter_for_pe " >>fastp.sh;done 
parallel -j 6 -a fastp.sh
rm -f fastp.sh

###比对
cd hisat2
ls -1 |while read i;do echo "hisat2 -x /path/to/reference/S_lycopersicum_chromosomes.4.00.fa -U ${i} -S ${i%_clean*}.sam -p 6;samtools sort -@ 6 ${i%_clean*}.sam -o ${i%_clean*}.sorted.bam;samtools index ${i%_clean*}.sorted.bam">>hisat2.sh;done
parallel -j 6 -a hisat2.sh
rm -f hisat2.sh

###记数
featureCounts -g gene_id -a /path/to/reference/ITAG4.0_gene_models.gtf --primary -T 28 -o CountMatrix.txt *.bam 

###bw可视化文件
mkdir bw
mv *.bam* bw
cd bw
ls -1 *.bam|while read i;do echo "bamCoverage -b $i -o ${i%.sorted*}.bw --binSize 10 --normalizeUsing RPKM --scaleFactor 1 --numberOfProcessors 28">>bw.sh;done;
parallel -j 6 -a bw.sh