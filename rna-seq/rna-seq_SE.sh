#######前提：已经构建好笔记中的分析环境，已经有参考序列的索引，已经有原始fastq文件#######
#######注意虚拟环境！！！！！！！
#######注意控制线程数！！！！！！
#######注意路径！！！！！！！！！

###质控
conda activate rna-seq
mkdir qc
mkdir hisat2
fastqc  /ifs1/User/rna-seq/tomato/AC_heat/raw_data/*.fastq -o qc  -t 28  #改路径
multiqc -d qc/ -o multiqc #生成整合化的质控报告

###过滤
ls -1 /ifs1/User/rna-seq/tomato/AC_heat/raw_data/| while read i ; do echo "fastp -i /ifs1/User/rna-seq/tomato/AC_heat/raw_data/${i} -w 3 -z 4 -q 20 -u 30 -n 10 -o hisat2/${i%.*}_clean.fq.gz" >>fastp.sh;done #改路径
parallel -j 12 -a fastp.sh 
rm -f fastp.sh 


###比对
cd hisat2
ls -1 |while read i;do echo "hisat2 -x /ifs1/User/my_work/rna-seq/ref/S_lycopersicum_chromosomes.4.00.fa -U ${i} -S ${i%_clean*}.sam -p 6;samtools sort -@ 6 ${i%_clean*}.sam -o ${i%_clean*}.sorted.bam;samtools index ${i%_clean*}.sorted.bam">>hisat2.sh;done
parallel -j 6 -a hisat2.sh 
rm -f hisat2.sh


###记数
featureCounts -g gene_id -a /ifs1/User/my_work/rna-seq/ref/ITAG4.0_gene_models.gtf --primary -T 28 -o CountMatrix.txt *.bam #改路径
