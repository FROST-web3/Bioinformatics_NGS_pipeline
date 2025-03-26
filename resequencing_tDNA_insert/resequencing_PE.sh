#给tDNA单独建立索引

#bowtie2比对
ls -1 *.gz|xargs -n2 |while read {i,j};do echo "bowtie2 -p 28 -x /path/to/reference/tDNA_index/tdna_GFP -S ${i%.gz}.sam -1 $i -2 $j --local --very-sensitive-local">>bowtie2.sh;done;
sh bowtie2.sh  

#提取比对上的reads
samtools view -F 4 *.sam -@ 28 > matched_reads.sam

#筛选部分是tDNA，部分是番茄序列的reads
awk '$6 ~ /^[0-9]+S[0-9]+M$/ || $6 ~ /^[0-9]+M[0-9]+S$/ {print $0}' matched_reads.sam > mix_reads.sam

#提取出可能来自番茄基因组的序列，保存为fa文件
awk '
{
    if ($6 ~ /^[0-9]+S/) {
        match($6, /^([0-9]+)S/, arr)
        soft_clip_length = arr[1]
        print ">" $1 "_left"
        print substr($10, 1, soft_clip_length)
    } else if ($6 ~ /[0-9]+S$/) {
        match($6, /([0-9]+)S$/, arr)
        soft_clip_length = arr[1]
        print ">" $1 "_right"
        print substr($10, length($10) - soft_clip_length + 1)
    }
}' mix_reads.sam >tomato_genome_parts.fa

#过滤过短的序列
awk '
BEGIN {FS = "\n"; RS = ">"}
NR > 1 {
    seq = $2;
    gsub(/\n/, "", seq);  # 移除序列中可能的换行符
    if (length(seq) >= 50) {
        print ">" $1 "\n" seq
    }
}' tomato_genome_parts.fa > filtered_tomato_parts.fa



#比对
nohup blastn -query filtered_tomato_parts.fa -db /path/to/reference/S_lycopersicum_chromosomes.4.00 -out blast_results.txt -outfmt 6 -evalue 1e-5 -num_threads 28 &

#只要完全的匹配
awk '$3==100{print($0)}' blast_results.txt >filtered_blast_results.txt

#看看都比对在哪里了
awk '{count[$2]++} END {print "Number of unique values: " length(count); for (value in count) print value ": " count[value]}' filtered_blast_results.txt | sort -k2 -nr >information.txt

#blast结果转化为bed，用于IGV可视乎
awk 'BEGIN{OFS="\t"} {if($9<$10) print $2,$9,$10,$1,".",$9<$10?"+":"-"; else print $2,$10,$9,$1,".",$9<$10?"+":"-"}' filtered_blast_results.txt > blast_results.bed
sort -k1,1 -k2,2n blast_results.bed > sorted_blast_results.bed






###############
#the other way#
###############

#blastn建立索引
makeblastdb -dbtype nucl -in SLM_r2.0.pmol.fasta -input_type fasta -parse_seqids -out SLM_r2.0.pmol

#bowtie建立索引
bowtie2-build -f hg19.fna human

#如果以番茄为参考基因组，用此命令筛选read
samtools view -F 4 -e 'cigar =~ "^[0-9]+S[0-9]+M$" || cigar =~ "^[0-9]+M[0-9]+S$"' *.sam -@ 28 > one_end_soft_clipped_reads.sam

#再截取出比对上的一端，过滤掉过长的，保证S足够长，然后blast
awk '
{
    if ($6 ~ /^[0-9]+S/) {
        match($6, /^([0-9]+)S/, arr)
        soft_clip_length = arr[1]
        print ">" $1 "_matched"
        print substr($10, soft_clip_length + 1)
    } else if ($6 ~ /[0-9]+S$/) {
        match($6, /([0-9]+)S$/, arr)
        soft_clip_length = arr[1]
        print ">" $1 "_matched"
        print substr($10, 1, length($10) - soft_clip_length)
    }
}' one_end_soft_clipped_reads.sam > matched_portions.fa

awk '
/^>/ {header=$0; next}
{
    if (length($0) <= 80) {
        print header
        print $0
    }
}' matched_portions.fa > filtered_80bp.fa


nohup blastn -query filtered_80bp.fa -db /path/to/reference/S_lycopersicum_chromosomes.4.00 -out blast_results.txt -outfmt 6 -evalue 1e-5 -num_threads 28 &

#只要完全的匹配
awk '$3==100{print($0)}' blast_results.txt >filtered_blast_results.txt

#看看都比对在哪里了
awk '{count[$2]++} END {print "Number of unique values: " length(count); for (value in count) print value ": " count[value]}' filtered_blast_results.txt | sort -k2 -nr >information.txt

#blast结果转化为bed，用于IGV可视乎
awk 'BEGIN{OFS="\t"} {if($9<$10) print $2,$9,$10,$1,".",$9<$10?"+":"-"; else print $2,$10,$9,$1,".",$9<$10?"+":"-"}' filtered_blast_results.txt > blast_results.bed
sort -k1,1 -k2,2n blast_results.bed > sorted_blast_results.bed