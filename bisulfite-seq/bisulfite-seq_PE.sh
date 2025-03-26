#!/bin/bash
#定位在sra文件目录中，提前进入虚拟环境bisulfite-seq

process_methylation_data() {
    local input_dir=$1
    local initial_dir=$(pwd) 
    
    echo "Processing directory: $input_dir"
    cd "$input_dir"

    # 提前质控,不包含质控
    fasterq-dump *.sra -e 28


    # fastp
    mkdir cleanData/
    ls -1 *.fastq |xargs -n2|while read i j; do 
        echo "fastp -i $i -I $j -o cleanData/${i%.fa*}_clean.fq -O cleanData/${j%.fa*}_clean.fq --detect_adapter_for_pe -f 6 -t 6 -F 6 -T 6 -w 28 -q 18 -u 30 -n 20 --length_required 30 ">>fastp.sh
    done
    sh fastp.sh 
    rm -rf *.fastq

    # bismark比对
    cd cleanData/
    ls -1 *.fq |xargs -n2|while read i j; do 
        echo "bismark --genome /path/to/reference/ref --parallel 5 -o bismark_alignment/ -1 $i -2 $j --bowtie2 --gzip">>bismark.sh
    done
    sh bismark.sh

    # 根据建库方法选择是否去重
    cd bismark_alignment/
    for i in `ls -1 *.bam`; do 
        echo "deduplicate_bismark --output_dir deduplicate/  $i">>deduplicate_bismark.sh
    done
    sh deduplicate_bismark.sh
    rm -rf *.bam

    # 统计甲基化信息,并生成全基因组C报告文件
    cd deduplicate/
    ulimit -n 4096
    ls -1 *.bam|while read i; do 
        echo "bismark_methylation_extractor --bedGraph --CX --cytosine_report --genome_folder /path/to/reference/ref --output_dir methylation_extractor/ --parallel 10 --buffer_size 10G $i">>methylation_extractor.sh
    done
    sh methylation_extractor.sh

    # 将报告文件转化为IGV格式用于可视化
    cd methylation_extractor/
    mkdir CX_report
    mv *.CX_report.txt CX_report/
    rm -f *
    cd CX_report/

    # 转化全基因组报告文件为bedgraph文件,并转化为bw可视化文件
    convert_to_bedgraph() {
        local input_file=$1
        local output_file=$2

        echo "track type=bedGraph" > "$output_file"

        awk '
        BEGIN {OFS="\t"}
        !/^#/ && NF == 7 && $4 != 0 {
            chrom = $1
            pos = $2
            strand = $3
            meth = $4
            unmeth = $5
            context = $6
            trinucleotide = $7
            total = meth + unmeth
            if (total > 0) {
                level = meth / total
                if (strand == "-") {
                    level = -level
                }
                printf "%s\t%d\t%d\t%.4f\n", chrom, pos-1, pos, level
            }
        }
        ' "$input_file" >> "$output_file"
    }

    ls -1 | while read i; do convert_to_bedgraph "$i" "$i.bedGraph"; done
    ls -1 *.txt.bedGraph|while read i; do bedGraphToBigWig $i /path/to/reference/size/S_lycopersicum_chromosomes.4.00.chrom.sizes $i.bw; done

    # MethPipe软件：计算甲基化区域并生成CSV表格
    cd ../..)

    # 给bam文件排序
    for i in `ls -1 *.bam`; do samtools sort $i -o $i.sorted.bam -@ 28; done

    # 调用MethPipe软件计算甲基化水平
    ls -1  *.sorted.bam| while read i; do methcounts -c /path/to/reference/ref/S_lycopersicum_chromosomes.4.00.fa $i > $i.meth; done
    

    # 拆分文件
    ls -1 *.meth | while read i; do
        awk -v filename="$i" '{
            if ($4 == "CpG" || $4 == "CpGx") print > (filename ".CpG.1.meth");
            else if ($4 == "CXG" || $4 == "CCG" || $4 == "CCGx" || $4 == "CXGx") print > (filename ".CHG.1.meth");
            else if ($4 == "CHH" || $4 == "CHHx") print > (filename ".CHH.1.meth");
        }' "$i"
    done

    # 寻找高甲基化区域并转化为CSV
    ls -1 *.1.meth|while read i; do echo "hypermr -o $i.hmr $i" >>hmr.sh; done
    parallel -j 6 -a hmr.sh
    ls -1 *.hmr| while read i; do echo "sed 's/\t/,/g' $i > $i.csv">>sed.sh; done
    parallel -j 6 -a sed.sh
    rm -rf *.hmr *.meth

    cd "$initial_dir"
}

# 主函数
main() {
    
    # 遍历所有输入的文件夹
    for folder in "$@"; do
        if [ -d "$folder" ]; then
            process_methylation_data "$folder"
        else
            echo "Warning: $folder is not a directory, skipping."
        fi
    done
}

# 运行主函数
main "$@"
