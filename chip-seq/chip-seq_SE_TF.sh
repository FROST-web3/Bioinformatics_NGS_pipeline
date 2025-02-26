#!/bin/bash
#提前进入bisulfite-seq虚拟环境，定位在sra文件夹所在的文件夹

shopt -s extglob

alignment_bowtie2() {
    local input_dir=$1
    local initial_dir=$(pwd)
    echo "start to do bowtie2_alignment in  $input_dir"
    cd "$input_dir"

    # 拆分sra，合并fastq,删除临时文件
    mv */* ./
    rm -rf */
    fasterq-dump -e 28 *.sra
    cat *.fastq > mix.fastq
    rm -rf S*.fastq

    # bowtie2比对
    ls -1 mi*.fastq | while read i; do
        echo "bowtie2 -x /ifs1/User/mahaifeng/my_work/chip-seq/ref/Solanum_lycopersicum.SL3.0.dna.toplevel -U $i -S $i.sam -p 28" >> bowtie2.sh
    done
    sh bowtie2.sh
    rm -rf *.fastq bowtie2.sh

    # sam2bam
    ls -1 *.sam | while read i; do
        echo "samtools view -@ 28 -b $i > $i.bam" >> sam2bam.sh
    done
    sh sam2bam.sh
    rm -rf sam2bam.sh 

    # sort bam
    ls -1 *.bam | while read i; do
        echo "samtools sort -@ 28 $i -o $i.sorted.bam" >> sort.sh
    done
    sh sort.sh
    rm -rf sort.sh 

    # deduplicated
    ls -1 *.sorted.bam | while read i; do
        echo "samtools markdup -@ 28  $i $i.markdup.bam" >> markdup-r.sh
    done
    sh markdup-r.sh
    rm -rf markdup-r.sh

    # 删除临时文件
    rm -rf !(*.markdup.bam|*.sra)

    # 返回初始文件夹
    cd "$initial_dir"
}

macs2_peak_calling() {
    local input_dir=$1
    local initial_dir=$(pwd)
    echo "start to do MACS2 peak calling in $input_dir "
    cd "$input_dir"

    # MACS2 peak calling for TF
    ls -1 *.markdup.bam | while read bam_file; do
    macs2 callpeak -t "$bam_file" -n "$bam_file.chip.out" --outdir ./ -g 828000000 -B --SPMR --keep-dup auto -q 0.05 --shift -75 --call-summits
    done
    

    #生成bw可视化文件
    ls -1 *.bdg|while read i ;do 
     sort -k1,1 -k2,2n "$i" > "$i.sorted"
     bedGraphToBigWig "$i.sorted" /ifs1/User/mahaifeng/my_work/bisulfite-seq/tomato/chromSize/chrom.sizes "$i.bw"
     rm -rf "$i.sorted"
    done



    # 返回初始文件夹
    cd "$initial_dir"
}




main() {
    for folder in "$@"; do
        if [ -d "$folder" ]; then
            alignment_bowtie2 "$folder"
            macs2_peak_calling "$folder"
        else
            echo "Warning: $folder is not a directory, skipping."
        fi
    done
}

main "$@"






