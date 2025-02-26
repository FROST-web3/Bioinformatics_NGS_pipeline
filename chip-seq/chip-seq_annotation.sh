##!/bin/bash
#提前进入bisulfite-seq虚拟环境
#提前构建homer的基因组库
#组蛋白用*.broadPeak，转录因子用*.narrowPeak

shopt -s extglob

annotation(){
    local input_dir=$1
    local initial_dir=$(pwd)
    if [ ! -d "$input_dir" ]; then
        echo "Error: Directory $input_dir does not exist\n\n\n"
        return 1
    fi
    cd "$input_dir"
    ls -1 *.narrowPeak|while read i; do 
        annotatePeaks.pl $i tomato > ${i%.narrowPeak}_annotated.txt 
    done

    cd "$initial_dir"

}


main(){
    for folder in "$@"; do
        annotation "$folder"
    done
}

main "$@"