#!/bin/bash
shopt -s extglob

hisat2() {
    local input_dir=$1
    local initial_dir=$(pwd)
    
    cd "$input_dir"
    ls -1 *.gz| xargs -n2 | while read r1 r2; do 
    sample_name=${r1%.fastq*}
    echo "hisat2 -x /ifs1/User/FROST/my_work/rna-seq/ref/S_lycopersicum_chromosomes.4.00.fa -1 ${r1} -2 ${r2} -S ${sample_name}.sam -p 28;
    samtools sort -@ 28 "${sample_name}.sam" -o "${sample_name}.sorted.bam";
    samtools index -@ 28 "${sample_name}.sorted.bam";
    rm -f "${sample_name}.sam";">> hisat.sh;

    done;
    sh hisat.sh;
    featureCounts -g gene_id -a /ifs1/User/mahaifeng/my_work/rna-seq/ref/ITAG4.0_gene_models.gtf --primary -T 28 -p -o CountMatrix.txt *.bam
    

    cd "$initial_dir"
}

main() {
    for folder in "$@"; do
        if [ -d "$folder" ]; then
            hisat2 "$folder"
        else
            echo "Error: $folder is not a directory"
            exit 1
        fi
    done
    
}

main "$@"

#featureCounts -g gene_id -a /ifs1/User/mahaifeng/my_work/rna-seq/ref/ITAG4.0_gene_models.gtf --primary -T 28 -p -o CountMatrix.txt *.bam
