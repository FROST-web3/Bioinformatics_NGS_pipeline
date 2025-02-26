#####对于甲基化数据画TSS图的额外处理(bisulfite-seq环境)

#改txt为bdg
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
                printf "%s\t%d\t%d\t%.4f\n", chrom, pos-1, pos, level
            }
        }
        ' "$input_file" >> "$output_file"

ls -1 *.bedGraph|while read i; do bedGraphToBigWig $i /ifs1/User/mahaifeng/my_work/bisulfite-seq/size/S_lycopersicum_chromosomes.4.00.chrom.sizes $i.bw; done






#进入deeptools 虚拟环境
ls -1 *.bw | while read i; do
    # 使用 scale-regions 替代 reference-point，这样可以分析整个基因区域
    computeMatrix scale-regions \
        --regionBodyLength 1000 \
        --beforeRegionStartLength 5000 \
        --afterRegionStartLength 3000 \
        --numberOfProcessors max/2 \
        --skipZeros \
        -R /ifs1/User/mahaifeng/my_work/bisulfite-seq/bed/gene.bed \
        -S $i \
        --missingDataAsZero \
        -o $i.gene_body.gz \
        --outFileSortedRegions $i.gene_body.bed

    # 绘制热图，使用蓝到红的颜色方案
    plotHeatmap \
        -m $i.gene_body.gz \
        -out $i.gene_body_Heatmap.png \
        --colorMap 'RdBu_r' \
        --whatToShow 'plot, heatmap and colorbar' \
        --heatmapHeight 15 \
        --heatmapWidth 4 \
        --legendLocation upper-right
done


#基因起始位点上游5K，终止位点下游3K
#基因体1K
#配色反转


