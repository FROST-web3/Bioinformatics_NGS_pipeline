#进入deeptools 虚拟环境
ls -1 *.bw | while read i; do
    # 使用 scale-regions 替代 reference-point，这样可以分析整个基因区域
    computeMatrix scale-regions \
        --regionBodyLength 1000 \
        --beforeRegionStartLength 5000 \
        --afterRegionStartLength 3000 \
        --numberOfProcessors 10 \
        --skipZeros \
        -R /ifs1/User/mahaifeng/my_work/bisulfite-seq/bed/gene.bed \
        -S $i \
        --missingDataAsZero \
        -o $i.gene_body.gz \
        --outFileNameMatrix matrix.tab # 输出文本格式的矩阵,用于后续处理

    # 绘制热图，使用蓝到红的颜色方案
    plotHeatmap \
        -m $i.gene_body.gz \
        -out $i.gene_body_Heatmap.png \
        --colorMap 'RdBu_r' \
        --legendLocation none \

done


#基因起始位点上游5K，终止位点下游3K
#基因体1K


