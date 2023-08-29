rmsk_gtf=/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/velocity/ref/GRCh38_rmsk.gtf # 从genome.ucsc.edu下载 
cellranger_outDir=/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/cellranger/counts #前面cellranger命令的outputs目录 
cellranger_gtf=/cluster/home/yzy_jh/reference/scell/cellranger/human/refdata-gex-GRCh38-2020-A/genes/genes.gtf # 这个是cellranger官网提供的
#ls -lh $rmsk_gtf  $cellranger_outDir $cellranger_gtfi
ls  $cellranger_outDir/*/outs/possorted_genome_bam.bam|while read id;do  new=${id/possorted_genome_bam.bam/cellsorted_possorted_genome_bam.bam}
echo $new 
nohup samtools sort -@ 4  -t CB -O BAM -o $new   $id  & 
done

ls -d *|while read cellranger_outDir;do 
nohup velocyto run10x -m $rmsk_gtf  $cellranger_outDir $cellranger_gtf & 
done 

