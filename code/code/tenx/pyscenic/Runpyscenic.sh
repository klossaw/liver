pyscenic grn \
--num_workers 20 \
--output adj.sample.tsv \
--method grnboost2 \
sample.loom \
hs_hgnc_tfs.txt #转录因子文件，1839  个基因的名字列表

pyscenic ctx \
adj.sample.tsv \
hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname sample.loom \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 30 \
--mask_dropouts

pyscenic aucell \
sample.loom \
reg.csv \
--output sample_SCENIC.loom \
--num_workers 30

