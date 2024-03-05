library(monocle)

HSMM_HNSCC_epi_tumor <- readRDS('/work/smy/Project/HNSCC_26sample/2.data/26sample_epi_tumor_monocleFigure.rds')
plot_cell_trajectory(HSMM_HNSCC_epi_tumor, color_by = "stage")
pdf('/work/brj/Collaboration/2022/scRNA/HNSCC/Result/Monocle/Monocle_Stage.pdf',4,4)
plot_cell_trajectory(HSMM_HNSCC_epi_tumor, 
                     theta = -10,
                     show_branch_points = F,
                     show_tree = TRUE, color_by = "stage", cell_size = 0.8) +
  scale_color_manual(values = .cluster_cols)+ NoLegend()
  #theme(legend.position = "bottom")
dev.off()

pdf('/work/brj/Collaboration/2022/scRNA/HNSCC/Result/Monocle/Monocle_Split.pdf',32,4)
plot_cell_trajectory(HSMM_HNSCC_epi_tumor, 
                     theta = -10,
                     show_branch_points = F,
                     show_tree = TRUE, color_by = "stage", cell_size = 0.8) +
  scale_color_manual(values = .cluster_cols)+
  NoLegend() + facet_wrap('~stage',nrow = 1)
dev.off()

sig_gene_names1 <- c("ELF3","GRHL1","MXD1","EHF","FOXA1","FOXM1","MYBL2","YBX1","TFDP1","MAZ",
                     "TP63","IRF6","KDM5B","BHLHE40","TFAP2A","ID1","PITX1","NR4A1","CREB3L1","SPDEF","XBP1","FOSB")

sig_gene_names2 <- c("CEBPD","BCL3","GRHL1","MXD1","RORC","MYBL2","POLE3","BRCA1","TFDP1",
                     "STAT1","TP63","KLF7","STAT1","HIF1A","MYC","GRHL3","RARG",
                     "SOX15","GRHL1","MAF","BATF","RUNX3","FLI1","IRF4","FOXP1")

## 使用BEAM函数进行分支表达建模分析 BEAM：branched expression analysis modeling
BEAM_res <- BEAM(HSMM_HNSCC_epi_tumor, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

readr::write_rds(BEAM_res,"/work/brj/Collaboration/2022/scRNA/HNSCC/Result/Monocle/BEAM_res.rds.gz")

BEAM_res <- readr::read_rds("/work/brj/Collaboration/2022/scRNA/HNSCC/Result/Monocle/BEAM_res.rds.gz")
pdf('/work/brj/Collaboration/2022/scRNA/HNSCC/Result/Monocle/branched_heatmap.pdf',width = 4.5, height = 15)
plot_genes_branched_heatmap(HSMM_HNSCC_epi_tumor[row.names(subset(BEAM_res,qval < 1e-4)),],
                            cluster_rows = F,
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()

pdf('/work/brj/Collaboration/2022/scRNA/HNSCC/Result/Monocle/branched_heatmapsig1.pdf',width = 4.5, height = 15)
plot_genes_branched_heatmap(HSMM_HNSCC_epi_tumor[row.names(subset(BEAM_res,gene_short_name %in% sig_gene_names1)),],
                            branch_point = 1,
                            cores = 1,
                            use_gene_short_name = T,
                            cluster_rows = F,
                            show_rownames = T)
dev.off()

pdf('/work/brj/Collaboration/2022/scRNA/HNSCC/Result/Monocle/branched_heatmapsig2.pdf',width = 4.5, height = 15)
plot_genes_branched_heatmap(HSMM_HNSCC_epi_tumor[row.names(subset(BEAM_res,gene_short_name %in% unique(sig_gene_names2))),],
                            cluster_rows = FALSE,
                            #branch_point = 1,
                            #branch_states = c(2,3),
                            #num_clusters = 1,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()

#https://alexthiery.github.io/otic-reprogramming/downstream/smartseq2_downstream/
pdf('/work/brj/Collaboration/2022/scRNA/HNSCC/Result/Monocle/branched_heatmapsig1clusters1.pdf',width = 4, height = 4)
beam_hm = plot_genes_branched_heatmap(HSMM_HNSCC_epi_tumor[sig_gene_names1,],
                                      branch_point = 1,
                                      cluster_rows=FALSE,
                                      num_clusters = 1,
                                      cores = 1,
                                      use_gene_short_name = T,
                                      show_rownames = T,
                                      return_heatmap=T,
                                      branch_colors=c("#66C2A5","#FC8D62","#8DA0CB"),
                                      branch_labels=c("S2",'S3'))
dev.off()
pdf('/work/brj/Collaboration/2022/scRNA/HNSCC/Result/Monocle/branched_heatmapsig2clusters1.pdf',width = 4, height = 4)
beam_hm = plot_genes_branched_heatmap(HSMM_HNSCC_epi_tumor[unique(sig_gene_names2),],
                                      branch_point = 1,
                                      cluster_rows=FALSE,
                                      num_clusters = 1,
                                      cores = 1,
                                      use_gene_short_name = T,
                                      show_rownames = T,
                                      return_heatmap=T,
                                      branch_colors=c("#66C2A5","#FC8D62","#8DA0CB"),
                                      branch_labels=c("S2",'S3'))
dev.off()


pdf('/work/brj/Collaboration/2022/scRNA/HNSCC/Result/Monocle/branched_heatmapsig1cluster3.pdf',width = 4, height = 4)
beam_hm = plot_genes_branched_heatmap(HSMM_HNSCC_epi_tumor[sig_gene_names1,],
                                      branch_point = 1,
                                      cluster_rows=FALSE,
                                      num_clusters = 3,
                                      cores = 1,
                                      use_gene_short_name = T,
                                      show_rownames = T,
                                      return_heatmap=T,
                                      branch_colors=c("#66C2A5","#FC8D62","#8DA0CB"),
                                      branch_labels=c("S2",'S3'))
dev.off()
pdf('/work/brj/Collaboration/2022/scRNA/HNSCC/Result/Monocle/branched_heatmapsig2cluster3.pdf',width = 4, height = 4)
beam_hm = plot_genes_branched_heatmap(HSMM_HNSCC_epi_tumor[unique(sig_gene_names2),],
                                      branch_point = 1,
                                      cluster_rows=FALSE,
                                      num_clusters = 3,
                                      cores = 1,
                                      use_gene_short_name = T,
                                      show_rownames = T,
                                      return_heatmap=T,
                                      branch_colors=c("#66C2A5","#FC8D62","#8DA0CB"),
                                      branch_labels=c("S2",'S3'))
dev.off()
