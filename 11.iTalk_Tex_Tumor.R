#devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)
library(nichenetr)
library(Seurat)
library(tidyverse)
library(circlize)
library(dplyr)
library(clusterProfiler)
library(RColorBrewer)
library(iTALK,lib.loc = "/home/yhdu/R/x86_64-pc-linux-gnu-library/4.1/")

# #load Seurat rds
HNSCC_Whole <- readRDS('/work/smy/Project/HNSCC/2.data/DeinfeTypes/HNSCC_26sampleAll_DefineTypes.rds.gz')
Idents(HNSCC_Whole) <- HNSCC_Whole$DefineTypes
HNSCC_Whole$celltype_RE <- HNSCC_Whole$DefineTypes

## add tumor in epi : metadata$Malignancy  == 'malignant'
metadata <- readRDS('/work/smy/Project/HNSCC_26sample/2.data/26sample_tumor_cell_metadata.rds')
malignant <- rownames(metadata[grep('malignant',metadata$Malignancy),])
HNSCC_Whole$celltype_RE[malignant] <- 'TumorCell'
HNSCC_TumorTex_Allstage <- subset(HNSCC_Whole,subset = celltype_RE %in% c("CD8 Tex","TumorCell"))
HNSCC_TumorTex <- subset(HNSCC_TumorTex_Allstage,subset = stage %in% c("LN-in","LN-out"))
Idents(HNSCC_TumorTex) <- HNSCC_TumorTex$celltype_RE
TumorTex_split <- SplitObject(HNSCC_TumorTex,split.by = 'stage')

top_genes <- "100"
# target genes = 50 , tumor to CD8 tex 只有 
TumorTex_iTALK <- lapply(names(TumorTex_split), function(name){
  merge_seu <- TumorTex_split[[name]]
  # iTALK 要求的矩阵: 行为细胞，列为基因
  iTALK_data <- as.data.frame(t(merge_seu@assays$RNA@counts))
  #iTALK_data需要加两列：cell_type和compare_group的metadata
  iTALK_data$cell_type <- merge_seu$celltype_RE
  iTALK_data$compare_group <- merge_seu$stage
  
  highly_exprs_genes <- rawParse(iTALK_data,top_genes = 100,stats = "mean")
  comm_list <- c("growth factor","other","cytokine","checkpoint")
  iTALK_res <- NULL
  for(i in comm_list){
    res <- FindLR(highly_exprs_genes, datatype="mean count", comm_type = i)
    iTALK_res <- rbind(iTALK_res,res)
  }
  outpath <- "/work/brj/Collaboration/2022/scRNA/HNSCC/Result/iTALK/"
  saveRDS(iTALK_res,paste(outpath,"/TextoTumor_",name,"iTALK_res.rds.gz",sep = ""))
  return(iTALK_res)
})
names(TumorTex_iTALK) <- names(TumorTex_split)

outpath <- "/work/brj/Collaboration/2022/scRNA/HNSCC/Result/iTALK/"
saveRDS(TumorTex_iTALK,paste(outpath,"/TexTumor_iTALK_res.rds.gz",sep = ""))

TumorTex_iTALK <- readr::read_rds(paste(outpath,"/TexTumor_iTALK_res.rds.gz",sep = ""))


iTALK_res_out <- TumorTex_iTALK[[1]]

checkpoint1 <- iTALK_res_out %>% dplyr::filter(cell_from=="TumorCell"&cell_to=="CD8 Tex"&comm_type=="checkpoint")
pdf(paste(outpath,"/TumortoTex_",top_genes,"iTALK_out_res.pdf",sep = ""),width = 5,height = 5)
LRPlot(checkpoint1,datatype='mean count',transparency=0.5,link.arr.lwd=(checkpoint1$cell_from_mean_exprs),link.arr.width=checkpoint1$cell_to_mean_exprs)
dev.off()


pdf(paste(outpath,"/TextoTumor_",top_genes,"iTALK_out_res.pdf",sep = ""),width = 5,height = 5)
LRPlot(checkpoint2,datatype='mean count',transparency=0.5,link.arr.lwd=(checkpoint2$cell_from_mean_exprs),link.arr.width=checkpoint2$cell_to_mean_exprs)
dev.off()

iTALK_res_in <- TumorTex_iTALK[[2]]

checkpoint3 <- iTALK_res_in %>% dplyr::filter(cell_from=="TumorCell"&cell_to=="CD8 Tex"&comm_type=="checkpoint")

pdf(paste(outpath,"/TumortoTex_",top_genes,"iTALK_in_res.pdf",sep = ""),width = 5,height = 5)
LRPlot(checkpoint3,datatype='mean count',transparency=0.5,link.arr.lwd=(checkpoint3$cell_from_mean_exprs),link.arr.width=checkpoint3$cell_to_mean_exprs)
dev.off()

checkpoint4 <- iTALK_res_in %>% dplyr::filter(cell_from=="CD8 Tex"&cell_to=="TumorCell"&comm_type=="checkpoint")
pdf(paste(outpath,"/TextoTumor_",top_genes,"iTALK_in_res.pdf",sep = ""),width = 5,height = 5)
LRPlot(checkpoint4,datatype='mean count',transparency=0.5,link.arr.lwd=(checkpoint4$cell_from_mean_exprs),link.arr.width=checkpoint4$cell_to_mean_exprs)
dev.off()

pdf(paste(outpath,"/TextoTumor_iTALK_res_targetn500.pdf",sep = ""),width = 5,height = 5)
cowplot::plot_grid(plotlist = list(LRPlot(checkpoint1,datatype='mean count',transparency=0.5,link.arr.lwd=(checkpoint1$cell_from_mean_exprs),link.arr.width=checkpoint1$cell_to_mean_exprs),
                                   LRPlot(checkpoint2,datatype='mean count',transparency=0.5,link.arr.lwd=(checkpoint2$cell_from_mean_exprs),link.arr.width=checkpoint2$cell_to_mean_exprs),
                                   LRPlot(checkpoint3,datatype='mean count',transparency=0.5,link.arr.lwd=(checkpoint3$cell_from_mean_exprs),link.arr.width=checkpoint3$cell_to_mean_exprs),
                                   LRPlot(checkpoint4,datatype='mean count',transparency=0.5,link.arr.lwd=(checkpoint4$cell_from_mean_exprs),link.arr.width=checkpoint4$cell_to_mean_exprs)))
dev.off()



### 计算组间  DEG
merge_seu1 <- HNSCC_TumorTex
merge_seu1[["RNA"]]@counts <- as.matrix(merge_seu1[["RNA"]]@counts)+1

# iTALK 要求的矩阵: 行为细胞，列为基因
iTalk_data <- as.data.frame(t(merge_seu1@assays$RNA@counts))
#iTALK_data需要加两列：cell_type和compare_group的metadata
iTalk_data$cell_type <- merge_seu1$celltype_RE
iTalk_data$compare_group <- merge_seu1$stage


deg_Tex <- DEG(iTalk_data %>% filter(cell_type=='CD8 Tex'),method='DESeq2',contrast=c('LN-out', 'LN-in'))
deg_TumorCell <- DEG(iTalk_data %>% filter(cell_type=='TumorCell'),method='DESeq2',contrast=c('LN-out', 'LN-in'))

res<-NULL
for(comm_type in comm_list){
  res_cat<-FindLR(deg_Tex,deg_TumorCell,datatype='DEG',comm_type=comm_type)
  res<-rbind(res,res_cat)
}

res_filter <- res[order(res$cell_from_logFC*res$cell_to_logFC,decreasing=T),]
res_filters <- res_filter %>% filter(comm_type %in% c("checkpoint"))

saveRDS(res,paste(outpath,"/TexTumor_DEGs_iTALK_res.rds.gz",sep = ""))


pdf(paste(outpath,"/TextoTumor_Checkpoint_DEG",top_genes,"iTALK_in_res.pdf",sep = ""),width = 5,height = 5)
LRPlot(res_filters,datatype='DEG',link.arr.lwd=res_filters$cell_from_logFC,link.arr.width=res_filters$cell_to_logFC)
dev.off()


### check iTALK DE L-R between out and in 
library(Seurat)
library(tidyverse)
library(dplyr)
setwd("/work/brj/Collaboration/2022/scRNA/HNSCC/Result/iTALK/")
# load data 
outpath <- "/work/brj/Collaboration/2022/scRNA/HNSCC/Result/iTALK/"
TexTumorDEG <- readr::read_rds(paste(outpath,"/TexTumor_DEGs_iTALK_res.rds.gz",sep = ""))
res_filters <- TexTumorDEG %>% filter(comm_type %in% c("checkpoint"))
TexLigand <- res_filters %>% filter(cell_from == "CD8 Tex") %>% pull(ligand) %>% unique()
TexReceptor <- res_filters %>% filter(cell_to == "CD8 Tex") %>% pull(receptor) %>% unique()
TumorLigand <- res_filters %>% filter(cell_from == "TumorCell") %>% pull(ligand) %>% unique()
TumorReceptor <- res_filters %>% filter(cell_to == "TumorCell") %>% pull(receptor) %>% unique()


DefaultAssay(HNSCC_TumorTex) <- "RNA"

Tex <- subset(HNSCC_TumorTex,subset = celltype_RE  == "CD8 Tex")
Tumor <- subset(HNSCC_TumorTex,subset = celltype_RE  == "TumorCell")

p1 <- VlnPlot(Tumor, features = c("TNFRSF14","CD274","LGALS9"),pt.size = 0.05, group.by = "stage",col = c("#882E72","#B17BA6"),ncol =1 )
p2 <- VlnPlot(Tex, features = c("BTLA","PDCD1","HAVCR2"),pt.size = 0.05, group.by = "stage",col =c("#882E72","#B17BA6"),ncol =1 )


pdf(paste0(outpath,"TexTumoriTalkLigandReceptorViolinPlot.pdf"),width = 6,height = 9)
cowplot::plot_grid(plotlist = list(p1,p2),ncol = 2)
dev.off()        
       





DefaultAssay(Tumor) <- "RNA"
vp_case1 <- function(gene_signature, file_name, test_sign){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(Tumor, features = signature,
            pt.size = 0.05, 
            group.by = "stage"
    ) + stat_compare_means(comparisons = test_sign, label = "p.signif")
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]])
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1) )
  }
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, "rplot.png")
  ggsave(file_name, width = 13, height = 8)
}

gene_sig <- TumorReceptor
comparisons <- list(c("LN-out","LN-in"))
vp_case1(gene_signature = gene_sig, file_name = "TumorReceptorLN-outLN-in", test_sign = comparisons)




### check iTALK DE L-R between out and in 
library(Seurat)
library(tidyverse)
library(dplyr)
setwd("/work/brj/Collaboration/2022/scRNA/HNSCC/Result/iTALK/")
# load data 
outpath <- "/work/brj/Collaboration/2022/scRNA/HNSCC/Result/iTALK/"
TexTumorDEG <- readr::read_rds(paste(outpath,"/TexTumor_DEGs_iTALK_res.rds.gz",sep = ""))
res_filters <- TexTumorDEG %>% filter(comm_type %in% c("checkpoint"))
TexLigand <- res_filters %>% filter(cell_from == "CD8 Tex") %>% pull(ligand) %>% unique()
TexReceptor <- res_filters %>% filter(cell_to == "CD8 Tex") %>% pull(receptor) %>% unique()
TumorLigand <- res_filters %>% filter(cell_from == "TumorCell") %>% pull(ligand) %>% unique()
TumorReceptor <- res_filters %>% filter(cell_to == "TumorCell") %>% pull(receptor) %>% unique()


Tex <- subset(HNSCC_TumorTex,subset = celltype_RE  == "CD8 Tex")
Tumor <- subset(HNSCC_TumorTex,subset = celltype_RE  == "TumorCell")

pdf(paste0(outpath,"TexLigandRecptorDotPlot.pdf"),width = 5,height = 3)
DotPlot(Tex,features = c(TexLigand[-c(1,4)],TexReceptor),group.by = "stage")+
        theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))+
        scale_color_distiller(palette = "RdBu") + NoLegend()   
dev.off()

pdf(paste0(outpath,"TumorLigandRecptorDotPlot.pdf"),width = 5,height = 3)
DotPlot(Tumor,features = c(TumorLigand[-c(3)],TumorReceptor),group.by = "stage")+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))+
  scale_color_distiller(palette = "RdBu") + NoLegend()   
dev.off()

VlnPlot(Tex , features = signature,
        pt.size = 0.05, 
        group.by = "stage", 
        y.max = 1.0
) + stat_compare_means(comparisons = test_sign, label = "p.signif")



DefaultAssay(Tumor) <- "RNA"
vp_case1 <- function(gene_signature, file_name, test_sign){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(Tex , features = signature,
            pt.size = 0.05, 
            group.by = "stage" #, 
            #y.max = 1.0
    ) + stat_compare_means(comparisons = test_sign, label = "p.signif")
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]])
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1) )
  }
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, "rplot.png")
  ggsave(file_name, width = 13, height = 8)
}

gene_sig <- TexReceptor
comparisons <- list(c("LN-out","LN-in"))
vp_case1(gene_signature = gene_sig, file_name = "TexReceptorLN-outLN-in", test_sign = comparisons)




### L-R pairs 
P1 <- VlnPlot(Tex , features = c("CXCL13"),pt.size = 0.05, group.by = "stage")
P2 <- VlnPlot(Tumor, features = c("ACKR4","CCR10","CXCR3","CXCR5","HTR2A","OPRD1"),pt.size = 0.05, group.by = "stage")


Tex <- subset(HNSCC_TumorTex,subset = celltype_RE %in% c("CD8 Tex"))
Tumor <- subset(HNSCC_TumorTex,subset = celltype_RE %in% c("TumorCell"))
P1 <- FeaturePlot(Tex , features = c("CXCL13"),split.by = "stage",ncol = 1)
P2 <- FeaturePlot(Tumor, features = c("CXCR5"),split.by = "stage",ncol = 1)

pdf(paste0(outpath,"Fig/","FeaturePlotCXCL13.pdf"),width = 12,height = 3)
cowplot::plot_grid(plotlist = list(P1,P2))    
dev.off()


