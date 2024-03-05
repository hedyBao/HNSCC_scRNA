# Epithelial features and function enrichment and regulon construction for each malignant cell type
HNSCC_epi_tumor <- readRDS('/work/smy/Project/HNSCC/2.data/26sample_tumorcells_Harmony_immuneReference_All.rds')
Idents(HNSCC_epi_tumor) <- HNSCC_epi_tumor$stage
p1 <- DimPlot(HNSCC_epi_tumor,cols = .cluster_cols)
Idents(HNSCC_epi_tumor) <- HNSCC_epi_tumor$copykat
p2 <- DimPlot(HNSCC_epi_tumor,cols = RColorBrewer::brewer.pal(10,'Dark2')[c(8,7)])
Idents(HNSCC_epi_tumor) <- HNSCC_epi_tumor$RNA_snn_res_final
p3 <- DimPlot(HNSCC_epi_tumor,cols = RColorBrewer::brewer.pal(10,'Dark2'))
### Fig2D
TumorSub_proportion <- read.delim('/work/smy/Project/HNSCC_26sample/2.data/CIBERSORTx_Results_TumorSubcluster.txt')
Tumorsub <- TumorSub_proportion %>% dplyr::select(c(Mixture,cluster_0,cluster_1,cluster_2,cluster_3,cluster_4))
HNSCC_proportion <- Tumorsub
HNSCC_proportion <- HNSCC_proportion[substr(HNSCC_proportion$Mixture,14,16) == "01A",]
ClinicalData <- read.delim("/server1_work/brj/GEO/CRC/multiCox/TCGA_ClinicalData_20180420.txt")
ClinicalData$OS <- as.numeric(as.character(ClinicalData$OS))
ClinicalData$OS.time <- as.numeric(as.character(ClinicalData$OS.time)) 
HNSCC_proportion$Mixture <- substr(HNSCC_proportion$Mixture,1,12)
HNSCC_proportion$Mixture <- gsub('[.]','-',HNSCC_proportion$Mixture)
HNSCC_proportion[,c(7,8)] <- ClinicalData[match(HNSCC_proportion$Mixture,ClinicalData$bcr_patient_barcode),c('OS','OS.time')] 
HNSCC_proportion <- HNSCC_proportion %>% dplyr::mutate(OS.time = OS.time/30) %>% dplyr::filter(!is.na(OS.time))
openxlsx::write.xlsx(HNSCC_proportion,file=paste(OutPath,"/",folder,"/Fig2D.xlsx",sep=""))
### Fig2E
source("/server1_work/yye/AliRstudio/Data/Public/ToolsData/ulcerative_colitis/scores.r")
HNSCC_epi_tumor <- readRDS('/work/smy/Project/HNSCC_26sample/2.data/26sample_tumorcells_Harmony_immuneReference_All.rds')
HMterms <- gmtPathways("/work/smy/Public/data/GeneList/h.all.v6.1.symbols.gmt") 

HMscores <- score_cells(seur=HNSCC_epi_tumor, names=HMterms, combine_genes='mean', groups=NULL, group_stat='mean', cells.use=NULL)
HMscoresMatrix <- as.matrix(HMscores)
HMscoresMatrix  <- t(HMscoresMatrix )
HMscoresMatrix <- apply(HMscoresMatrix,2,function(x)signif(x,digits = 3))
#rownames(HMscoresMatrix ) <- names(nHMPathways)
colnames(HMscoresMatrix) <- rownames(HNSCC_epi_tumor@meta.data)
HMscoresMatrixSeurat <- CreateSeuratObject(counts = HMscoresMatrix )
HMscoresMatrixSeurat@meta.data <- HNSCC_epi_tumor@meta.data
HNSCC_epi_tumor@assays$HM <- HMscoresMatrixSeurat@assays$RNA
DefaultAssay(HNSCC_epi_tumor) <- "HM"
Idents(HNSCC_epi_tumor) <- "RNA_snn_res_final"
HNSCC_epi_tumor_HM_FeatureDiff <- FindAllMarkers(HNSCC_epi_tumor,min.pct=0.1,logfc.threshold = 0.05,pseudocount.use = 0.1,only.pos = T)
topN <- 5
HNSCC_epi_tumor_HM_FeatureDiff$cluster <- as.character(HNSCC_epi_tumor_HM_FeatureDiff$cluster)
HNSCC_epi_tumor_HM_FeatureDiff_Order <- HNSCC_epi_tumor_HM_FeatureDiff %>% dplyr::arrange(cluster,desc(avg_log2FC)) %>% dplyr::group_by(cluster) %>% dplyr::top_n(topN,avg_log2FC) # 
HNSCC_epi_tumor_HM_FeatureDiff_Order$gene <- gsub("[.]"," ",HNSCC_epi_tumor_HM_FeatureDiff_Order$gene)
HNSCC_epi_tumor_HM_FeatureDiff_Order$gene <- gsub("   "," ",HNSCC_epi_tumor_HM_FeatureDiff_Order$gene)
TT <- data.frame(Subtypes=rep(unique(HNSCC_epi_tumor_HM_FeatureDiff_Order$cluster),times=length(unique(HNSCC_epi_tumor_HM_FeatureDiff_Order$gene))),
                 gene=rep(unique(HNSCC_epi_tumor_HM_FeatureDiff_Order$gene),each=length(unique(HNSCC_epi_tumor_HM_FeatureDiff_Order$cluster))))
#图片处理 美观
HNSCC_epi_tumor_HM_FeatureDiff_OrderF <- HNSCC_epi_tumor_HM_FeatureDiff_Order %>% dplyr::select(gene,cluster,avg_log2FC) %>%
  tidyr::spread(cluster,avg_log2FC,fill=0)
HNSCC_epi_tumor_HM_FeatureDiff_OrderF <- HNSCC_epi_tumor_HM_FeatureDiff_OrderF[do.call(order,HNSCC_epi_tumor_HM_FeatureDiff_OrderF[,2:ncol(HNSCC_epi_tumor_HM_FeatureDiff_OrderF)]),]
HNSCC_epi_tumor_HM_FeatureDiff_Order$cluster <- paste0("C",HNSCC_epi_tumor_HM_FeatureDiff_Order$cluster)
openxlsx::write.xlsx(HNSCC_epi_tumor_HM_FeatureDiff_Order,file=paste(OutPath,"/",folder,"/Fig2E.xlsx",sep=""))

### Fig2F
MBterms <- gmtPathways("/work/smy/Public/data/GeneList/KEGG_metabolism.gmt") 
MBscores <- score_cells(seur=HNSCC_epi_tumor, names=MBterms, combine_genes='mean', groups=NULL, group_stat='mean', cells.use=NULL)
MBscoresMatrix <- as.matrix(MBscores)
MBscoresMatrix  <- t(MBscoresMatrix )
MBscoresMatrix <- apply(MBscoresMatrix,2,function(x)signif(x,digits = 3))
#rownames(MBscoresMatrix ) <- names(nMBPathways)
colnames(MBscoresMatrix) <- rownames(HNSCC_epi_tumor@meta.data)
MBscoresMatrixSeurat <- CreateSeuratObject(counts = MBscoresMatrix )
MBscoresMatrixSeurat@meta.data <- HNSCC_epi_tumor@meta.data
HNSCC_epi_tumor@assays$MB <- MBscoresMatrixSeurat@assays$RNA
DefaultAssay(HNSCC_epi_tumor) <- "MB"
Idents(HNSCC_epi_tumor) <- "RNA_snn_res_final"
HNSCC_epi_tumor_MB_FeatureDiff <- FindAllMarkers(HNSCC_epi_tumor,min.pct=0.1,logfc.threshold = 0.05,pseudocount.use = 0.1,only.pos = T)
topN <- 6
HNSCC_epi_tumor_MB_FeatureDiff$cluster <- as.character(HNSCC_epi_tumor_MB_FeatureDiff$cluster)
HNSCC_epi_tumor_MB_FeatureDiff_Order <- HNSCC_epi_tumor_MB_FeatureDiff %>% dplyr::arrange(cluster,desc(avg_log2FC)) %>% dplyr::group_by(cluster) %>% dplyr::top_n(topN,avg_log2FC) # 
HNSCC_epi_tumor_MB_FeatureDiff_Order$gene <- gsub("[.]"," ",HNSCC_epi_tumor_MB_FeatureDiff_Order$gene)
HNSCC_epi_tumor_MB_FeatureDiff_Order$gene <- gsub("   "," ",HNSCC_epi_tumor_MB_FeatureDiff_Order$gene)
TT <- data.frame(Subtypes=rep(unique(HNSCC_epi_tumor_MB_FeatureDiff_Order$cluster),times=length(unique(HNSCC_epi_tumor_MB_FeatureDiff_Order$gene))),
                 gene=rep(unique(HNSCC_epi_tumor_MB_FeatureDiff_Order$gene),each=length(unique(HNSCC_epi_tumor_MB_FeatureDiff_Order$cluster))))
HNSCC_epi_tumor_MB_FeatureDiff_OrderF <- HNSCC_epi_tumor_MB_FeatureDiff_Order %>% dplyr::select(gene,cluster,avg_log2FC) %>%
  tidyr::spread(cluster,avg_log2FC,fill=0)
HNSCC_epi_tumor_MB_FeatureDiff_OrderF <- HNSCC_epi_tumor_MB_FeatureDiff_OrderF[do.call(order,HNSCC_epi_tumor_MB_FeatureDiff_OrderF[,2:ncol(HNSCC_epi_tumor_MB_FeatureDiff_OrderF)]),]

ggplot(HNSCC_epi_tumor_MB_FeatureDiff_Order,aes(x=gene,y=cluster))+
  geom_point(aes(color=avg_log2FC,size=-log10(p_val_adj+10^(-100))))+
  scale_color_gradientn(colors= rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"),space="rgb")(50)))+
  theme(panel.background=element_rect(color="black",fill=NA),legend.key =  element_blank(),axis.ticks=element_blank(),
        panel.grid=element_blank(),axis.text.x=element_text(size=8,angle=90,vjust = 0.5,hjust = 1,color="black"),
        axis.text.y=element_text(size=8,color="black"),axis.title=element_blank(),legend.position = "bottom")+
  scale_x_discrete(limit=rev(HNSCC_epi_tumor_MB_FeatureDiff_OrderF$gene))+
  #scale_y_discrete(limit = c('0','1','2'))+
  coord_flip()+
  geom_tile(data=TT,aes(y=Subtypes,x=gene),fill=NA,color="lightgray")

HNSCC_epi_tumor_MB_FeatureDiff_Order$cluster <- paste0("C",HNSCC_epi_tumor_MB_FeatureDiff_Order$cluster)
openxlsx::write.xlsx(HNSCC_epi_tumor_MB_FeatureDiff_Order,file=paste(OutPath,"/",folder,"/Fig2F.xlsx",sep=""))

### Fig2H left
### server1:/work/smy/Project/HNSCC/1.code/13.26Sample_tumorcell_SCENIC.R
Tumor_Anno <- readRDS('/server1_work/smy/Project/HNSCC/2.data/SCENIC/TumorAnno_AddTFs.rds.gz')
Tumor_Anno@meta.data$Resolution_final <- HNSCC_epi_tumor@meta.data[match(rownames(Tumor_Anno@meta.data),rownames(HNSCC_epi_tumor@meta.data)),'RNA_snn_res_final']
Sub_Data <- Tumor_Anno

CandidateFeature <- openxlsx::read.xlsx("/server1_work/smy/Project/HNSCC/2.data/SCENIC/Tumor_Regulon_Will.xlsx")
CandidateFeature  <- CandidateFeature[-grep('extended',CandidateFeature$gene),]
CandidateFeature$cluster <- factor(CandidateFeature$cluster,levels=unique(CandidateFeature$cluster))

CandidateFeature <- CandidateFeature[(order(CandidateFeature$cluster,sort(CandidateFeature$p_val))),]
CandidateFeature <- CandidateFeature %>% dplyr::group_by(cluster) %>% dplyr::arrange(cluster,p_val,desc(avg_logFC)) 
CandidateFeature <- CandidateFeature %>% dplyr::top_n(5,rev(p_val)) #select top10 TF regulons

CandidateFeatureList <- lapply(levels(CandidateFeature$cluster),function(cluster){
  ls1 <- CandidateFeature[CandidateFeature$cluster == cluster,]$gene
  ls2 <- CandidateFeature[CandidateFeature$cluster != cluster,]$gene
  ls1_unique <- ls1[!ls1 %in% ls2]
  if (length(ls1_unique) > 10){
    sub <- CandidateFeature[CandidateFeature$cluster == cluster & CandidateFeature$gene %in% ls1_unique,]
    ls1_unique <- sort(sub$p_val,decreasing = F)[1:10]
  }else(ls1_unique <- ls1_unique)
  return(ls1_unique)
}) 
names(CandidateFeatureList) <- levels(CandidateFeature$cluster)

PlotGenes <- CandidateFeature %>% dplyr::top_n(10,rev(p_val))
PlotGenesF <- PlotGenes$gene
PlotGenes_Can <- CandidateFeature %>% dplyr::group_by(cluster) %>% dplyr::top_n(10,rev(p_val))
#PlotGenes_CanF <- PlotGenes_Can[PlotGenes_Can$cluster %in% "FAP+ Tumorroblasts",]$gene

#get TF matrix
SubDatam <- Sub_Data@assays$TF@data%>% as.matrix
colnames(SubDatam) <- gsub("\\.","-",colnames(SubDatam))
SubDatam <- SubDatam[apply(SubDatam,1,sum) > 0,]
ClusterMean <- t(apply(SubDatam,1,function(x){
  lapply(Sub_ClusterID,function(ID){
    mean(x[ID]) 
  }) %>% unlist
}))
ClusterMean <- ClusterMean[,levels(CandidateFeature$cluster)]
ClusterMean <- ClusterMean[apply(ClusterMean,1,sum)>0,]

ClusterMeanM <- ClusterMean[match(CandidateFeature$gene,rownames(ClusterMean)),]
ClusterMeanM[,1:ncol(ClusterMeanM)] <- t(apply(ClusterMeanM,1,scale))

PlotDF <- ClusterMeanM %>% as.data.frame()

openxlsx::write.xlsx(PlotDF,file=paste(OutPath,"/","Fig2","/Fig2Hleft.xlsx",sep=""),row.names = T)

### right
TF_candidate <- read.table('/work/smy/Project/HNSCC_26sample/2.data/SCENIC_Tumorcell_TF_filter_MeanExpr.txt')
names(TF_candidate) <- gsub("C","X",names(TF_candidate))
openxlsx::write.xlsx(TF_candidate,file=paste(OutPath,"/",folder,"/Fig2Hright.xlsx",sep=""),row.names = T)
### Fig2H right
Sub_Data <- HNSCC_epi_tumor
#####plot TF relative heatmap
sig_gene_names1 <- c("CEBPD","BCL3","GRHL1",
                     "MXD1","RORC","MYBL2","POLE3",
                     "BRCA1","TFDP1","STAT1","TP63","KLF7",
                     "EHF","FOXA1","FOXM1","MYBL2","YBX1","TFDP1","MAZ",
                     "TP63","IRF6","KDM5B","BHLHE40","TFAP2A","ID1","PITX1","NR4A1","CREB3L1","SPDEF","XBP1","FOSB")
PlotGenes <- sig_gene_names1
Sub_TF_Matrix <- Sub_Data@assays$RNA@data[PlotGenes,] %>% as.matrix()
Sub_TF_Matrix <- Sub_TF_Matrix[apply(Sub_TF_Matrix,1,sum) > 0,]

Sub_ClusterID <- lapply(split(Sub_Data@meta.data,list(Sub_Data$RNA_snn_res_final)), function(x)rownames(x))

Sub_TF_ClusterMean <- t(apply(Sub_TF_Matrix,1,function(x){
  lapply(Sub_ClusterID,function(ID){
    mean(x[ID])
  }) %>% unlist
}))
Sub_TF_ClusterMean <- Sub_TF_ClusterMean[apply(Sub_TF_ClusterMean,1,sum)>0,]
Sub_TF_ClusterMean[,1:ncol(Sub_TF_ClusterMean)] <- t(apply(Sub_TF_ClusterMean,1,scale)) 
PlotDF <- Sub_TF_ClusterMean %>% as.data.frame()

openxlsx::write.xlsx(PlotDF,file=paste(OutPath,"/","Fig2","/Fig2Hright.xlsx",sep=""))

##2K : SECNIC TF survival=======
TF_candidate <- read.table('/work/smy/Project/HNSCC_26sample/2.data/SCENIC_Tumorcell_TF_filter_MeanExpr.txt')
rownames(TF_candidate)
HNSC_FPKM <- readr::read_rds("/work/smy/Public/data/RNA-seq/TCGA/HNSCC_FPKM.rds")
ClinicalData <- read.delim("/work/smy/Public/data/RNA-seq/TCGA/TCGA_ClinicalData_20180420.txt")
ClinicalData$PFI.time <- as.numeric(as.character(ClinicalData$PFI.time))
ClinicalData$OS.time <- as.numeric(as.character(ClinicalData$OS.time))
ClinicalDataF <- ClinicalData[,c("bcr_patient_barcode",'OS','OS.time')]
HNSC_FPKM_F <- HNSC_FPKM[rownames(TF_candidate),]
HNSC_FPKM_F <- t(HNSC_FPKM_F) %>% as.data.frame() 
HNSC_FPKM_F[,1:22] <- apply(HNSC_FPKM_F[,1:22],2,function(x)log2(x+1))

HNSC_FPKM_F$bcr_patient_barcode <- substr(rownames(HNSC_FPKM_F),0,12)
HNSC_FPKM_FF <- HNSC_FPKM_F[substr(rownames(HNSC_FPKM_F),14,16) %in% c("01A"),]

HNSC_FPKM_clinical <- merge(HNSC_FPKM_FF,ClinicalDataF,by = "bcr_patient_barcode") %>% 
  dplyr::filter(!is.na(OS.time)) %>% 
  dplyr::mutate(OS.time = OS.time/30) %>%
  dplyr::select(c(bcr_patient_barcode,TFDP1,OS,OS.time))
openxlsx::write.xlsx(HNSC_FPKM_clinical,file=paste(OutPath,"/",folder,"/Fig2K.xlsx",sep=""))
