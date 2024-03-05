### malignant features
HNSCC_Whole <- readRDS('/work/smy/Project/HNSCC/2.data/DeinfeTypes/HNSCC_26sampleAll_DefineTypes.rds.gz')
Idents(HNSCC_Whole) <- 'Maincluster'
##Epi细胞=====
HNSCC_epi <- subset(HNSCC_Whole,idents = c('Epithelial cells')) ##22751
### Copykat 结果
copykat_epi <- readRDS('/work/smy/Project/HNSCC_26sample/2.data/copykatN4/26sample_copykat_epi_result.rds')
copykat_epi$prediction <- as.data.frame(copykat_epi$prediction)
HNSCC_epi$copykat <- copykat_epi$prediction[match(rownames(HNSCC_epi@meta.data),copykat_epi$prediction$cell.names),2]

PPercentage <- table(HNSCC_epi@meta.data[,c('orig.ident',"copykat")]) %>% data.frame %>%
  #Por：用一群细胞silent/expression占比(total:100%)
  dplyr::mutate(Por= Freq/apply(table(HNSCC_epi@meta.data[,c('orig.ident',"copykat")]),1,sum))


p <- ggplot(PPercentage,aes(x=orig.ident,y=Por*100,fill=copykat))+
  geom_bar(stat="identity",width = 0.7)+labs(x="",y="Percentage (%)")+
  scale_y_continuous(expand = c(0,0))+
  #scale_fill_manual(values=c('#85A3C2','#3A7699'),name="")+
  #geom_text(aes(y=Por*100,label=signif(Por*100,digits = 2)),  vjust=1.6, color="black", size=5)+labs(title='BCAT1')+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,size=10),legend.position = "top",panel.background = element_blank(),
        axis.line = element_line())

### SFig2B
HNSCC_epi_tumor <- readRDS('/work/smy/Project/HNSCC_26sample/2.data/26sample_tumorcells_Harmony_immuneReference_All.rds')
cellcyle <- readxl::read_xlsx('/work/smy/Public/data/GeneList/cellcyle_signature.xlsx')
colnames(cellcyle) <- cellcyle[1,]
cellcyle <- cellcyle[-1,]
cellcyle$`CELL CYCLE: G1/S` <- as.character(cellcyle$`CELL CYCLE: G1/S`)
cellcyle$`CELL CYCLE: G2/M` <- as.character(cellcyle$`CELL CYCLE: G2/M`)

DefaultAssay(HNSCC_epi_tumor) <- 'RNA'
HNSCC_epi_tumor <- CellCycleScoring(HNSCC_epi_tumor,s.features = cellcyle$`CELL CYCLE: G1/S`,g2m.features = cellcyle$`CELL CYCLE: G2/M`)

meta <- HNSCC_epi_tumor@meta.data[,c('orig.ident','Phase')]

HNSCC_epi_tumor$cycle <- ifelse(HNSCC_epi_tumor$Phase == 'G1','non-cycling cell','cycling cell')
PPercentage <- table(HNSCC_epi_tumor@meta.data[,c("stage","cycle")])%>% data.frame %>%
  #Por：用一群细胞silent/expression占比(total:100%)
  dplyr::mutate(Por= Freq/apply(table(HNSCC_epi_tumor@meta.data[,c("stage","cycle")]),1,sum))
PPercentage$cycle <- as.character(PPercentage$cycle)
openxlsx::write.xlsx(PPercentage,file=paste(OutPath,"/",folder,"/SFig2B.xlsx",sep=""))


##Sfig2c
HNSCC_Whole <- readRDS('/work/smy/Project/HNSCC/2.data/DeinfeTypes/HNSCC_26sampleAll_DefineTypes.rds.gz')
Idents(HNSCC_Whole) <- 'Maincluster'
##Epi细胞=====
HNSCC_epi <- subset(HNSCC_Whole,idents = c('Epithelial cells')) ##22751
### copykat result
copykat_epi <- readRDS('/work/smy/Project/HNSCC_26sample/2.data/copykatN4/26sample_copykat_epi_result.rds')
Prediction <- copykat_epi$prediction %>% as.data.frame()
### HNSCC_epi
SplDf <- HNSCC_epi@meta.data %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("cellid") %>% 
  dplyr::select(c(cellid,orig.ident,stage)) %>%
  inner_join(Prediction,by = c("cellid" = "cell.names")) %>%
  dplyr::filter(orig.ident %in% c("NT1","Pre-Ca1","E-OSCC1",
                                  "NT2","Pre-Ca2","E-OSCC2",
                                  "NT3","Pre-Ca3","E-OSCC3"))

CNA.test <- copykat_epi$CNAmat
expr_1 <- CNA.test[,-c(1:3)]
cnv_score <- as.data.frame(colSums(expr_1 * expr_1))
colnames(cnv_score)="cnv_score"
cnv_score <- tibble::rownames_to_column(cnv_score, var='cellid')
cnv_score$cellid <- gsub('[.]','-',cnv_score$cellid) #不包含参考的样本在内

spl <- c("NT3","Pre-Ca3","E-OSCC3")
CNV_mtx <- cnv_score %>%
  inner_join(SplDf) %>% 
  dplyr::filter(orig.ident %in% spl) %>%
  dplyr::mutate(group = paste(.$stage, .$copykat.pred,sep = "_"))
CNV_mtxs <- CNV_mtx %>% dplyr::filter(copykat.pred == "aneuploid")
sobjlists <- CNV_mtxs %>% dplyr::select(c(stage,cnv_score))
sobjlists$stage <- factor(sobjlists$stage,levels = c("NT","Pre","E"))
sobjlists <- sobjlists %>% arrange(stage)
openxlsx::write.xlsx(sobjlists,file=paste(OutPath,"/",folder,"/SFig2C.xlsx",sep=""))

##F
features <- c("NR4A1","CCL4","CXCL2","CXCL3","IL7R","CXCL10","CXCL14","IL1RN",
              "IL1R2","IL18","TYMP","CXCL9","TNFRSF12A","INHBA","TNC","PLAU","IL36G","SDC1",
              "ACKR3","EDN2","CXCL8","CXCL1","HTN3","EGFR","SAA2","SAA1","DEFB1","IL16",
              "IL4R","CD40","IL32","ANGPTL4","IL20RB","SEMA4B")
              Idents(HNSCC_epi_tumor) <- HNSCC_epi_tumor$stage
p <- DotPlot(HNSCC_epi_tumor,features = features,cols = 'RdBu') +
  theme(axis.text.x = element_text(size = 10,angle=90,hjust = 1)) +
  scale_y_discrete(limits = rev(c('NT','Pre','E','A','LN-in','LN-out','LN-normal','R')))+
  scale_colour_gradientn(colours =rev(RColorBrewer::brewer.pal(12,'RdBu')),limits = c(-2.5,2.5),name="")
PlotDf <- p$data       

openxlsx::write.xlsx(PlotDf,file=paste(OutPath,"/",folder,"/SFig2F.xlsx",sep=""),row.names = F)

##g 
PPercentage <- table(HNSCC_epi_tumor@meta.data[,c('RNA_snn_res_final',"stage")]) %>% data.frame %>%
  #Por：用一群细胞silent/expression占比(total:100%)
  dplyr::mutate(Por= Freq/apply(table(HNSCC_epi_tumor@meta.data[,c('RNA_snn_res_final',"stage")]),1,sum))

openxlsx::write.xlsx(PPercentage,file=paste(OutPath,"/",folder,"/SFig2G.xlsx",sep=""))

##j
library(stringr)
TumorSub_proportion <- read.delim('/work/smy/Project/HNSCC_26sample/2.data/CIBERSORTx_Results_TumorSubcluster.txt')
TumorSub_proportion$Mixture <- gsub('[.]','-',TumorSub_proportion$Mixture)
TumorSub_proportion <- TumorSub_proportion[str_sub(TumorSub_proportion$Mixture,14,16) %in% c("01A"),]
TumorSub_proportion$Mixture <- str_sub(TumorSub_proportion$Mixture,1,12)
TumorSub <- TumorSub_proportion %>% dplyr::select(c(Mixture,cluster_0,cluster_1,cluster_2,cluster_3,cluster_4))

##clinical info
surviva <- read.delim("/server1_work/brj/GEO/CRC/multiCox/TCGA_ClinicalData_20180420.txt")#,","gender"
survival <- surviva %>% dplyr::select("bcr_patient_barcode",'OS','OS.time') %>% 
  dplyr::rename("Mixture" = "bcr_patient_barcode") %>%
  inner_join(TumorSub) %>%
  dplyr::mutate(OS.time = as.numeric(as.numeric(as.character(.$OS.time))/30)) %>% 
  dplyr::filter(!is.na(OS.time))
openxlsx::write.xlsx(survival,file=paste(OutPath,"/",folder,"/SFig2j.xlsx",sep=""))


###k 
## /work/brj/Collaboration/2022/scRNA/HNSCC/Revise2/Code/6.TumorPuramSignatureScore.R
p1$data
openxlsx::write.xlsx(p1$data,file=paste(OutPath,"/",folder,"/SFig2k.xlsx",sep=""))

###l-m
HNSCC_epi_tumorF <- subset(HNSCC_epi_tumor,idents = c('1','2'))
Markers <- FindAllMarkers(HNSCC_epi_tumorF,only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

Markers$Class <- ifelse( Markers$avg_log2FC > 0.4 & Markers$p_val < 0.05 & Markers$cluster == '1',"Cluster_1",
                         ifelse(Markers$avg_log2FC > 0.4 & Markers$p_val < 0.05 & Markers$cluster == '2',"Cluster_2",'None'))
Markers$avg_log2FC <- ifelse(Markers$cluster == '2',-Markers$avg_log2FC,Markers$avg_log2FC)
Markers$pvalue_L <- -log10(Markers$p_val+10^-200)
sub <- Markers

Gene_Labels <- sub %>% dplyr::filter(Class %in% c("Cluster_1","Cluster_2") ) %>% dplyr::group_by(Class)  %>% 
               dplyr::top_n(15,abs(avg_log2FC))

GO_BP_Fun <- function(CompareGroup,groups=c("Cluster_1","Cluster_2"),OrgData="org.Hs.eg.db",pcut=0.05,qcut=0.2){
  library(clusterProfiler)
  if(OrgData=="org.Hs.eg.db"){library(org.Hs.eg.db)}else{library(org.Hs.eg.db)}
  CompareGroup <- CompareGroup %>% dplyr::mutate(GO_Compare=purrr::map(.x=CompareResult,function(.x){
    SigSub <- .x[.x$Class %in% groups,] %>% as.data.frame()
    if(length(unique(SigSub$gene))>80){
      formula_res <- clusterProfiler::compareCluster(
        keyType="SYMBOL",
        gene ~ Class,
        data=SigSub,
        fun="enrichGO", 
        OrgDb=OrgData,
        ont		   = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = pcut,
        qvalueCutoff  = qcut
      )
      
      # Run GO enrichment test and merge terms that are close to each other to remove result redundancy
      #   Simple_ego <- clusterProfiler::simplify(
      #     formula_res, 
      #     cutoff=0.7, 
      #     by="p.adjust", 
      #    select_fun=min
      #   )
      return(formula_res ) # 
    }else{return(NA)}    # 
  })) 
  
}

library(KEGG.db,lib.loc = '/home/xyding/R/x86_64-pc-linux-gnu-library/4.1/')
KeggPathway <- readRDS('/work/xyding/2022/metabolism/KEGG_hsapathway.rds.gz')

KEGG_Fun <- function(CompareGroup,groups=c("Cluster_1","Cluster_2"),OrgSp="hsa",pcut=0.05,qcut=0.2){
  CompareGroup <- CompareGroup %>% dplyr::mutate(KEGG_Compare= purrr::map(.x=CompareResult,function(.x){
    
    SigSub <- .x[.x$Class %in% c("Cluster_2","Cluster_1"),] %>% as.data.frame()
    genelist <- clusterProfiler::bitr(SigSub$gene,fromType = "SYMBOL",
                                      toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
    
    SigSub <- merge( SigSub,genelist,by.x="gene",by.y="SYMBOL")
    if(length(unique(SigSub$gene))>80){
      formula_res <- clusterProfiler::compareCluster(
        #keyType="SYMBOL",
        ENTREZID~Class, 
        data=SigSub, 
        organism=OrgSp,
        fun="enrichKEGG", 
        pAdjustMethod = "BH",
        pvalueCutoff  = pcut,
        qvalueCutoff  = qcut,
        use_internal_data = TRUE #设置使用内部database
        
      )
      formula_res@compareClusterResult$Description <- lapply(formula_res@compareClusterResult$ID,function(x)KeggPathway[KeggPathway$ID %in% x,]$Description) %>% unlist()
      
      return(formula_res)
    }else{return(NA)}
  })) 
  
}

table(Markers$Class)
CompareList = c("Cluster_1 vs Cluster_2")
CompareGroup <- tibble::tibble(CompareList=c("Cluster_1 vs Cluster_2"),CompareResult = list(V1=Markers))

#GO
CompareGroup <- GO_BP_Fun(CompareGroup,groups=c("Cluster_1","Cluster_2"),OrgData="org.Hs.eg.db",pcut=0.05,qcut=0.2)
GOresult <- CompareGroup$GO_Compare$V1@compareClusterResult
p <- dotplot(CompareGroup$GO_Compare$V1, showCategory = 10)+
  scale_color_gradientn(colors= colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"),space="rgb")(50))

openxlsx::write.xlsx(p$data,file=paste(OutPath,"/",folder,"/SFig2l.xlsx",sep=""))

#KEGG
CompareGroup <- KEGG_Fun(CompareGroup,groups=c("Cluster_1","Cluster_2"),OrgSp="hsa",pcut=0.05,qcut=0.2)
KEGGresult <- CompareGroup$KEGG_Compare$V1@compareClusterResult

p <- dotplot(CompareGroup$KEGG_Compare$V1, showCategory = 10)+
  scale_color_gradientn(colors= colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"),space="rgb")(50))

openxlsx::write.xlsx(p$data,file=paste(OutPath,"/",folder,"/SFig2m.xlsx",sep=""))

### n
library(monocle,lib.loc = "/home/yhdu/R/x86_64-pc-linux-gnu-library/4.1/")
HSMM_HNSCC_epi_tumor <- readRDS('/work/smy/Project/HNSCC_26sample/2.data/26sample_epi_tumor_monocleFigure.rds')

PPercentage <- table(pData(HSMM_HNSCC_epi_tumor)[,c('State2',"stage")]) %>% data.frame %>%
               #Por：用一群细胞silent/expression占比(total:100%)
               dplyr::mutate(Por= Freq/apply(table(pData(HSMM_HNSCC_epi_tumor)[,c('State2',"stage")]),1,sum))
PPercentage$State2 <- paste0("S",PPercentage$State2)
names(PPercentage)[1] <- "State"
openxlsx::write.xlsx(PPercentage,file=paste(OutPath,"/",folder,"/SFig2n.xlsx",sep=""))

### o 没做
library(velocyto.R,lib.loc = "/home/yhdu/R/x86_64-pc-linux-gnu-library/4.1/")
Tb_emb <- readRDS('/work/smy/Project/HNSCC_26sample/2.data/26Sample_epi_velocity_Tb_emb.rds')
Tb_rvel.cd <- readRDS('/work/smy/Project/HNSCC_26sample/2.data/26Sample_epi_velocity_Tb_rvel_cd.rds')
Tb_cell.colors <- readRDS('/work/smy/Project/HNSCC_26sample/2.data/26Sample_epi_velocity_Tb_cell_colors.rds')
  
show.velocity.on.embedding.cor(Tb_emb,Tb_rvel.cd,n = 2000, scale = "sqrt", cell.colors = ac(x = Tb_cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)
### p
#/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Code/40.Monocle_P2.R
p1 <- plot_cell_trajectory(HSMM_HNSCC_epi_tumor, 
                           theta = -70,
                           show_branch_points = F,
                           show_tree = TRUE, color_by = "Pseudotime", cell_size = 0.8) +
  scale_color_gradientn(colors= colorRampPalette(RColorBrewer::brewer.pal(20,"YlOrRd"),space="rgb")(50),name = 'Pseudotime')+
  theme(legend.position = "bottom")
openxlsx::write.xlsx(p1$data,file=paste(OutPath,"/",folder,"/SFig2p1.xlsx",sep=""))

df <- pData(HSMM_HNSCC_epi_tumor)
df$State2 <- ifelse(df$State %in% c("3","4","6","5","7"),"S2",ifelse(df$State=="1","S1","S3"))
pData(HSMM_HNSCC_epi_tumor) <- df
p2 <- plot_cell_trajectory(HSMM_HNSCC_epi_tumor, 
                     theta = -70,
                     show_branch_points = F,
                     show_tree = TRUE, color_by = "State2", cell_size = 0.8) +
      scale_color_manual(values = RColorBrewer::brewer.pal(3,'Set2'))+
      theme(legend.position = "bottom")
openxlsx::write.xlsx(p2$data,file=paste(OutPath,"/",folder,"/SFig2p2.xlsx",sep=""))



PPercentage <- table(pData(HSMM_HNSCC_epi_tumor)[,c('State2',"RNA_snn_res_final")]) %>% data.frame %>%
  #Por：用一群细胞silent/expression占比(total:100%)
  dplyr::mutate(Por= Freq/apply(table(pData(HSMM_HNSCC_epi_tumor)[,c('State2',"RNA_snn_res_final")]),1,sum))
openxlsx::write.xlsx(PPercentage,file=paste(OutPath,"/",folder,"/SFig2p3.xlsx",sep=""))


p4 <- plot_cell_trajectory(HSMM_HNSCC_epi_tumor, 
                     theta = -70,
                     show_branch_points = F,
                     show_tree = TRUE, color_by = "RNA_snn_res_final", cell_size = 0.8) +
  scale_color_manual(values = RColorBrewer::brewer.pal(10,'Dark2'))+
  theme(legend.position = "bottom") + facet_wrap('~RNA_snn_res_final',nrow = 1)

openxlsx::write.xlsx(p4$data,file=paste(OutPath,"/",folder,"/SFig2p4.xlsx",sep=""))

### q
# /work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Code/41.Monocle_P10Re.R
p1 <- plot_cell_trajectory(HSMM_HNSCC_epi_tumor, 
                           theta = -70,
                           show_branch_points = F,
                           show_tree = TRUE, color_by = "Pseudotime", cell_size = 0.8) +
  scale_color_gradientn(colors= colorRampPalette(RColorBrewer::brewer.pal(20,"YlOrRd"),space="rgb")(50),name = 'Pseudotime')+
  theme(legend.position = "bottom")
openxlsx::write.xlsx(p1$data,file=paste(OutPath,"/",folder,"/SFig2q1.xlsx",sep=""))

df <- pData(HSMM_HNSCC_epi_tumor)
df$State2 <- ifelse(df$State == "1","S1",ifelse(df$State=="7","S2","S3"))
pData(HSMM_HNSCC_epi_tumor) <- df
p2 <- plot_cell_trajectory(HSMM_HNSCC_epi_tumor, 
                           theta = -70,
                           show_branch_points = F,
                           show_tree = TRUE, color_by = "State2", cell_size = 0.8) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3,'Set2'))+
  theme(legend.position = "bottom")
openxlsx::write.xlsx(p2$data,file=paste(OutPath,"/",folder,"/SFig2q2.xlsx",sep=""))



PPercentage <- table(pData(HSMM_HNSCC_epi_tumor)[,c('State2',"RNA_snn_res_final")]) %>% data.frame %>%
  #Por：用一群细胞silent/expression占比(total:100%)
  dplyr::mutate(Por= Freq/apply(table(pData(HSMM_HNSCC_epi_tumor)[,c('State2',"RNA_snn_res_final")]),1,sum))
openxlsx::write.xlsx(PPercentage,file=paste(OutPath,"/",folder,"/SFig2q3.xlsx",sep=""))


p4 <- plot_cell_trajectory(HSMM_HNSCC_epi_tumor, 
                           theta = -70,
                           show_branch_points = F,
                           show_tree = TRUE, color_by = "RNA_snn_res_final", cell_size = 0.8) +
  scale_color_manual(values = RColorBrewer::brewer.pal(10,'Dark2'))+
  theme(legend.position = "bottom") + facet_wrap('~RNA_snn_res_final',nrow = 1)

openxlsx::write.xlsx(p4$data,file=paste(OutPath,"/",folder,"/SFig2q4.xlsx",sep=""))

### r
p1 <- plot_cell_trajectory(HSMM_HNSCC_epi_tumor, 
                           theta = -70,
                           show_branch_points = F,
                           show_tree = TRUE, color_by = "Pseudotime", cell_size = 0.8) +
  scale_color_gradientn(colors= colorRampPalette(RColorBrewer::brewer.pal(20,"YlOrRd"),space="rgb")(50),name = 'Pseudotime')+
  theme(legend.position = "bottom")
openxlsx::write.xlsx(p1$data,file=paste(OutPath,"/",folder,"/SFig2r1.xlsx",sep=""))

df <- pData(HSMM_HNSCC_epi_tumor)
df$State2 <- ifelse(df$State == "1","S1",ifelse(df$State=="2","S2","S3"))
pData(HSMM_HNSCC_epi_tumor) <- df
p2 <- plot_cell_trajectory(HSMM_HNSCC_epi_tumor, 
                           theta = -70,
                           show_branch_points = F,
                           show_tree = TRUE, color_by = "State2", cell_size = 0.8) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3,'Set2'))+
  theme(legend.position = "bottom")
openxlsx::write.xlsx(p2$data,file=paste(OutPath,"/",folder,"/SFig2r2.xlsx",sep=""))



PPercentage <- table(pData(HSMM_HNSCC_epi_tumor)[,c('State2',"RNA_snn_res_final")]) %>% data.frame %>%
  #Por：用一群细胞silent/expression占比(total:100%)
  dplyr::mutate(Por= Freq/apply(table(pData(HSMM_HNSCC_epi_tumor)[,c('State2',"RNA_snn_res_final")]),1,sum))
openxlsx::write.xlsx(PPercentage,file=paste(OutPath,"/",folder,"/SFig2r3.xlsx",sep=""))


p4 <- plot_cell_trajectory(HSMM_HNSCC_epi_tumor, 
                           theta = -70,
                           show_branch_points = F,
                           show_tree = TRUE, color_by = "RNA_snn_res_final", cell_size = 0.8) +
  scale_color_manual(values = RColorBrewer::brewer.pal(10,'Dark2'))+
  theme(legend.position = "bottom") + facet_wrap('~RNA_snn_res_final',nrow = 1)
