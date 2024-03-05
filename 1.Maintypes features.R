# Maintypes features at different stage
HNSCC_Whole <- readRDS('/work/smy/Project/HNSCC/2.data/DeinfeTypes/HNSCC_26sampleAll_DefineTypes.rds.gz')
HNSCC_Whole$Subcluster2 <- ifelse(HNSCC_Whole$Subcluster == 'CD4 Tcells'|HNSCC_Whole$Subcluster == 'CD8 T cells','T cells',HNSCC_Whole$Subcluster)
###Fig1=======================================================
##figure1B 
Idents(HNSCC_Whole) <- HNSCC_Whole$Subcluster2
DimPlot(HNSCC_Whole,cols = c('#DE901C','#0D5EA4','#C9C9C9','#C64C1E','#BF6196','#158F61','#6B297B','#A8AC1C','#AC142E'))

HNSCC_Whole$Subcluster2 <- factor(HNSCC_Whole$Subcluster2,levels = unique(HNSCC_Whole$Subcluster2))
p <- DimPlot(HNSCC_Whole,cols = RColorBrewer::brewer.pal(10,'Set1'),raster=FALSE)
Umap <- HNSCC_Whole@reductions$umap@cell.embeddings[,c('UMAP_1','UMAP_2')] %>% as.data.frame()
Umap$Subcluster <- HNSCC_Whole@meta.data$Subcluster2

openxlsx::write.xlsx(p$data,file=paste(OutPath,"/",folder,"/Fig1B.xlsx",sep=""))

##figure1C
p <- DotPlot(HNSCC_Whole,features = c('CDH5','PECAM1',"COL1A1","COL3A1",'GPM6B','S100B','MCAM','RGS5','CD19','MS4A1','SDC1','MZB1','CD14','FCGR3A',
                                 'CD3E','CD3G','CD3D','EPCAM','CDH1'),cols = 'Spectral') + theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1)) 

openxlsx::write.xlsx(p$data,file=paste(OutPath,"/",folder,"/Fig1C.xlsx",sep=""))

##figure1D
p <- DimPlot(HNSCC_Whole,cols = RColorBrewer::brewer.pal(10,'Set1'),split.by = "stage")

openxlsx::write.xlsx(p$data,file=paste(OutPath,"/",folder,"/Fig1D.xlsx",sep=""))

##figure1E(left panel)
PPercentage <- table(HNSCC_Whole@meta.data[,c('Subcluster2',"stage")]) %>% data.frame %>%
  #Por：用一群细胞silent/expression占比(total:100%)
  dplyr::mutate(Por= Freq/apply(table(HNSCC_Whole@meta.data[,c('Subcluster2',"stage")]),1,sum))
openxlsx::write.xlsx(PPercentage,file=paste(OutPath,"/",folder,"/Fig1Ekleft.xlsx",sep=""))

##figure1E(right panel)
number <- table(HNSCC_Whole$Subcluster2) %>% as.data.frame()
number <- number[order(number$Freq,decreasing = T),]
openxlsx::write.xlsx(number,file=paste(OutPath,"/",folder,"/Fig1Eright.xlsx",sep=""))

##figureS1B ======
subF <- HNSCC_Whole@meta.data
SubCellsDisA <- table(subF[,c("stage","Subcluster2")]) %>% 
  data.frame %>% set_colnames(c("Stage","CellTypes","Number"))

SubCellsDisA_Tissue <- lapply(split(SubCellsDisA,SubCellsDisA$Stage),function(X){
  X%>% dplyr::mutate(Per=100*Number/sum(Number))
}) %>% dplyr::bind_rows(.)

SubCellsDisA_Tissue$Stage <- as.character(SubCellsDisA_Tissue$Stage)
SubCellsDisA_Tissue$CellTypes <- factor(SubCellsDisA_Tissue$CellTypes,levels = names(table(HNSCC_Whole$Subcluster2)))

SubCellsDisA_Tissue$CellTypes <- factor(SubCellsDisA_Tissue$CellTypes,levels = unique(HNSCC_Whole$Subcluster2))

SubCellsDisA_Tissues <- arrange(SubCellsDisA_Tissue, CellTypes)

openxlsx::write.xlsx(SubCellsDisA_Tissues,file=paste(OutPath,"/",folder,"/SFig1B.xlsx",sep=""))
