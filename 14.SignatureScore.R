### CAF subtypes siganture 
library(dplyr)
WorkPath <- "/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Fig/Signature/"
### Reference 
CAFs <- openxlsx::read.xlsx("/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Data/Clinical Cancer Research, 2021_CAF subtypes.xlsx",sheet = "Table S3")
names(CAFs) <- CAFs[1,]
CAFsub <- CAFs[-1,]

Fib_color_panel <- c("Fibroblast" = '#A8373D',
                     "RSPO1+ fibroblast" = '#F1B998',
                     "POSTN+ fibroblast" = '#43739F',
                     "SFRP1+ fibroblast" = '#374D74',
                     "SEMA4A+ fibroblast" = '#947A7B',
                     "CCL19+ fibroblast"  = '#91553D',
                     "DES+ myofibroblast" = '#D48054',
                     "Proliferating fibroblast" = '#E09194')

### Fibroblast HNSC
HNSCC_fibro <- readRDS("/work/yye/Project/Collaboration/HNSC/Stroma/HNSC_Fibro1_DefineTypes.rds.gz")
Idents(HNSCC_fibro) <- HNSCC_fibro$DefineTypes

cols. = c('#A8373D','#F1B998','#43739F','#374D74','#947A7B','#91553D','#D48054','#E09194')

pdf(paste0(WorkPath,'Dimplot_fibroblast.pdf'),6,4)
DimPlot(HNSCC_fibro,reduction = 'umap',cols = cols.)
dev.off()

### 换成文章里面热图的gene list 
#Puram.CAF1 <- c("CTHRC1","COL1A1","COL3A1","POSTN","MFAP2")
#Puram.CAF2<- c("CXCL12","NDUFA4L2")
PuramCAF <- openxlsx::read.xlsx("/server1_work/brj/Collaboration/scRNA/2022/HNSCC/Data/Puram_CAF.xlsx",sheet = "Sheet1")
Puram.CAF1 <- PuramCAF$CAF1
Puram.CAF2 <- PuramCAF$CAF2
Puram.MyoFib <- c("ACTA2","MYLK","MYL9","MCAM","IL6","PDGFA")

HNSCC_fibro@meta.data$Puram.CAF1 <- apply(HNSCC_fibro@assays$RNA@data[rownames(HNSCC_fibro@assays$RNA@data) %in% Puram.CAF1,],2,mean)
HNSCC_fibro@meta.data$Puram.CAF2 <- apply(HNSCC_fibro@assays$RNA@data[rownames(HNSCC_fibro@assays$RNA@data) %in% Puram.CAF2,],2,mean)
HNSCC_fibro@meta.data$Puram.MyoFib <- apply(HNSCC_fibro@assays$RNA@data[rownames(HNSCC_fibro@assays$RNA@data) %in% Puram.MyoFib,],2,mean)


# 1.dotplot =======================================================
PuramCAFdf <- as.data.frame(HNSCC_fibro@meta.data) %>% 
              dplyr::select(c(DefineTypes,Puram.CAF1,Puram.CAF2,Puram.MyoFib)) 

PuramCAF <- PuramCAFdf %>% tidyr::gather(key = CellTypes, value = Score,-DefineTypes)


p1 <- ggplot(PuramCAF,aes(x=DefineTypes,y=CellTypes))+
      geom_point(aes(color=Score))+
      scale_color_gradientn(colors= rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"),space="rgb")(50))) +
      theme(panel.background=element_rect(color="black",fill=NA),legend.key =  element_blank(),axis.ticks=element_blank(),
      panel.grid=element_blank(),axis.text.x=element_text(size=12,angle=45,vjust = 1,hjust = 1,color="black"),
      axis.text.y=element_text(size=10,color="black"),axis.title=element_blank(),legend.position = "bottom")+
      geom_tile(data=PuramCAF,aes(x=DefineTypes,y=CellTypes),fill=NA,color="lightgray")

pdf(paste0(WorkPath,'Puram.CAFsub_DotPlot.pdf'),5,4)
p1
dev.off()

# 2.Vlnplot with p  =======================================================
sobjlists = FetchData(object = HNSCC_fibro, vars = names(PuramCAFdf))

pp <- lapply(c(2:length(names(PuramCAFdf))), function(sub){
  p1 <- ggplot(sobjlists,aes(x= DefineTypes, y = sobjlists[,sub])) + 
        geom_violin(aes(fill=DefineTypes),alpha=0.3) +   
        scale_fill_manual(values = Fib_color_panel,guide = F)+
        geom_boxplot(aes(fill=DefineTypes),width=0.2, alpha=1,outlier.colour = NA) + 
        labs(y= as.character(names(sobjlists)[sub]), x = "DefineTypes") + 
        theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5),
              axis.text.x=element_text(size=12,angle=45,vjust = 1,hjust = 1,color="black"),
              axis.text.y=element_text(size=9,color="black"))+
        ggpubr::stat_compare_means(method = "anova",data = sobjlists) # Add pairwise comparisons p-value
  return(p1)
})


pdf(paste0(OutPath,'Puram.CAFsub_fibroblast_VlnPlot.pdf'),width = 15,height = 5)
cowplot::plot_grid(plotlist = pp,ncol = 3)
dev.off()

### FeaturePlot
pdf(paste0(WorkPath,'Puram.CAFsub_fibroblast_FeaturePlot.pdf'),15,5)
FeaturePlot(HNSCC_fibro,reduction = 'umap',feature = c('Puram.CAF1','Puram.CAF2','Puram.MyoFib'),cols = c("gray","orange","red"),ncol = 3)
dev.off()


### Vln
pdf(paste0(WorkPath,'Puram.CAFsub_fibroblast_VlnPlot.pdf'),15,5)
VlnPlot(HNSCC_fibro,feature = c('Puram.CAF1','Puram.CAF2','Puram.MyoFib'),
        cols = c('#A8373D','#F1B998','#43739F','#374D74','#947A7B','#91553D','#D48054','#E09194'),pt.size = 0,ncol = 3)+
        ggpubr::compare_means(method = "anova") + labs(x="") + geom_boxplot(width=0.1) 
dev.off()


VlnPlot(HNSCC_fibro,feature = c('Puram.CAF1'),cols = c('#A8373D','#F1B998','#43739F','#374D74','#947A7B','#91553D','#D48054','#E09194'),pt.size = 0,ncol = 3)+
ggpubr::stat_compare_means(comparisons = list(c("POSTN+ fibroblast","CCL19+ fibroblast")),label = "p.format") + labs(x="")+geom_boxplot(width=0.1) 



### pan_myCAF
pan_myCAF = CAFsub$`pan-myCAF`
pan_dCAF = CAFsub$`pan-dCAF`
pan_iCAF = CAFsub$`pan-iCAF`
pan_iCAF_2 = CAFsub$`pan-iCAF-2`
pan_pCAF = CAFsub$`pan-pCAF`


HNSCC_fibro@meta.data$pan_myCAF <- apply(HNSCC_fibro@assays$RNA@data[rownames(HNSCC_fibro@assays$RNA@data) %in% pan_myCAF,],2,mean)
HNSCC_fibro@meta.data$pan_dCAF <- apply(HNSCC_fibro@assays$RNA@data[rownames(HNSCC_fibro@assays$RNA@data) %in% pan_dCAF,],2,mean)
HNSCC_fibro@meta.data$pan_iCAF <- apply(HNSCC_fibro@assays$RNA@data[rownames(HNSCC_fibro@assays$RNA@data) %in% pan_iCAF,],2,mean)
HNSCC_fibro@meta.data$pan_iCAF_2 <- apply(HNSCC_fibro@assays$RNA@data[rownames(HNSCC_fibro@assays$RNA@data) %in% pan_iCAF_2,],2,mean)
HNSCC_fibro@meta.data$pan_pCAF <- apply(HNSCC_fibro@assays$RNA@data[rownames(HNSCC_fibro@assays$RNA@data) %in% pan_pCAF,],2,mean)

# 1.dotplot =======================================================
PuramCAFdf <- as.data.frame(HNSCC_fibro@meta.data) %>% 
              dplyr::select(c(DefineTypes,pan_myCAF,pan_dCAF,pan_iCAF,pan_iCAF_2,pan_pCAF)) 
xorder <- names(PuramCAFdf)[-1]
#PuramCAFdf[,xorder] <- t(apply(PuramCAFdf[,xorder],1,scale))
#PuramCAFdf[,xorder] <- apply(PuramCAFdf[,xorder],2,scale)

PuramCAF <- PuramCAFdf %>% tidyr::gather(key = CellTypes, value = Score,-DefineTypes)
p1 <- ggplot(PuramCAF,aes(x=DefineTypes,y=CellTypes))+
  geom_point(aes(color=Score))+
  scale_color_gradientn(colors= rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"),space="rgb")(50))) +
  theme(panel.background=element_rect(color="black",fill=NA),legend.key =  element_blank(),axis.ticks=element_blank(),
        panel.grid=element_blank(),axis.text.x=element_text(size=12,angle=45,vjust = 1,hjust = 1,color="black"),
        axis.text.y=element_text(size=10,color="black"),axis.title=element_blank(),legend.position = "bottom")+
  geom_tile(data=PuramCAF,aes(x=DefineTypes,y=CellTypes),fill=NA,color="lightgray")

pdf(paste0(WorkPath,'Five.CAFsub_DotPlot.pdf'),5,5)
p1
dev.off()

# 2.Vlnplot with p  =======================================================
sobjlists = FetchData(object = HNSCC_fibro, vars = names(PuramCAFdf))

pp <- lapply(c(2:length(names(PuramCAFdf))), function(sub){
  p1 <- ggplot(sobjlists,aes(x= DefineTypes, y = sobjlists[,sub])) + 
    geom_violin(aes(fill=DefineTypes),alpha=0.3) +   
    scale_fill_manual(values = Fib_color_panel,guide = F)+    
    geom_boxplot(aes(fill=DefineTypes),width=0.1,alpha=1,alpha=1,outlier.colour = NA) + # alpha=0 可以让中间的boxplot的颜色变为跟vlnplot一样的
    labs(y= as.character(names(sobjlists)[sub]), x = "DefineTypes") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5),
          axis.text.x=element_text(size=12,angle=45,vjust = 1,hjust = 1,color="black"),
          axis.text.y=element_text(size=9,color="black"))+
    ggpubr::stat_compare_means(method = "anova",data = sobjlists) # Add pairwise comparisons p-value
  return(p1)
})


pdf(paste0(OutPath,'Five.CAFsub_fibroblast_VlnPlot.pdf'),width = 25,height = 5)
cowplot::plot_grid(plotlist = pp,ncol = 5)
dev.off()

### FeaturePlot
pdf(paste0(WorkPath,'5CAFsub_fibroblast_FeaturePlot.pdf'),25,5)
FeaturePlot(HNSCC_fibro,reduction = 'umap',feature = c('pan_myCAF','pan_dCAF','pan_iCAF',"pan_iCAF_2","pan_pCAF"),cols = c("gray","orange","red"),ncol = 5)
dev.off()


### Vln
pdf(paste0(WorkPath,'5CAFsub_fibroblast_VlnPlot.pdf'),25,5)
VlnPlot(HNSCC_fibro,feature = c('pan_myCAF','pan_dCAF','pan_iCAF',"pan_iCAF_2","pan_pCAF"),cols = c('#A8373D','#F1B998','#43739F','#374D74','#947A7B','#91553D','#D48054','#E09194'),pt.size = 0,ncol = 5)
dev.off()


### iCAF, myCAF, ApCAF, M1 or M2 macrophages  
CAF3s <- openxlsx::read.xlsx("/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Data/mouseCAF.xlsx",sheet = "Sheet1")
Homo <- read.csv("/server1_work/yye/Project/Data/Public/Reference/mouse.human.homology.genes.csv")[,c("mgi_symbol","hgnc_symbol")] %>% unique
iCAFs <- CAF3s[,1] %>% as.data.frame() %>% set_colnames("gene")
iCAFs <- merge(iCAFs,Homo,by.x="gene",by.y="mgi_symbol")

myCAFs <- CAF3s[,3] %>% as.data.frame() %>% set_colnames("gene")
myCAFs <- merge(myCAFs,Homo,by.x="gene",by.y="mgi_symbol")

ApCAFs <- CAF3s[,2] %>% as.data.frame() %>% set_colnames("gene")
ApCAFs <- merge(ApCAFs,Homo,by.x="gene",by.y="mgi_symbol")

iCAF <- iCAFs$hgnc_symbol
myCAF <- myCAFs$hgnc_symbol
ApCAF <- ApCAFs$hgnc_symbol


HNSCC_fibro@meta.data$iCAF <- apply(HNSCC_fibro@assays$RNA@data[rownames(HNSCC_fibro@assays$RNA@data) %in% iCAF,],2,mean)
HNSCC_fibro@meta.data$myCAF <- apply(HNSCC_fibro@assays$RNA@data[rownames(HNSCC_fibro@assays$RNA@data) %in% myCAF,],2,mean)
HNSCC_fibro@meta.data$ApCAF <- apply(HNSCC_fibro@assays$RNA@data[rownames(HNSCC_fibro@assays$RNA@data) %in% ApCAF,],2,mean)



### FeaturePlot
pdf(paste0(WorkPath,'3CAFsub_fibroblast_FeaturePlot_FromMouse.pdf'),16.5,5)
FeaturePlot(HNSCC_fibro,reduction = 'umap',feature = c('iCAF','myCAF','ApCAF'),order = T,cols = c("gray","orange","red"),ncol = 3)
dev.off()

### Vln
pdf(paste0(WorkPath,'3CAFsub_fibroblast_VlnPlot_FromMouse.pdf'),15,6)
VlnPlot(HNSCC_fibro,feature = c('iCAF','myCAF','ApCAF'),cols = c('#A8373D','#F1B998','#43739F','#374D74','#947A7B','#91553D','#D48054','#E09194'),pt.size = 0,ncol = 3)
dev.off()

##  CANCER DISGOVERY marker
MHCII <- grep("HLA-DR|HLA-DQ|HLA-DP",rownames(HNSCC_fibro@assays$RNA@data),value = T)

iCAFh <- c("IL6","PDGFRA","CXCL12","CFD","DPT","LMNA","AGTR","HAS1","CXCL1","CXCL2","CCL2","IL8") 
myCAFh <- c("ACTA2","TAGLN","MMP11","MYL9","HOPX","POSTN","TPM1","TPM2")
apCAFh <- c("CD74","HLA-DRA","HLA-DPA1","HLA-DQA1","SLPI")


HNSCC_fibro@meta.data$iCAF_Human <- apply(HNSCC_fibro@assays$RNA@data[rownames(HNSCC_fibro@assays$RNA@data) %in% iCAFh,],2,mean)
HNSCC_fibro@meta.data$myCAF_Human <- apply(HNSCC_fibro@assays$RNA@data[rownames(HNSCC_fibro@assays$RNA@data) %in% myCAFh,],2,mean)
HNSCC_fibro@meta.data$ApCAF_Human <- apply(HNSCC_fibro@assays$RNA@data[rownames(HNSCC_fibro@assays$RNA@data) %in% apCAFh,],2,mean)

# 1.dotplot =======================================================
PuramCAFdf <- as.data.frame(HNSCC_fibro@meta.data) %>% 
              dplyr::select(c(DefineTypes,iCAF_Human,myCAF_Human,ApCAF_Human)) 
#xorder <- names(PuramCAFdf)[-1]
#PuramCAFdf[,xorder] <- t(apply(PuramCAFdf[,xorder],1,scale))
PuramCAF <- PuramCAFdf %>% tidyr::gather(key = CellTypes, value = Score,-DefineTypes)
p1 <- ggplot(PuramCAF,aes(x=DefineTypes,y=CellTypes))+
  geom_point(aes(color=Score))+
  scale_color_gradientn(colors= rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"),space="rgb")(50))) +
  theme(panel.background=element_rect(color="black",fill=NA),legend.key =  element_blank(),axis.ticks=element_blank(),
        panel.grid=element_blank(),axis.text.x=element_text(size=12,angle=45,vjust = 1,hjust = 1,color="black"),
        axis.text.y=element_text(size=10,color="black"),axis.title=element_blank(),legend.position = "bottom")+
  geom_tile(data=PuramCAF,aes(x=DefineTypes,y=CellTypes),fill=NA,color="lightgray")

pdf(paste0(WorkPath,'ThreeHuman.CAFsub_DotPlot.pdf'),5,4)
p1
dev.off()


# 2.Vlnplot with p  =======================================================
sobjlists = FetchData(object = HNSCC_fibro, vars = names(PuramCAFdf))
pp <- lapply(c(2:length(names(PuramCAFdf))), function(sub){
  p1 <- ggplot(sobjlists,aes(x= DefineTypes, y = sobjlists[,sub])) + 
    geom_violin(aes(fill=DefineTypes),alpha=0.3) +   
    scale_fill_manual(values = Fib_color_panel,guide = F)+        
    geom_boxplot(aes(fill=DefineTypes),width=0.1,alpha=1,outlier.colour = NA) + 
    labs(y= as.character(names(sobjlists)[sub]), x = "DefineTypes") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5),
          axis.text.x=element_text(size=12,angle=45,vjust = 1,hjust = 1,color="black"),
          axis.text.y=element_text(size=9,color="black"))+
    ggpubr::stat_compare_means(method = "anova",data = sobjlists) # Add pairwise comparisons p-value
  return(p1)
})
pdf(paste0(OutPath,'Three.CAFsub_fibroblast_VlnPlot.pdf'),width =15,height = 5)
cowplot::plot_grid(plotlist = pp,ncol = 3)
dev.off()



### FeaturePlot
pdf(paste0(WorkPath,'3CAFsub_fibroblast_FeaturePlot.pdf'),15,5)
FeaturePlot(HNSCC_fibro,reduction = 'umap',feature = c('iCAF_Human','myCAF_Human','ApCAF_Human'),order = T,cols = c("gray","orange","red"),ncol = 3)
dev.off()

### Vln
pdf(paste0(WorkPath,'3CAFsub_fibroblast_VlnPlot.pdf'),15,5)
VlnPlot(HNSCC_fibro,feature = c('iCAF_Human','myCAF_Human','ApCAF_Human'),cols = c('#A8373D','#F1B998','#43739F','#374D74','#947A7B','#91553D','#D48054','#E09194'),pt.size = 0,ncol = 3)
dev.off()

### Macrophages 
HNSCC_myeloids <- readRDS('/work/smy/Project/HNSCC/2.data/DeinfeTypes/26Sample_Myeloid_DefineTypesFinal.rds')
HNSCC_myeloid <- subset(HNSCC_myeloids,subset = DefineTypes %in% c("CXCL10+ macrophages","C1QC+MRC-macrophages","SPP1+ macorphages","FOLR2+ macorphages") )
HNSCC_myeloid$DefineTypes <- factor(HNSCC_myeloid$DefineTypes,levels = c("CXCL10+ macrophages","C1QC+MRC-macrophages","SPP1+ macorphages","FOLR2+ macorphages"))
Idents(HNSCC_myeloid) <- HNSCC_myeloid$DefineTypes

pdf(paste0(WorkPath,'Dimplot_Myeloid.pdf'),7,4)
DimPlot(HNSCC_myeloid,reduction = 'umap',cols = c('#79545C','#FAD181','#483D4D','#A44E41'))
dev.off()

M1 <- c("IL23","TNF","CXCL9","CXCL10","CXCL11","CD86","IL1A","IL1B","IL6","CCL5","IRF5","IRF1","CD40","IDO1","KYNU","CCR7")
M2 <- c("IL4R","CCL4","CCL13","CCL20","CCL17","CCL18","CCL22","CCL24","LYVE1","VEGFA","VEGFB","VEGFC","VEGFD","EGF","CTSA","CTSB","CTSC","CTSD","TGFB1","TGFB2","TGFB3","MMP14","MMP19","MMP9","CLEC7A","WNT7B","FASL","TNFSF12","TNFSF8","CD276","VTCN1","MSR1","FN1","IRF4")

HNSCC_myeloid@meta.data$M1 <- apply(HNSCC_myeloid@assays$RNA@data[rownames(HNSCC_myeloid@assays$RNA@data) %in% M1,],2,mean)
HNSCC_myeloid@meta.data$M2 <- apply(HNSCC_myeloid@assays$RNA@data[rownames(HNSCC_myeloid@assays$RNA@data) %in% M2,],2,mean)
HNSCC_myeloid@meta.data$M1_divide_M2 <- as.numeric(HNSCC_myeloid@meta.data$M1/HNSCC_myeloid@meta.data$M2)



marcotypes <- grep("macro|macor",HNSCC_myeloid@meta.data$DefineTypes,value = T) 
HNSCC_marco <- subset(HNSCC_myeloid, subset = DefineTypes %in% marcotypes )

MacroDF <- HNSCC_marco@meta.data %>% as.data.frame()

### save
readr::write_rds(MacroDF,"/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Data/Signature/MacScore.rds")


library(harmony)
HNSCC_marco <- HNSCC_marco %>% 
               RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(HNSCC_marco, 'harmony')
HNSCC_marco <- HNSCC_marco %>% 
               RunUMAP(reduction = "harmony", dims = 1:30) %>% 
               FindNeighbors(reduction = "harmony", dims = 1:30)

pdf(paste0(WorkPath,'Dimplot_Macro.pdf'),7,4)
DimPlot(HNSCC_marco,reduction = 'umap',cols = c('#79545C','#FAD181','#483D4D','#A44E41'))
dev.off()

FeaturePlot(HNSCC_marco,reduction = 'umap',feature = c( 'CXCL10','CXCL9','C1QC','MRC1','SPP1','CD68','FOLR2','CCL18'),cols = c("gray","orange","red"),ncol = 2)



### FeaturePlot
pdf(paste0(WorkPath,'M1M2_Macro_FeaturePlot.pdf'),11,5)
FeaturePlot(HNSCC_marco,reduction = 'umap',feature = c('M1','M2'),cols = c("gray","orange","red"),ncol = 2)
dev.off()


### Vln
pdf(paste0(WorkPath,'M1M2_Macro_VlnPlot.pdf'),10,6)
VlnPlot(HNSCC_marco,feature = c('M1','M2'),cols =c('#79545C','#FAD181','#483D4D','#A44E41'),pt.size = 0,ncol =2)
dev.off()

###M1,M2的Density
data = HNSCC_marco@meta.data[,c("DefineTypes","M1","M2")]
PP <- lapply(unique(marcotypes), function(marco){
  dataFilter2 <- data %>% dplyr::filter(M1 != 0 & M2 != 0 & DefineTypes== as.character(marco))
  M1M2out <-ggplot(dataFilter2,aes(x=M1,y=M2)) +
    geom_point(alpha=0)+stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +xlim(0,1)+ylim(0,3)+
    scale_fill_distiller(palette="YlOrRd", direction=1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank(),# x,y刻度去掉
          axis.line.x = element_blank(),axis.line.y = element_blank()) # x, y 轴线去掉
})
pdf(paste0(WorkPath,'M1M2_Density.pdf'),20,4)
cowplot::plot_grid(plotlist = PP,labels = unique(marcotypes),rel_widths = c(5,5),ncol = 4)
dev.off()

###M1,M2的heatmap
data = HNSCC_marco@meta.data[,c("DefineTypes","M1","M2","M1_divide_M2")]
data <- data %>% dplyr::filter(!is.infinite(M1_divide_M2)&!is.na(M1_divide_M2))
DefineTypes = unique(data$DefineTypes)
test = t(as.data.frame(c("SPP1+ HNSCC_marcophages",colMeans(data[data$DefineTypes=="SPP1+ HNSCC_marcophages",2:4]))))

for (i in DefineTypes) {
  test1 = t(as.data.frame(c(i,colMeans(data[data$DefineTypes==i,2:4]))))
  test = rbind(test,test1) 
}
colMeans(data[,2:4])
test=test[-1,]
test <- as.data.frame(test)
test[,c(2:4)] <- apply(test[,c(2:4)],2,as.numeric) 
colnames(test)=c("DefineTypes","M1","M2","M1_divide_M2")
rownames(test) <- test$DefineTypes
plotdf <- test[,-1]
plotdfe <- plotdf[,-3]  
# 按mac order and out first 
library(pheatmap)
p1 <- pheatmap(plotdfe,scale = "column",
               cluster_rows = F,
               show_rownames = T,
               cluster_cols = F,
               color = colorRampPalette(RColorBrewer::brewer.pal(9, "OrRd"))(50))

pdf(paste0(WorkPath,"M1_M2Heatmap.pdf"),width = 4,height = 3)
print(p1)
dev.off()

# 1.dotplot =======================================================
Macrophagedf <- as.data.frame(HNSCC_myeloid@meta.data) %>% 
  dplyr::select(c(DefineTypes,M1,M2)) 

#xorder <- names(Macrophagedf)[-1]
#Macrophagedf[,xorder] <- t(apply(Macrophagedf[,xorder],1,scale))
Macrophage <- Macrophagedf %>% tidyr::gather(key = CellTypes, value = Score,-DefineTypes)
p1 <- ggplot(Macrophage,aes(x=DefineTypes,y=CellTypes))+
  geom_point(aes(color=Score))+
  scale_color_gradientn(colors= rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"),space="rgb")(50)))+
  theme(panel.background=element_rect(color="black",fill=NA),legend.key =  element_blank(),axis.ticks=element_blank(),
        panel.grid=element_blank(),axis.text.x=element_text(size=12,angle=45,vjust = 1,hjust = 1,color="black"),
        axis.text.y=element_text(size=10,color="black"),axis.title=element_blank(),legend.position = "bottom")+
  geom_tile(data=Macrophage,aes(x=DefineTypes,y=CellTypes),fill=NA,color="lightgray")

pdf(paste0(WorkPath,'M1M2_DotPlot.pdf'),5,4)
p1
dev.off()


# 2.Vlnplot with p  =======================================================
sobjlists = FetchData(object = HNSCC_myeloid, vars = names(Macrophagedf))

pp <- lapply(c(2:length(names(Macrophagedf))), function(sub){
  p1 <- ggplot(sobjlists,aes(x= DefineTypes, y = sobjlists[,sub])) + 
    geom_violin(aes(fill=DefineTypes), alpha=0.3 ) +   
    scale_fill_manual(limits = unique(sobjlists$DefineTypes),
                      values = c('#79545C','#FAD181','#483D4D','#A44E41'),guide = F)+
    geom_boxplot(aes(fill=DefineTypes),width=0.1, alpha=1,outlier.colour = NA) + # 不要离群值
    labs(y= as.character(names(sobjlists)[sub]), x = "DefineTypes") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5),
          axis.text.x=element_text(size=12,angle=45,vjust = 1,hjust = 1,color="black"),
          axis.text.y=element_text(size=9,color="black"))+
    ggpubr::stat_compare_means(method = "anova",data = sobjlists) # Add pairwise comparisons p-value
  return(p1)
})

pdf(paste0(WorkPath,'M1M2_VlnPlot.pdf'),width = 5,height = 5)
cowplot::plot_grid(plotlist = pp,ncol = 2)
dev.off()



### epi
HNSCC_epi_tumor <- readRDS('/work/smy/Project/HNSCC_26sample/2.data/26sample_tumorcells_Harmony_immuneReference_All.rds')
Idents(HNSCC_epi_tumor) <- HNSCC_epi_tumor$RNA_snn_res_final


TumorSignature <- openxlsx::read.xlsx("/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Data/Puram.Tumor.xlsx")

Puram.Cell.Cycle <- TumorSignature$Cell.Cycle
Puram.p.EMT <- TumorSignature$`p-EMT`
Puram.Epi.Dif1 <- TumorSignature$Epi.dif..1
Puram.Epi.Dif2 <- TumorSignature$Epi.dif..2
Puram.Stress <- TumorSignature$Stress
Puram.Hypoxia <- TumorSignature$Hypoxia


HNSCC_epi_tumor@meta.data$Puram.Cell.Cycle <- apply(HNSCC_epi_tumor@assays$RNA@data[rownames(HNSCC_epi_tumor@assays$RNA@data) %in% Puram.Cell.Cycle,],2,mean)
HNSCC_epi_tumor@meta.data$Puram.p.EMT <- apply(HNSCC_epi_tumor@assays$RNA@data[rownames(HNSCC_epi_tumor@assays$RNA@data) %in% Puram.p.EMT,],2,mean)
HNSCC_epi_tumor@meta.data$Puram.Epi.Dif1 <- apply(HNSCC_epi_tumor@assays$RNA@data[rownames(HNSCC_epi_tumor@assays$RNA@data) %in% Puram.Epi.Dif1,],2,mean)
HNSCC_epi_tumor@meta.data$Puram.Epi.Dif2 <- apply(HNSCC_epi_tumor@assays$RNA@data[rownames(HNSCC_epi_tumor@assays$RNA@data) %in% Puram.Epi.Dif2,],2,mean)
HNSCC_epi_tumor@meta.data$Puram.Stress <- apply(HNSCC_epi_tumor@assays$RNA@data[rownames(HNSCC_epi_tumor@assays$RNA@data) %in% Puram.Stress,],2,mean)
HNSCC_epi_tumor@meta.data$Puram.Hypoxia <- apply(HNSCC_epi_tumor@assays$RNA@data[rownames(HNSCC_epi_tumor@assays$RNA@data) %in% Puram.Hypoxia,],2,mean)

TumorDF <- HNSCC_epi_tumor@meta.data %>% as.data.frame()
### save
readr::write_rds(TumorDF,"/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Data/Signature/TumorScore.rds")


pdf(paste0(WorkPath,'Tumor_DimPlot.pdf'),5,5)
DimPlot(HNSCC_epi_tumor,cols = RColorBrewer::brewer.pal(10,'Dark2'))
dev.off()

### FeaturePlot
pdf(paste0(WorkPath,'Puram_tumor_FeaturePlot.pdf'),30,5)
FeaturePlot(HNSCC_epi_tumor,reduction = 'umap',feature = c('Puram.Cell.Cycle','Puram.p.EMT',"Puram.Epi.Dif1","Puram.Epi.Dif2","Puram.Stress","Puram.Hypoxia"),cols = c("gray","orange","red"),ncol = 6)
dev.off()


### Vln
pdf(paste0(WorkPath,'Puram_tumor_VlnPlot.pdf'),30,5)
VlnPlot(HNSCC_epi_tumor,feature = c('Puram.Cell.Cycle','Puram.p.EMT',"Puram.Epi.Dif1","Puram.Epi.Dif2","Puram.Stress","Puram.Hypoxia"),cols = RColorBrewer::brewer.pal(10,'Dark2'),pt.size = 0,ncol = 6)
dev.off()

# 1.dotplot =======================================================
TumorSubdf <- as.data.frame(HNSCC_epi_tumor@meta.data) %>% 
              dplyr::select(c(RNA_snn_res_final,Puram.Cell.Cycle,Puram.p.EMT,Puram.Epi.Dif1,Puram.Epi.Dif2,Puram.Stress,Puram.Hypoxia))
TumorSubdf$RNA_snn_res_final <- paste0("C",TumorSubdf$RNA_snn_res_final)
xorder <- names(TumorSubdf)[-1]
TumorSubdf[,xorder] <- t(apply(TumorSubdf[,xorder],1,scale))
#TumorSubdf[,xorder] <- apply(TumorSubdf[,xorder],2,scale)
TumorSub <- TumorSubdf %>% tidyr::gather(key = CellTypes, value = Score,-RNA_snn_res_final)

p1 <- ggplot(TumorSub,aes(x=RNA_snn_res_final,y=CellTypes))+
  geom_point(aes(color=Score))+
  scale_color_gradientn(colors= rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"),space="rgb")(50)))+
  theme(panel.background=element_rect(color="black",fill=NA),legend.key =  element_blank(),axis.ticks=element_blank(),
        panel.grid=element_blank(),axis.text.x=element_text(size=12,angle=45,vjust = 1,hjust = 1,color="black"),
        axis.text.y=element_text(size=10,color="black"),axis.title=element_blank(),legend.position = "bottom")+
  geom_tile(data=TumorSub,aes(x=RNA_snn_res_final,y=CellTypes),fill=NA,color="lightgray")

pdf(paste0(WorkPath,'Puram.Tumor_DotPlot.pdf'),5,5)
p1
dev.off()

# 2.Vlnplot with p  =======================================================
sobjlists = FetchData(object = HNSCC_epi_tumor, vars = names(TumorSubdf))
sobjlists$RNA_snn_res_final <- paste0("C",sobjlists$RNA_snn_res_final)
sobjlists$RNA_snn_res_final <- factor(sobjlists$RNA_snn_res_final,levels = c("C0","C1","C2","C3","C4"))

pp <- lapply(c(2:length(names(TumorSubdf))), function(sub){
  p1 <- ggplot(sobjlists,aes(x= RNA_snn_res_final, y = sobjlists[,sub])) + 
    geom_violin(aes(fill=RNA_snn_res_final),alpha=0.3) +   
    scale_fill_manual(limits = c("C0","C1","C2","C3","C4"),
                      values = RColorBrewer::brewer.pal(10,'Dark2'),guide = F)+
    geom_boxplot(aes(fill=RNA_snn_res_final),width=0.1, alpha=1,outlier.colour = NA) + 
    labs(y= as.character(names(sobjlists)[sub]), x = "RNA_snn_res_final") + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5),
          axis.text.x=element_text(size=12,angle=45,vjust = 1,hjust = 1,color="black"),
          axis.text.y=element_text(size=9,color="black"))+
    ggpubr::stat_compare_means(method = "anova",data = sobjlists) # Add pairwise comparisons p-value
  return(p1)
})
pdf(paste0(OutPath,'Puram.TumorVlnPlot.pdf'),width = 30,height = 5)
cowplot::plot_grid(plotlist = pp,ncol = 6)
dev.off()

### Heatmap
TumorDF <- readr::read_rds("/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Data/Signature/TumorScore.rds")
CandiCol <- c("Puram.Cell.Cycle","Puram.p.EMT","Puram.Epi.Dif1","Puram.Epi.Dif2","Puram.Stress","Puram.Hypoxia")

data <- TumorDF %>% dplyr::mutate(DefineTypes = paste0("C",RNA_snn_res_final)) %>%
  dplyr::select(c(DefineTypes,CandiCol))

test <- t(as.data.frame(c("C0",colMeans(data[data$DefineTypes=="C0",2:ncol(data)]))))
DefineTypes <- unique(data$DefineTypes)
for (i in DefineTypes) {
  test1 <- t(as.data.frame(c(i,colMeans(data[data$DefineTypes==i,2:ncol(data)]))))
  test <- rbind(test,test1) 
}
test <- test[-1,]
test <- as.data.frame(test)

test[,2:ncol(data)] <- apply(test[,c(2:ncol(data))],2,as.numeric) 

colnames(test)=c("DefineTypes",CandiCol)
rownames(test) <- test$DefineTypes
plotdf <- test[c("C0","C1","C2","C3","C4"),-1]
plotdff <- t(plotdf)

# 按mac order and out first 
library(pheatmap)
p1 <- pheatmap(plotdff,scale = "row",
               cluster_rows = F,
               show_rownames = T,
               cluster_cols = F,
               color = colorRampPalette(RColorBrewer::brewer.pal(9, "OrRd"))(50))

pdf(paste0(OutPath,"Puram_Tumor_Heatmap.pdf"),width = 4,height = 2.5)
print(p1)
dev.off()


