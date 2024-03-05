.cluster_cols <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

OutPath <- "/work/brj/Collaboration/2022/scRNA/HNSCC/Revise3/"
folder <- "Fig3"

##fibro cell=======
HNSCC_fibro <- readRDS("/work/yye/Project/Collaboration/HNSC/Stroma/HNSC_Fibro1_DefineTypes.rds.gz")
Idents(HNSCC_fibro) <- HNSCC_fibro$DefineTypes

##A
p <- DimPlot(HNSCC_fibro,reduction = 'umap',cols = c('#A8373D','#F1B998','#43739F','#374D74','#947A7B','#91553D','#D48054','#E09194'))
openxlsx::write.xlsx(p$data,file=paste(OutPath,"/",folder,"/Fig3A.xlsx",sep=""))

##B
p <-DotPlot(HNSCC_fibro,feature = c('COL1A1','COL1A2','RSPO1','CRABP1','POSTN','LAMP5','SFRP1','PLA2G2A','SEMA4A','SOD2','CCL19','CCL2','DES','MYF5','MKI67','TOP2A'),cols = 'Spectral')+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1)) 
openxlsx::write.xlsx(p$data,file=paste(OutPath,"/",folder,"/Fig3B.xlsx",sep=""))


##C
HNSCC_fibro$stage <- ifelse(HNSCC_fibro$orig.ident == 'A-OSCC1'|HNSCC_fibro$orig.ident == 'A-OSCC3'|HNSCC_fibro$orig.ident == 'A-OSCC4'|
                              HNSCC_fibro$orig.ident == 'A-OSCC5'|HNSCC_fibro$orig.ident == 'A-OSCC6'|HNSCC_fibro$orig.ident == 'A-OSCC7'|HNSCC_fibro$orig.ident == 'E-OSCC4','A',
                            ifelse(HNSCC_fibro$orig.ident == 'E-OSCC1'|HNSCC_fibro$orig.ident == 'E-OSCC2'|HNSCC_fibro$orig.ident == 'E-OSCC3','E',
                                   ifelse(HNSCC_fibro$orig.ident == 'LN1'|HNSCC_fibro$orig.ident == 'LN6','LN-out',
                                          ifelse(HNSCC_fibro$orig.ident == 'LN2'|HNSCC_fibro$orig.ident == 'LN3'|HNSCC_fibro$orig.ident == 'LN8','LN-in',
                                                 ifelse(HNSCC_fibro$orig.ident == 'LN4'|HNSCC_fibro$orig.ident == 'LN5','LN-normal',
                                                        ifelse(HNSCC_fibro$orig.ident == 'NT1'|HNSCC_fibro$orig.ident == 'NT2'|HNSCC_fibro$orig.ident == 'NT3','NT',
                                                               ifelse(HNSCC_fibro$orig.ident == 'Pre-Ca1'|HNSCC_fibro$orig.ident == 'Pre-Ca2'|HNSCC_fibro$orig.ident == 'Pre-Ca3','Pre',
                                                                      ifelse(HNSCC_fibro$orig.ident == 'R-OSCC1'|HNSCC_fibro$orig.ident == 'R-OSCC2'|
                                                                               HNSCC_fibro$orig.ident == 'R-OSCC3'|HNSCC_fibro$orig.ident == 'R-OSCC4','R','NO'))))))))
HNSCC_fibro$stage <- factor(HNSCC_fibro$stage,levels = c('NT','Pre','E','A','LN-in','LN-out','LN-normal','R'))
HNSCC_fibroSub <- subset(HNSCC_fibro,subset = stage %in% c('NT','Pre','E','A','R'))

##比例变化
subF <- HNSCC_fibroSub@meta.data
SufibrosDisA <- table(subF[,c("stage","DefineTypes")]) %>% 
  data.frame %>% set_colnames(c("Stage","CellTypes","Number"))

SufibrosDisA_Tissue <- lapply(split(SufibrosDisA,SufibrosDisA$Stage),function(X){
  X%>% dplyr::mutate(Per=100*Number/sum(Number))
}) %>% dplyr::bind_rows(.)

SufibrosDisA_Tissue$Stage <- factor(SufibrosDisA_Tissue$Stage,levels = c('NT','Pre','E','A','R'))
SufibrosDisA_Tissue$CellTypes <- factor(SufibrosDisA_Tissue$CellTypes,levels = names(table(HNSCC_fibro$DefineTypes)))

SufibrosDisA_TissueP <- arrange(SufibrosDisA_Tissue,CellTypes) %>% dplyr::filter(Number != 0 )
openxlsx::write.xlsx(SufibrosDisA_TissueP,file=paste(OutPath,"/",folder,"/Fig3C.xlsx",sep=""))

###D 
### PSOTN fib, RSPO 
GSVAscore <- readr::read_rds("/server1_work/brj/Collaboration/scRNA/2022/HNSCC/Revise/Result/HNSC_POSTN_SPP1.rds.gz")
GSVAF <- GSVAscore$POSTN_SPP1_score[[1]] %>% t() %>% as.data.frame() 
GSVAF$Mixture <- substr(rownames(GSVAF),0,12)
HNSC_SignatureF <- GSVAF[substr(rownames(GSVAF),14,16) %in% c("01A"),]
HNSC_SignatureF$Mixture <- gsub("[.]","-",HNSC_SignatureF$Mixture)

HNSC_POSTN_SPP1 <- HNSC_SignatureF %>% dplyr::select(c(Mixture,POSTNTop50,SPP1Top50))


GSVAscore1 <- readr::read_rds("/server1_work/brj/Collaboration/scRNA/2022/HNSCC/Revise/Result/HNSC_RSPO1_FOLR2.rds.gz")
GSVAF1 <- GSVAscore1$RSPO1_FOLR2_score[[1]] %>% t() %>% as.data.frame() 
GSVAF1$Mixture <- substr(rownames(GSVAF1),0,12)
HNSC_SignatureF1 <- GSVAF1[substr(rownames(GSVAF1),14,16) %in% c("01A"),]
HNSC_SignatureF1$Mixture <- gsub("[.]","-",HNSC_SignatureF1$Mixture)

HNSC_RSPO1_FOLR2 <- HNSC_SignatureF1

HNSCC_proportion <- HNSC_SignatureF %>% inner_join(HNSC_SignatureF1) %>%
  dplyr::select(c(Mixture,RSPO1Top50,POSTNTop50,FOLR2Top50,SPP1Top50))
names(HNSCC_proportion) <- gsub("Top50","",names(HNSCC_proportion))
##clinical info
surviva <- read.delim("/server1_work/brj/GEO/CRC/multiCox/TCGA_ClinicalData_20180420.txt")#,","gender"
survival <- surviva %>% dplyr::select("bcr_patient_barcode",'OS','OS.time') %>% 
  rename("Mixture" = "bcr_patient_barcode") %>%
  dplyr::filter(Mixture %in% HNSCC_proportion$Mixture) %>% 
  dplyr::mutate(OS.time = as.numeric(as.numeric(as.character(.$OS.time))/30)) %>% 
  dplyr::filter(!is.na(OS.time)) %>% inner_join(HNSCC_proportion) %>%
  rename("Sample" = "Mixture")

openxlsx::write.xlsx(survival,file=paste(OutPath,"/",folder,"/FigDH.xlsx",sep=""))



## I
CompareGroup <- readRDS('/work/smy/Project/HNSCC_26sample/2.data/POSTN vs RSPO1_GO_KEGGenrichmentResult.rds')
p <- dotplot(CompareGroup$GO_Compare$V1, showCategory = 10)+
  scale_color_gradientn(colors= colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"),space="rgb")(50))

Plot <- p$data
openxlsx::write.xlsx(Plot,file=paste(OutPath,"/",folder,"/Fig3I.xlsx",sep=""))


## E
##mye cell=======
HNSCC_myeloid <- readRDS('/work/smy/Project/HNSCC/2.data/DeinfeTypes/26Sample_Myeloid_DefineTypesFinal.rds')
HNSCC_myeloid$DefineTypes <- factor(HNSCC_myeloid$DefineTypes,levels = c("cDC1","cDC2","pDC","LAMP3+ DCs","Monocytes","Mast cells",
                                                                         "CXCL10+ macrophages","C1QC+MRC-macrophages","SPP1+ macorphages","FOLR2+ macorphages","Proliferating myeloid cells"))
Idents(HNSCC_myeloid) <- HNSCC_myeloid$DefineTypes
p <- DimPlot(HNSCC_myeloid,reduction = 'umap',cols = c('#FEE0C4','#EE948B','#F4B3B1','#5C6479','#334873','#BBA399',
                                                  '#79545C','#FAD181','#483D4D','#A44E41','#F1AC73'))
openxlsx::write.xlsx(p$data,file=paste(OutPath,"/",folder,"/Fig3E.xlsx",sep=""))
## G
HNSCC_myeloid$stage <- ifelse(HNSCC_myeloid$orig.ident == 'A-OSCC1'|HNSCC_myeloid$orig.ident == 'A-OSCC3'|HNSCC_myeloid$orig.ident == 'A-OSCC4'|
                                HNSCC_myeloid$orig.ident == 'A-OSCC5'|HNSCC_myeloid$orig.ident == 'A-OSCC6'|HNSCC_myeloid$orig.ident == 'A-OSCC7'|HNSCC_myeloid$orig.ident == 'E-OSCC4','A',
                              ifelse(HNSCC_myeloid$orig.ident == 'E-OSCC1'|HNSCC_myeloid$orig.ident == 'E-OSCC2'|HNSCC_myeloid$orig.ident == 'E-OSCC3','E',
                                     ifelse(HNSCC_myeloid$orig.ident == 'LN1'|HNSCC_myeloid$orig.ident == 'LN6','LN-out',
                                            ifelse(HNSCC_myeloid$orig.ident == 'LN2'|HNSCC_myeloid$orig.ident == 'LN3'|HNSCC_myeloid$orig.ident == 'LN8','LN-in',
                                                   ifelse(HNSCC_myeloid$orig.ident == 'LN4'|HNSCC_myeloid$orig.ident == 'LN5','LN-normal',
                                                          ifelse(HNSCC_myeloid$orig.ident == 'NT1'|HNSCC_myeloid$orig.ident == 'NT2'|HNSCC_myeloid$orig.ident == 'NT3','NT',
                                                                 ifelse(HNSCC_myeloid$orig.ident == 'Pre-Ca1'|HNSCC_myeloid$orig.ident == 'Pre-Ca2'|HNSCC_myeloid$orig.ident == 'Pre-Ca3','Pre',
                                                                        ifelse(HNSCC_myeloid$orig.ident == 'R-OSCC1'|HNSCC_myeloid$orig.ident == 'R-OSCC2'|
                                                                                 HNSCC_myeloid$orig.ident == 'R-OSCC3'|HNSCC_myeloid$orig.ident == 'R-OSCC4','R','NO'))))))))

HNSCC_myeloidSub <- subset(HNSCC_myeloid,subset = stage %in% c('NT','Pre','E','A','R'))
HNSCC_myeloidSub$stage <- factor(HNSCC_myeloidSub$stage,levels = c('NT','Pre','E','A','R'))

##比例变化
subF <- HNSCC_myeloidSub@meta.data
SufibrosDisA <- table(subF[,c("stage","DefineTypes")]) %>% 
  data.frame %>% set_colnames(c("Stage","CellTypes","Number"))

SufibrosDisA_Tissue <- lapply(split(SufibrosDisA,SufibrosDisA$Stage),function(X){
  X%>% dplyr::mutate(Per=100*Number/sum(Number))
}) %>% dplyr::bind_rows(.)

SufibrosDisA_Tissue$Stage <- factor(SufibrosDisA_Tissue$Stage,levels = c('NT','Pre','E','A','R'))
SufibrosDisA_Tissue$CellTypes <- factor(SufibrosDisA_Tissue$CellTypes,levels = names(table(HNSCC_myeloid$DefineTypes)))

SufibrosDisA_TissueP <- arrange(SufibrosDisA_Tissue,CellTypes) %>% dplyr::filter(Number != 0 )
openxlsx::write.xlsx(SufibrosDisA_TissueP,file=paste(OutPath,"/",folder,"/Fig3G.xlsx",sep=""))

## F
p <- DotPlot(HNSCC_myeloid,feature = c('CLEC9A','BATF3','CD1C','FCER1A','LILRA4','GZMB','LAMP3','FSCN1','FCN1','VCAN','TPSB2','CMA1',
                                     'CXCL10','CXCL9','C1QC','MRC1','SPP1','CD68','FOLR2','CCL18','MKI67','TOP2A'),cols = 'Spectral')+
     theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1)) 
Plot <- p$data

openxlsx::write.xlsx(Plot,file=paste(OutPath,"/",folder,"/Fig3F.xlsx",sep=""))

## J
CompareGroup <- readRDS('/work/smy/Project/HNSCC_26sample/2.data/SPP1 vs FOLR2_GO_KEGGenrichmentResult.rds')
p <- dotplot(CompareGroup$GO_Compare$V1, showCategory = 10)+
  scale_color_gradientn(colors= colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"),space="rgb")(50))

Plot <- p$data
openxlsx::write.xlsx(Plot,file=paste(OutPath,"/",folder,"/Fig3J.xlsx",sep=""))
