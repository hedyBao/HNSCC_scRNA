##Figure 5

##LN-out LN-in上皮细胞差异基因======
HNSCC_epi_tumor <- readRDS('/work/smy/Project/HNSCC_26sample/2.data/26sample_tumorcells_Harmony_immuneReference_All.rds')

DimPlot(HNSCC_epi_tumor)
HNSCC_epi_tumor$stage <- as.character(HNSCC_epi_tumor$stage)
HNSCC_epi_tumor$stage[grep('LN1|LN6|LN7',HNSCC_epi_tumor$orig.ident)] <- 'LN-out'
HNSCC_epi_tumor$stage[grep('LN2|LN3|LN8',HNSCC_epi_tumor$orig.ident)] <- 'LN-in'
VlnPlot(HNSCC_epi_tumor,features = c('APOE','FN1'))

Idents(HNSCC_epi_tumor) <- HNSCC_epi_tumor$stage
Epi <- subset(HNSCC_epi_tumor,idents = c('LN-in','LN-out'))

Markers <- FindAllMarkers(Epi,only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

Markers$pvalue_L <- -log10(Markers$p_val+10^-200)
Markers$Class <- ifelse( Markers$avg_log2FC > 0.25 & Markers$p_val < 0.05 & Markers$cluster == 'LN-out',"Up",
                         ifelse(Markers$avg_log2FC > 0.25 & Markers$p_val < 0.05 & Markers$cluster == 'LN-in' ,"Down","None" ))

Markers$avg_log2FC <- ifelse(Markers$cluster == 'LN-in',-Markers$avg_log2FC,Markers$avg_log2FC)
Gene_Labels <- Markers %>% dplyr::filter(Class %in% c("Up","Down")) %>% dplyr::group_by(Class)  %>% 
  dplyr::top_n(10,abs(avg_log2FC))

outpath <- "/work/brj/Collaboration/2022/scRNA/HNSCC/Result/"
openxlsx::write.xlsx(Markers,paste0(outpath,"Fig5/","LN-in vs LN-out diff genes.xlsx"))


PP <- ggplot(Markers,aes(x=avg_log2FC,y= pvalue_L))+
  geom_point(aes(color=Class),size=0.8) +
  # geom_point(color=NA)+
  labs(x="log2(Foldchange)",y="-log10(p-value)")+
  scale_color_manual(limits=c("Up","Down"),values=c("#e41a1c","#3182bd"),#values=c(NA,NA,NA),#
                     labels=c(paste("LN-out (n=",table(Markers$Class)['Up'],")",sep=""),
                              paste("LN-in (n=",table(Markers$Class)['Down'],")",sep="")
                              ),name="")+
  geom_vline(xintercept = c(-1,1),linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme(panel.background = element_rect(fill=NA,color="black"),axis.title= element_text(color="black",size=12),panel.grid.major = element_blank(),
        axis.text= element_text(color="black",size=10),legend.position = "bottom",
        legend.background = element_blank(),legend.key = element_blank())+
  #ggrepel::geom_text_repel(data =Markers[Markers$Gene %in% TF$V1 & Markers$Class !="None",],aes(label=Gene),force=1,color="red")+
  ggrepel::geom_text_repel(data =Gene_Labels,aes(label=gene),force=1,color="black")+  
  geom_point(data =Markers[Markers$gene %in% c(Gene_Labels$gene) & Markers$Class !="None",],shape=1)

outpath <- "/work/brj/Collaboration/2022/scRNA/HNSCC/Result/"
pdf(paste0(outpath,"Fig5/","Epithial_0_1_diffGeneVolcano.pdf"),width = 6,height = 6)
print(PP)
dev.off()

##LN-out LN-in上皮细胞通路富集======
GO_BP_Fun <- function(CompareGroup,groups=c("LN-in","LN-out"),OrgData="org.Hs.eg.db",pcut=0.05,qcut=0.2){
  library(clusterProfiler)
  if(OrgData=="org.Hs.eg.db"){library(org.Hs.eg.db)}else{library(org.Hs.eg.db)}
  CompareGroup <- CompareGroup %>% dplyr::mutate(GO_Compare=purrr::map(.x=CompareResult,function(.x){
    SigSub <- .x[.x$Class %in% groups,] %>% as.data.frame()
    if(length(unique(SigSub$gene))>10){
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
KEGG_Fun <- function(CompareGroup,groups=c("LN-in","LN-out"),OrgSp="hsa",pcut=0.05,qcut=0.2){
  CompareGroup <- CompareGroup %>% dplyr::mutate(KEGG_Compare= purrr::map(.x=CompareResult,function(.x){
    
    SigSub <- .x[.x$Class %in% c("LN-out","LN-in"),] %>% as.data.frame()
    genelist <- clusterProfiler::bitr(SigSub$gene,fromType = "SYMBOL",
                                      toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
    
    SigSub <- merge( SigSub,genelist,by.x="gene",by.y="SYMBOL")
    if(length(unique(SigSub$gene))>10){
      formula_res <- clusterProfiler::compareCluster(
        #keyType="SYMBOL",
        ENTREZID~Class,
        data=SigSub,
        organism=OrgSp,
        fun="enrichKEGG", 
        # OrgDb="org.Mm.eg.db",
        pAdjustMethod = "BH",
        pvalueCutoff  = pcut,
        qvalueCutoff  = qcut
      )
      
      return(formula_res)
    }else{return(NA)}
  })) 
  
}

Markers <- FindAllMarkers(Epi,only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
Markers$Class <- ifelse(Markers$avg_log2FC > 0.25 & Markers$cluster == 'LN-in',"LN-in",
                        ifelse(Markers$avg_log2FC > 0.25 & Markers$cluster == 'LN-out',"LN-out","None" ))
table(Markers$Class)
CompareList = c("LN-in vs LN-out")
CompareGroup <- tibble::tibble(CompareList=c("LN-in vs LN-out"),CompareResult = list(V1=Markers))

#GO
CompareGroup <- GO_BP_Fun(CompareGroup,groups=c("LN-in","LN-out"),OrgData="org.Hs.eg.db",pcut=0.05,qcut=0.2)
GOresult <- CompareGroup$GO_Compare$V1@compareClusterResult
readr::write_csv(GOresult,paste0(outpath,"Fig5/","LNin vs LNout_GOenrichmentResult.csv"))

pdf(paste0(outpath,"Fig5/","LNin vs LNout_GO.pdf"),width = 6,height = 8)
dotplot(CompareGroup$GO_Compare$V1, showCategory = 10)+
  scale_color_gradientn(colors= colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"),space="rgb")(50))
dev.off()



### 不画共有的
GOresult <- CompareGroup$GO_Compare$V1@compareClusterResult
INres <- GOresult %>% dplyr::filter(Cluster == "LN-in") %>% dplyr::select(Description)
Outres <- GOresult %>% dplyr::filter(Cluster == "LN-out") %>% dplyr::select(Description)
SpecificIN <- setdiff(INres$Description,Outres$Description) # 取第一个参数唯一的
SpecificOuts <- setdiff(Outres$Description,INres$Description) # 取第一个参数唯一的
Virus <- grep("virus|viral",SpecificOuts,value = T)
SpecificOut <- setdiff(SpecificOuts,Virus)

topN <- 10
GOresults <- GOresult %>% dplyr::mutate(Ratio = round(Count/584,4)) %>% 
  dplyr::select(c(Cluster,Description,Ratio,p.adjust)) %>%
  dplyr::arrange(Cluster,p.adjust) #%>% dplyr::group_by(Cluster) %>% dplyr::top_n(-topN,p.adjust)  
GOresultsOut <- GOresults %>% dplyr::filter(Description %in% SpecificOut) %>% dplyr::top_n(-topN,p.adjust) 
GOresultsIn <- GOresults %>% dplyr::filter(Description %in% SpecificIN) %>% dplyr::top_n(-topN,p.adjust) 


DE_DescriptionFun_Groups_BPAll <- bind_rows(GOresultsOut,GOresultsIn)

TT2 <- data.frame(Cluster=rep(unique(DE_DescriptionFun_Groups_BPAll$Cluster),times=length(unique(DE_DescriptionFun_Groups_BPAll$Description))),
                  Description=rep(unique(DE_DescriptionFun_Groups_BPAll$Description),each=length(unique(DE_DescriptionFun_Groups_BPAll$Cluster))))

#图片处理 美观
DE_DescriptionFun_Groups_BPAllF <- DE_DescriptionFun_Groups_BPAll %>% dplyr::select(Description,Cluster,Ratio) %>%
  tidyr::spread(Cluster,Ratio,fill=0)
DE_DescriptionFun_Groups_BPAllF <- DE_DescriptionFun_Groups_BPAllF[do.call(order,DE_DescriptionFun_Groups_BPAllF[,2:ncol(DE_DescriptionFun_Groups_BPAllF)]),]


#DE_DescriptionFun_Groups_BPAll$Ratio <- ifelse(DE_DescriptionFun_Groups_BPAll$Ratio > 1,1,DE_DescriptionFun_Groups_BPAll$Ratio)

pdf(paste0(outpath,"Fig5/","LNin vs LNout_GOspecific.pdf"),width = 7,height = 6)
P <- ggplot(DE_DescriptionFun_Groups_BPAll,aes(x=Cluster,y=Description))+
  geom_point(aes(color=-log10(p.adjust+10^-27),size=Ratio))+
  labs(color = "-log10(Adjust P)",size="GeneRatio")+
  scale_color_gradientn(limit=c(min(-log10(DE_DescriptionFun_Groups_BPAll$p.adjust)),max(-log10(DE_DescriptionFun_Groups_BPAll$p.adjust))),
                        colors= colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")),space="rgb")(50))+
  theme(panel.background=element_rect(color="black",fill=NA),legend.key =  element_blank(),axis.ticks=element_blank(),
        panel.grid=element_blank(),axis.text.x=element_text(size=12,angle=90,vjust = 0.5,hjust = 1,color="black"),
        axis.text.y=element_text(size=15,color="black"),axis.title=element_blank(),legend.position = "bottom")+
  scale_y_discrete(limit=(DE_DescriptionFun_Groups_BPAllF$Description))+
  geom_tile(data=TT2,aes(x=Cluster,y=Description),fill=NA,color="lightgray")
dev.off()

#Kegg
CompareGroup <- KEGG_Fun(CompareGroup,groups=c("LN-in","LN-out"),OrgSp="hsa",pcut=0.05,qcut=0.2)
KEGGresult <- CompareGroup$KEGG_Compare$V1@compareClusterResult
readr::write_csv(KEGGresult,paste0(outpath,"Fig5/","LNin vs LNout_KEGGenrichmentResult.csv"))

dotplot(CompareGroup$KEGG_Compare$V1, showCategory = 10)+
  scale_color_gradientn(colors= colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"),space="rgb")(50))

##cellchat互作强度======
#上皮细胞只保留肿瘤细胞
#1.取出LN的样本
#2.上皮细胞只保留肿瘤细胞
#3.分stage跑cellchat
Idents(HNSCC_Whole) <- HNSCC_Whole$stage
HNSCC_LN <- subset(HNSCC_Whole,idents = c('LN-in','LN-out','LN-normal'))

HNSCC_LN$Subcluster_RE <- HNSCC_LN$Subcluster
Tumorcells <- rownames(HNSCC_epi_tumor@meta.data)
HNSCC_LN$Subcluster_RE[Tumorcells] <- 'Tumor cells'
unique(HNSCC_LN$Subcluster_RE)

HNSCC_LN_split <- SplitObject(HNSCC_LN,split.by = 'stage')


HNSCC_LN_cellchat <- lapply(names(HNSCC_LN_split), function(name){
  sub <- HNSCC_LN_split[[name]]
  ##cellchat
  cellchat <- createCellChat(object = sub,group.by = 'Subcluster_RE')
  
  #数据库信息
  CellChatDB <- CellChatDB.human
  
  ##pre-process
  cellchat@DB <- CellChatDB
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  #cellchat@LR$LRsig
  cellchat <- projectData(cellchat,PPI.human)
  
  ##infer interaction
  ##liagnd receptor
  cellchat <- computeCommunProb(cellchat,raw.use = FALSE,population.size = TRUE)
  #filter cell<10
  cellchat <- filterCommunication(cellchat,min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  
  cellchat <- computeCommunProbPathway(cellchat)
  df.netp <- subsetCommunication(cellchat,slot.name = 'netP')
  
  cellchat <- aggregateNet(cellchat)
  
  readr::write_rds(cellchat,paste('/work/smy/Project/HNSCC_26sample/2.data/26sample_LNcell_',name,'_cellchat_Subcluster.rds',sep = ''))
})
names(HNSCC_LN_cellchat) <- names(HNSCC_LN_split)
readr::write_rds(HNSCC_LN_cellchat,'26Sample_HNSCC_LN_cellchat.rds')


##figure
HNSCC_LN_cellchat <- readRDS('/work/smy/Project/HNSCC_26sample/2.data/26Sample_HNSCC_LN_cellchat.rds')

NetWeightAll <- lapply(c(1:3), function(sub){
  mat <- HNSCC_LN_cellchat[[sub]]@net$weight
  mat <- as.data.frame(mat)
  mat$stage <- names(HNSCC_LN_cellchat)[sub]
  return(mat)
})

NetWeightAll <- bind_rows(NetWeightAll)
NetWeightAll_melt <- melt(NetWeightAll)
NetWeightAll_melt$cluster <- rep(colnames(NetWeightAll)[1:10],times = 3)
NetWeightAll_melt$cellchat <- paste(NetWeightAll_melt$cluster,'to',NetWeightAll_melt$variable,sep = ' ')

unique(NetWeightAll_melt$cellchat)

NetWeightAll_meltF <- NetWeightAll_melt %>% filter(cluster == 'Tumor cells')
NetWeightAll_meltF$cellchat2 <- paste(NetWeightAll_meltF$cellchat,NetWeightAll_meltF$stage)
NetWeightAll_meltF <- NetWeightAll_meltF[c(7:9,19:21,4:6,1:3,10:12,16:18,22:27),]
NetWeightAll_meltF$cellchat2 <- stringr::str_remove(NetWeightAll_meltF$cellchat2,"Tumor cells to ")

pdf(paste0(outpath,"Fig5/","TumortoSubclusterCellchatLolliplot.pdf"),width = 4,height = 6)
P <- ggplot(NetWeightAll_meltF,aes(x=value,y=cellchat2))+
  ylab('')+
  xlab('weight')+
  #ggtitle(sub)+
  geom_segment(aes(yend=cellchat2),xend=0,colour='grey50')+  ###绘制以数据点为端点的线段
  geom_point(size=3,aes(colour=cellchat))+   ###此处我们将以正负相关(postive  negative)映射其颜色
  scale_colour_brewer(palette = 'Set1')+ ###颜色加深  
  # scale_y_discrete(limits= rev(c('NT','Pre','E','A','R')))+
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_y_discrete(limits= rev(NetWeightAll_meltF$cellchat2))+
  NoLegend()###删除网格线
dev.off()

names(PP) <- unique(NetWeightAll_meltF$cellchat)
cowplot::plot_grid(plotlist = PP,ncol = 4)




### cell chat venn plot olp gene 
#### 提取stage A 特异的信号通路中的基因
FibToMacGeneList <- c("EGF","FGF","GDNF","HGF","IL6","NRG","PARS","PDGF")
MacToFibGeneList <- c("CD23","COMPLEMENT","DESMOSOME","EPHB","FGF","FLT3","GDNF","NGF","NRG","PDGF")
FibToTumorGeneList <- c("CD99","CDH5","CSF3","DESMOSOME","EPO","GDNF","IL4","LIFR","NGF","OCLN","OSM","PDGF","PRL","PSAP","SEMA5","VEGF")
MacToTumorGeneList <- c("ACTIVIN","ADGRE5","CDH5","EPHB","HGF","L1CAM","LIFR","MPZ","MSTN","NPR2","OCLN","TGFb","VEGF")


## cellchat myefib
cellchatMyeFib <- readRDS('/work/brj/Collaboration/2022/scRNA/HNSCC/Data/26sample_definetype_MyeFib_cellchat_stageSplit.rds')
cellchatMyeFib <- cellchatMyeFib[c(6,7,2,1,4,3,5,8)]

A.Fib.net <- subsetCommunication(cellchatMyeFib$A,slot.name = 'net',sources.use = "POSTN+ fibroblast",targets.use = "SPP1+ macorphages" ) %>% 
  dplyr::filter(pathway_name %in% FibToMacGeneList ) %>% dplyr::select(ligand) %>% unique()

A.Mac.net <- subsetCommunication(cellchatMyeFib$A,slot.name = 'net',sources.use = "SPP1+ macorphages",targets.use = "POSTN+ fibroblast" ) %>% 
  dplyr::filter(pathway_name %in% MacToFibGeneList ) %>% dplyr::select(ligand) %>% unique()

## cellchat myefibtumor
cellchatMFT <- readRDS('/work/brj/Collaboration/2022/scRNA/HNSCC/Data/26sample_MyeFibTumor_cellchat_stageSplit2.rds')
cellchatMFT <- cellchatMFT[c(6,7,2,1,4,3,5,8)]

FibToTumor <- subsetCommunication(cellchatMFT$A,slot.name = 'net',sources.use = "POSTN+ fibroblast",targets.use = "TumorCell_malignant") %>% 
  dplyr::filter(pathway_name %in% FibToTumorGeneList) %>% dplyr::select(ligand) %>% unique()

MacToTumor <- subsetCommunication(cellchatMFT$A,slot.name = 'net',sources.use = "SPP1+ macorphages",targets.use = "TumorCell_malignant") %>% 
            dplyr::filter(pathway_name %in% MacToTumorGeneList) %>% dplyr::select(ligand) %>% unique()

TwoOlp <- Reduce(intersect, list(A.Fib.net$ligand, A.Mac.net$ligand)) 
TwoOlpdf <- data.frame(ligand = TwoOlp)

TwoUnion <- Reduce(union, list(A.Fib.net$ligand, A.Mac.net$ligand)) 
TwoUniondf <- data.frame(ligand = TwoUnion)

openxlsx::write.xlsx(TwoUniondf,"/work/brj/Collaboration/2022/scRNA/HNSCC/Result/Fig4/TwoUniondfGeneList.xlsx")




FourUnion <- Reduce(union, list(A.Fib.net$ligand, A.Mac.net$ligand,FibToTumor$ligand,MacToTumor$ligand)) 
Fourdf <- data.frame(ligand = FourUnion)

openxlsx::write.xlsx(Fourdf,"/work/brj/Collaboration/2022/scRNA/HNSCC/Result/Fig4/FourGeneList.xlsx")




