OutPath <- "/work/brj/Collaboration/2022/scRNA/HNSCC/Revise3/"
folder <- "Sfig6"

## b /work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Code/26.inferCNV_DE.R
library(data.table)
WorkPath <- "/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Data/inferCNV/"
finalDEG <- fread(paste0(WorkPath,"DEG_CNV.txt"))
finalUp <- finalDEG[finalDEG$class == "Up",]
finalDn <- finalDEG[finalDEG$class == "Down",]

## LB DEG 
LBgene <- fread("/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Data/LB/Gene_malignant-Sample_RvsA.txt",data.table = F)

LBgene$class <- ifelse(LBgene$p_val < 0.05 & LBgene$avg_log2FC > 0.25,"Up",ifelse(LBgene$p_val < 0.05 & LBgene$avg_log2FC < -0.25,"Down","None"))
LBUp <- LBgene[LBgene$class == "Up",]
LBDn <- LBgene[LBgene$class == "Down",]
CNV_UP = finalUp$symbol
scRNA_UP = LBUp$Gene

length_diff <- abs(length(CNV_UP) - length(scRNA_UP))
if(length(CNV_UP) < length(scRNA_UP)) {
  CNV_UP <- c(CNV_UP, rep(NA, length_diff))
} else if(length(CNV_UP) > length(scRNA_UP)) { # 如果scRNA_UP较短
  scRNA_UP <- c(scRNA_UP, rep(NA, length_diff))
}
Up <- data.frame(CNV_UP, scRNA_UP)
openxlsx::write.xlsx(Up,file=paste(OutPath,"/",folder,"/sFig6bUp.xlsx",sep=""))

CNV_DOWN = finalDn$symbol
scRNA_DOWN = LBDn$Gene

length_diff <- abs(length(CNV_DOWN) - length(scRNA_DOWN))
if(length(CNV_DOWN) < length(scRNA_DOWN)) {
  CNV_DOWN <- c(CNV_DOWN, rep(NA, length_diff))
} else if(length(CNV_DOWN) > length(scRNA_DOWN)) { # 如果scRNA_DOWN较短
  scRNA_DOWN <- c(scRNA_DOWN, rep(NA, length_diff))
}
DOWN <- data.frame(CNV_DOWN, scRNA_DOWN)
openxlsx::write.xlsx(DOWN,file=paste(OutPath,"/",folder,"/sFig6bDOWN.xlsx",sep=""))



## d GSE234933
#/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Code/80.GSE234933_OC_statistic.R
## GSE234933 TUMOR
PRtumors <- readr::read_rds(paste0("/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Fig/Recurrent/GSE234933/OC/FeaturePlot/","/PRtumorsOC_Sele.rds.gz"))
Sample <- PRtumors@meta.data %>% tibble::rownames_to_column("Sample") %>% dplyr::select(c(Sample,Tissue)) %>% set_colnames(c("Sample","Group"))

## fanjia cell proliferation
PRtumors_Proliferation <- readr::read_rds("/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Fig/Recurrent/GSE234933/OC/Score/ssGSEA_PR_Proliferation.rds.gz")
ProExp <- t(PRtumors_Proliferation) %>% as.data.frame %>% tibble::rownames_to_column("Sample")

PRtumors_KEGG <- readr::read_rds("/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Fig/Recurrent/GSE234933/OC/Score/ssGSEA_PR_OC_KEGG.rds.gz")
Candidatapwy <- c("Glycolysis / Gluconeogenesis","Oxidative phosphorylation","Cell cycle","Antigen processing and presentation")
sobjlists <- t(PRtumors_KEGG) %>% as.data.frame %>% dplyr::select(Candidatapwy) %>%
  tibble::rownames_to_column("Sample") %>% inner_join(Sample) %>%
  inner_join(ProExp) %>%
  dplyr::select(c(Sample,Group,Fanjia_Proliferation,'Oxidative phosphorylation'))%>%
  dplyr::rename("Cell proliferation" = "Fanjia_Proliferation")
openxlsx::write.xlsx(sobjlists,file=paste(OutPath,"/",folder,"/sFig6d.xlsx",sep=""))

pp <- lapply(c(3:4), function(sub){
  p1 <- ggplot(sobjlists,aes(x= Group, y = sobjlists[,sub])) + 
    geom_violin(aes(fill=Group)) +     
    geom_boxplot(aes(fill=Group),width=0.1)+#可用于将中位数点添加到箱线图中
    labs(y= as.character(names(sobjlists)[sub]), x = "Group") + 
    scale_fill_manual(limits=c("P","R"),values=c("#4072B5","#2CAE65"),guide=F)+
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5))+
    ggpubr::stat_compare_means(comparisons = list(c("P","R")), method = "wilcox.test",na.rm = T) # Add pairwise comparisons p-value
  return(p1)
})

## e  GSE173955
Vpath2 <- "/server1_work//brj/Collaboration/scRNA/2022/HNSCC/Revise/Fig/GEO/"
Markers_GSE173855 <- fread(paste0(Vpath2,"/GSE173855_OC_DEG.txt"),data.table = F)
names(Markers_GSE173855)[1] <- "gene"
subRank <- Markers_GSE173855 %>% dplyr::arrange(desc(logFC))
subRankList <- subRank$logFC
names(subRankList) <- subRank$gene
subRankList  <- rev(sort(subRankList))

MetabolismPathways <- read.gmt("/server1_work/yye/Project/Data/Public/GeneList/KEGG_metabolism.gmt")
NonMetabolismPathways <- read.gmt("/server1_work/yye/Project/Data/Public/GeneList/nonMetabolic_KEGG.gmt")
KEGG <- rbind(MetabolismPathways,NonMetabolismPathways)


fgseaRes_KEGG <- clusterProfiler::GSEA(geneList = subRankList,
                                       nPerm=100,TERM2GENE=KEGG,
                                       minGSSize=5,
                                       maxGSSize=1000,
                                       pvalueCutoff = 1)
GSEAR_KEGG <- fgseaRes_KEGG@result 
openxlsx::write.xlsx(GSEAR_KEGG, paste0(Opath,"GSE173855_GSEAresult_KEGG.xlsx"))


candidateGeneSets <- c("Glycolysis / Gluconeogenesis","Oxidative phosphorylation","Cell cycle","Antigen processing and presentation")

NES_can <- signif(fgseaRes_KEGG[fgseaRes_KEGG$ID %in% candidateGeneSets,c("NES","pvalue")],digits = 2)

library(enrichplot)
GSEAplot <- lapply(candidateGeneSets,function(candidateGeneSet){
  NES_can <- signif(fgseaRes_KEGG[fgseaRes_KEGG$ID %in% candidateGeneSet,c("NES","pvalue")],digits = 2)
  
  p <- gseaplot2(fgseaRes_KEGG,##KEGG gse的对象zw
                 geneSetID = candidateGeneSet,color = "green",#TcellRelated$ID,
                 pvalue_table=F,subplots = 1:2)+
    annotate("text",x=0.5,y=0.8,label = paste(candidateGeneSet,"\nNES = ",NES_can[1],"; p = ",NES_can[2],sep=""),vjust=1)
  return(p)
})

## f 
#/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Code/87.GSE234933_cellphonedb_Plot.R
cpdbresult
openxlsx::write.xlsx(cpdbresult,file=paste(OutPath,"/",folder,"/sFig6f.xlsx",sep=""))

## g
PRtumors <- readr::read_rds(paste0("/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Fig/Recurrent/GSE234933/OC/FeaturePlot/","/PRtumorsOC_Sele.rds.gz"))
Idents(PRtumors) <- "Tissue"
p1 <- VlnPlot(PRtumors,features = "CDKN2A",split.by = "Tissue")
p2 <- VlnPlot(PRtumors,features = "EGFR",split.by = "Tissue")
p3 <- VlnPlot(PRtumors,features = "VEGFA",split.by = "Tissue")
p4 <- VlnPlot(PRtumors,features = 'TGFB1',split.by = "Tissue")

PP <- p1$data %>% bind_cols(p2$data[,1]) %>% bind_cols(p3$data[,1])%>%bind_cols(p4$data[,1]) %>% 
      dplyr::select(c(ident,split),everything())
names(PP)[4:6] <- c("EGFR","VEGFA",'TGFB1')
openxlsx::write.xlsx(PP,file=paste(OutPath,"/",folder,"/sFig6g.xlsx",sep=""),row.names = T)
