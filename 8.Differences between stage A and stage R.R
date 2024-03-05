
#.libPaths(c('~/miniconda3/envs/R4.1/lib/R/library',"/home/smy/R/x86_64-pc-linux-gnu-library/4.1","/usr/local/lib64/R/library"))
library(rlang)
library(infercnv)
library(Seurat)
library(future)
setwd('/work/brj/Collaboration/2022/scRNA/HNSCC/inferCNV/')

HNSCC_Whole <- readRDS('/work/smy/Project/HNSCC/2.data/DeinfeTypes/HNSCC_26sampleAll_DefineTypes.rds.gz')
Idents(HNSCC_Whole) <- 'Maincluster'
DimPlot(HNSCC_Whole)
HNSCC_epi <- subset(HNSCC_Whole,subset = Maincluster ==  'Epithelial cells' & stage %in% c("A","R","NT")) 


dfcount = as.data.frame(HNSCC_epi@assays$RNA@counts)


groupinfo= data.frame(cellId = colnames(dfcount))
groupinfo$stage = HNSCC_epi$stage

library(AnnoProbe)
geneInfor=annoGene(rownames(dfcount),"SYMBOL",'human')
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]


dfcount =dfcount [rownames(dfcount ) %in% geneInfor[,1],]
dfcount =dfcount [match( geneInfor[,1], rownames(dfcount) ),] 


expFile='expFile.txt'
write.table(dfcount ,file = expFile,sep = '\t',quote = F)
groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)


infercnv_obj = CreateInfercnvObject(raw_counts_matrix='expFile.txt',
                                    annotations_file='groupFiles.txt',
                                    delim="\t",
                                    gene_order_file='geneFile.txt',
                                    ref_group_names = c("NT"))
future::plan("multiprocess",workers=12)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir1",  # 输出文件夹
                             cluster_by_groups=T,   # 聚类
                             denoise=T, #去噪
                             HMM=T,# 是否基于HMM预测CNV
                             BayesMaxPNormal = 0,
                             num_threads=10) 

readr::write_rds(infercnv_obj,'infercnv_obj',compress = 'gz')

expr <- read.table("/work/brj/Collaboration/2022/scRNA/HNSCC/inferCNV/output_dir1/infercnv.observations.txt", header=T) %>% as.matrix()
expr.scale <- scale(t(expr))
tmp1 <- sweep(expr.scale, 2, apply(expr.scale, 2, min),'-')
tmp2 <- apply(expr.scale, 2, max) - apply(expr.scale,2,min)
expr_1 <- t(2*sweep(tmp1, 2, tmp2, "/")-1)

cnv_score <- as.data.frame(colSums(expr_1 * expr_1))
colnames(cnv_score)="cnv_score"
cnv_score <- tibble::rownames_to_column(cnv_score, var='cell')
cnv_score$cell <- gsub('[.]','-',cnv_score$cell) #不包含参考的样本在内

HNSCC_Whole <- readRDS('/work/smy/Project/HNSCC/2.data/DeinfeTypes/HNSCC_26sampleAll_DefineTypes.rds.gz')
Idents(HNSCC_Whole) <- 'Maincluster'
#DimPlot(HNSCC_Whole)
HNSCC_epi <- subset(HNSCC_Whole,subset = Maincluster ==  'Epithelial cells' & stage %in% c("A","R")) # 不包括参考样本

HNSCC_epi$cnv_score <- cnv_score[match(rownames(HNSCC_epi@meta.data),cnv_score$cell),2]
HNSCC_epi@meta.data[grep('TRUE',is.na(HNSCC_epi$cnv_score)),]
cell <- rownames(HNSCC_epi@meta.data[grep('FALSE',is.na(HNSCC_epi$cnv_score)),])
HNSCC_epiF <- subset(HNSCC_epi,cells = cell)


#value <- as.numeric(HNSCC_epiF$cnv_score)
#thresthold <- median(value) - 1 * sd(value) # 这里是一个标准差

#meta %>% mutate(condition=if_else(.$cnv_score <= thresthold,'non-malignant','malignant')) -> TumorEpi@meta.data # 分出了恶性和非恶性

#HNSCC_epiF$cnv_group <- ifelse(HNSCC_epiF$cnv_score <=thresthold,'non-malignant','malignant')
#Idents(HNSCC_epiF) <- 'cnv_group'
#table(HNSCC_epiF$cnv_group) 

#HNSCC_epiMalignant <- subset(HNSCC_epiF,idents = 'malignant')
#table(HNSCC_epiMalignant$stage)

#
#Idents(HNSCC_epi) <- HNSCC_epi$stage
#HNSCC_epi_non_Malignant <- subset(HNSCC_epi,idents = c('NT','LN-normal'))
#HNSCC_epi_non_Malignant$cnv_group <- 'non-malignant'

#HNSCC_epi_inferCNV <- merge(HNSCC_epiMalignant,HNSCC_epi_non_Malignant)
#readr::write_rds(HNSCC_epi_inferCNV,'/work/smy/Project/HNSCC_26sample/2.data/HNSCC_epi_inferCNV.rds')

#violin plot =======================================================
sobjlists = FetchData(object = HNSCC_epiF, vars = c("stage","cnv_score"))
openxlsx::write.xlsx(sobjlists,file=paste(OutPath,"/",folder,"/Fig6A.xlsx",sep=""))

my_comparisons <- list(c("A","R"))
# A: "#E78AC3" R:"#E39A35"
p1 <- ggplot(sobjlists,aes(x= stage, y = cnv_score)) + 
      geom_violin(aes(fill=stage)) + #scale_color_manual(values = c("#E78AC3","#E39A35")) +  
      scale_fill_manual(limits=c("A","R"),
                    values= c("#E78AC3","#58A4C3"),name="")+
      stat_summary(fun.y = median, geom = "point", shape = 23, size=4)+#可用于将中位数点添加到箱线图中
      labs(y= "cnv_score", x = "stage") + 
      theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5))+
      ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",na.rm = T)

pdf("/work/brj/Collaboration/2022/scRNA/HNSCC/inferCNV/CNV_stageRAviolin.pdf",width = 5,height = 5)
p1
dev.off()

expr <- read.table("/work/brj/Collaboration/2022/scRNA/HNSCC/inferCNV//output_dir1/infercnv.observations.txt", header=T) %>% as.matrix()

Candigenes <- c("CCND1","CDKN2A","EGFR","MYC","VEGFA","TGFB1")
Expr <- expr %>% as.data.frame() %>%  dplyr::filter(rownames(.) %in% Candigenes) %>% 
        t() %>% as.data.frame() %>% tibble::rownames_to_column("cell") 
Expr$cell <- gsub('[.]','-',Expr$cell) 

### 
HNSCC_Whole <- readRDS('/work/smy/Project/HNSCC/2.data/DeinfeTypes/HNSCC_26sampleAll_DefineTypes.rds.gz')
HNSCC_epi <- subset(HNSCC_Whole,subset = Maincluster ==  'Epithelial cells' & stage %in% c("A","R")) 
HNSCC_epi$CCND1_cnv <- Expr[match(rownames(HNSCC_epi@meta.data),Expr$cell),"CCND1"]
HNSCC_epi$CDKN2A_cnv <- Expr[match(rownames(HNSCC_epi@meta.data),Expr$cell),"CDKN2A"]
HNSCC_epi$EGFR_cnv <- Expr[match(rownames(HNSCC_epi@meta.data),Expr$cell),"EGFR"]
HNSCC_epi$MYC_cnv <- Expr[match(rownames(HNSCC_epi@meta.data),Expr$cell),"MYC"]
HNSCC_epi$VEGFA_cnv <- Expr[match(rownames(HNSCC_epi@meta.data),Expr$cell),"VEGFA"]
HNSCC_epi$TGFB1_cnv <- Expr[match(rownames(HNSCC_epi@meta.data),Expr$cell),"TGFB1"]


sobjlists = FetchData(object = HNSCC_epi, vars = c("stage","CCND1_cnv","CDKN2A_cnv","EGFR_cnv","MYC_cnv","VEGFA_cnv","TGFB1_cnv"))
my_comparisons <- list(c("A","R"))
# A: "#E78AC3" R:"#E39A35"
pp <- lapply(c(2:7), function(sub){
   p1 <- ggplot(sobjlists,aes(x= stage, y = sobjlists[,sub])) + 
      geom_violin(aes(fill=stage)) + 
      scale_fill_manual(limits=c("A","R"),
                        values= c("#E78AC3","#58A4C3"),name="")+
      stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)+#可用于将中位数点添加到箱线图中
      labs(y= as.character(names(sobjlists)[sub]), x = "stage") + 
      theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5))+
      ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",na.rm = T) # Add pairwise comparisons p-value
   return(p1)
})

pdf("/work/brj/Collaboration/2022/scRNA/HNSCC/inferCNV/CNV_6Genes_Vln.pdf",width = 15,height = 10)
cowplot::plot_grid(plotlist = pp,ncol = 3)
dev.off()

expr <- read.table("/work/brj/Collaboration/2022/scRNA/HNSCC/inferCNV//output_dir1/infercnv.observations.txt", header=T) %>% as.matrix()
### 所有的值-1后取绝对值
Candigenes <- c("CCND1","CDKN2A","EGFR","MYC","VEGFA","TGFB1")
### 减1
expr.impute <- sweep(expr, 2, 1,'-')
### abs
expr.abs <- apply(expr.impute, 2, abs)

Expr <- expr.abs %>% as.data.frame() %>%  dplyr::filter(rownames(.) %in% Candigenes) %>% 
        t() %>% as.data.frame() %>% tibble::rownames_to_column("cell") 
Expr$cell <- gsub('[.]','-',Expr$cell) 

### 
HNSCC_Whole <- readRDS('/work/smy/Project/HNSCC/2.data/DeinfeTypes/HNSCC_26sampleAll_DefineTypes.rds.gz')
HNSCC_epi <- subset(HNSCC_Whole,subset = Maincluster ==  'Epithelial cells' & stage %in% c("A","R")) 
HNSCC_epi$CCND1_cnv <- Expr[match(rownames(HNSCC_epi@meta.data),Expr$cell),"CCND1"]
HNSCC_epi$CDKN2A_cnv <- Expr[match(rownames(HNSCC_epi@meta.data),Expr$cell),"CDKN2A"]
HNSCC_epi$EGFR_cnv <- Expr[match(rownames(HNSCC_epi@meta.data),Expr$cell),"EGFR"]
HNSCC_epi$MYC_cnv <- Expr[match(rownames(HNSCC_epi@meta.data),Expr$cell),"MYC"]
HNSCC_epi$VEGFA_cnv <- Expr[match(rownames(HNSCC_epi@meta.data),Expr$cell),"VEGFA"]
HNSCC_epi$TGFB1_cnv <- Expr[match(rownames(HNSCC_epi@meta.data),Expr$cell),"TGFB1"]


sobjlists = FetchData(object = HNSCC_epi, vars = c("stage","CCND1_cnv","CDKN2A_cnv","EGFR_cnv","MYC_cnv","VEGFA_cnv","TGFB1_cnv"))
my_comparisons <- list(c("A","R"))
# A: "#E78AC3" R:"#E39A35"
pp <- lapply(c(2:7), function(sub){
   p1 <- ggplot(sobjlists,aes(x= stage, y = sobjlists[,sub])) + 
      geom_violin(aes(fill=stage)) + 
      scale_fill_manual(limits=c("A","R"),
                        values= c("#E78AC3","#58A4C3"),name="")+
      stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)+#可用于将中位数点添加到箱线图中
      labs(y= as.character(names(sobjlists)[sub]), x = "stage") + 
      theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5))+
      ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",na.rm = T) # Add pairwise comparisons p-value
   return(p1)
})




pdf("/work/brj/Collaboration/2022/scRNA/HNSCC/inferCNV/CNV_6Genes_scale_Vln.pdf",width = 15,height = 10)
cowplot::plot_grid(plotlist = pp,ncol = 3)
dev.off()



expr <- read.table("/work/brj/Collaboration/2022/scRNA/HNSCC/inferCNV//output_dir1/infercnv.observations.txt", header=T) %>% as.matrix()
### 归一化到1
Candigenes <- c("CCND1","CDKN2A","EGFR","MYC","VEGFA","TGFB1")
expr.scale <- scale(t(expr))
tmp1 <- sweep(expr.scale, 2, apply(expr.scale, 2, min),'-')
tmp2 <- apply(expr.scale, 2, max) - apply(expr.scale,2,min)
expr_1 <- t(2*sweep(tmp1, 2, tmp2, "/")-1)

Expr <- expr_1 %>% as.data.frame() %>%  dplyr::filter(rownames(.) %in% Candigenes) %>% 
   t() %>% as.data.frame() %>% tibble::rownames_to_column("cell") 
Expr$cell <- gsub('[.]','-',Expr$cell) #不包含参考的样本在内

### 
HNSCC_Whole <- readRDS('/work/smy/Project/HNSCC/2.data/DeinfeTypes/HNSCC_26sampleAll_DefineTypes.rds.gz')
HNSCC_epi <- subset(HNSCC_Whole,subset = Maincluster ==  'Epithelial cells' & stage %in% c("A","R")) 
HNSCC_epi$CCND1_cnv <- Expr[match(rownames(HNSCC_epi@meta.data),Expr$cell),"CCND1"]
HNSCC_epi$CDKN2A_cnv <- Expr[match(rownames(HNSCC_epi@meta.data),Expr$cell),"CDKN2A"]
HNSCC_epi$EGFR_cnv <- Expr[match(rownames(HNSCC_epi@meta.data),Expr$cell),"EGFR"]
HNSCC_epi$MYC_cnv <- Expr[match(rownames(HNSCC_epi@meta.data),Expr$cell),"MYC"]
HNSCC_epi$VEGFA_cnv <- Expr[match(rownames(HNSCC_epi@meta.data),Expr$cell),"VEGFA"]
HNSCC_epi$TGFB1_cnv <- Expr[match(rownames(HNSCC_epi@meta.data),Expr$cell),"TGFB1"]


sobjlists = FetchData(object = HNSCC_epi, vars = c("stage","CCND1_cnv","CDKN2A_cnv","EGFR_cnv","MYC_cnv","VEGFA_cnv","TGFB1_cnv"))
my_comparisons <- list(c("A","R"))
# A: "#E78AC3" R:"#E39A35"
pp <- lapply(c(2:7), function(sub){
   p1 <- ggplot(sobjlists,aes(x= stage, y = sobjlists[,sub])) + 
      geom_violin(aes(fill=stage)) + 
      scale_fill_manual(limits=c("A","R"),
                        values= c("#E78AC3","#58A4C3"),name="")+
      stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)+#可用于将中位数点添加到箱线图中
      labs(y= as.character(names(sobjlists)[sub]), x = "stage") + 
      theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5))+
      ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",na.rm = T) # Add pairwise comparisons p-value
   return(p1)
})

pdf("/work/brj/Collaboration/2022/scRNA/HNSCC/inferCNV/CNV_6Genes_scale1_Vln.pdf",width = 15,height = 10)
cowplot::plot_grid(plotlist = pp,ncol = 3)
dev.off()
