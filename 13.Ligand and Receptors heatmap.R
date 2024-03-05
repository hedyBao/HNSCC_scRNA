### ligand and receptor heatmap
rm(list = ls())
### GSE182227
WorkPath <- "/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Data/GSE182227/NichenetrN/"
nichenet_output <- readr::read_rds(paste0(WorkPath,"1.SenderPOSTN_ReciverSPP1_nichenet_output.rds",sep = ''))
TargetLigand_GSE182227 <- rownames(nichenet_output$ligand_target_matrix)
Targetgenes_GSE182227 <- colnames(nichenet_output$ligand_target_matrix)
Ligands_GSE182227 <- colnames(nichenet_output$ligand_receptor_matrix) # 该矩阵 行是receptors;列是ligands
Receptors_GSE182227 <- rownames(nichenet_output$ligand_receptor_matrix)

### 自己的数据
WorkPath <- "/work/brj/Collaboration/2022/scRNA/HNSCC/Result/Fig4/"
nichenet_output_Ref <- readr::read_rds(paste0(WorkPath,"1.SenderPOSTNFibSPP1Reciver_nichenet_output.rds",sep = ''))
TargetLigand_Ref <- rownames(nichenet_output_Ref$ligand_target_matrix)
Targetgenes_Ref <- colnames(nichenet_output_Ref$ligand_target_matrix)
Receptors_Ref <- rownames(nichenet_output_Ref$ligand_receptor_matrix) # row is receptor
Ligands_Ref <- colnames(nichenet_output_Ref$ligand_receptor_matrix) # column is ligand

### ligand and target genes heatmap : vis_ligand_target 是行为ligand，列为靶基因的矩阵
ligand_target <- nichenet_output$ligand_target_matrix
inter_ligand <- intersect(Ligands_Ref,TargetLigand_GSE182227)#
inter_target <- intersect(Targetgenes_Ref,Targetgenes_GSE182227)#,Targetgenes_GSE182227
inter_receptor <- intersect(Receptors_Ref,Receptors_GSE182227)#,Targetgenes_GSE182227

order_ligand <- c("HAS2","CCL2","IL15","ICAM1","COL18A1","ADAM17","HMGB1","CSF1","GAS6","TGFB1","ANGPT1","CXCL12")


### ligand activity 
ligand_pearson_df <- nichenet_output$ligand_activities[nichenet_output$ligand_activities$test_ligand %in% inter_ligand,c("test_ligand","pearson")] %>% as.data.frame()
vis_ligand_pearson <- ligand_pearson_df %>% dplyr::select(-test_ligand) %>% dplyr::arrange(pearson) %>%
  as.matrix() %>% magrittr::set_rownames(rev(ligand_pearson_df$test_ligand))
p_ligand_pearson = vis_ligand_pearson %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Ligand\nActivity") + theme(legend.text = element_text(size = 9))

library(magrittr)
ileft <- vis_ligand_pearson %>% as.data.frame()
ilefts <- ileft[order_ligand,] %>% as.data.frame() %>% set_colnames("pearson") %>% set_rownames(order_ligand)
openxlsx::write.xlsx(ilefts,file=paste(OutPath,"/",folder,"/Fig4ileft.xlsx",sep=""),row.names = T)


plot_target <- inter_target
order_exp <- ligand_target[rev(order_ligand),plot_target]
p_ligand_target_network = order_exp %>% 
  nichenetr::make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic"))

iright <- ligand_target[order_ligand,plot_target] %>% as.data.frame()
openxlsx::write.xlsx(iright,file=paste(OutPath,"/",folder,"/Fig4iright.xlsx",sep=""),row.names = T)


pdf(paste0(WorkPath,'/1.ligand_activity_target_heatmapPOSTNSPP1.pdf'),width = 12,height = 5)
cowplot::plot_grid(plotlist = list(p_ligand_pearson,p_ligand_target_network),rel_widths = c(1.25,10.75))
dev.off()

### ligand_receptor_heatmap
ligand_receptor <- nichenet_output$ligand_receptor_matrix %>% as.data.frame() %>% dplyr::select(inter_ligand)

p_ligand_receptor_network = ligand_receptor[inter_receptor,] %>% as.matrix() %>% t() %>%
  nichenetr::make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")

j <- ligand_receptor[inter_receptor,] %>% as.data.frame()%>% t() %>% as.data.frame()
openxlsx::write.xlsx(j,file=paste(OutPath,"/",folder,"/Fig4j.xlsx",sep=""),row.names = T)


pdf(paste0(WorkPath,'/2.ligand_receptor_heatmapPOSTNSPP1.pdf'),width = 12,height =  5)
p_ligand_receptor_network
dev.off()


# 4. GO and KEGG pwy 
OurGO <- openxlsx::read.xlsx("/work/brj/Collaboration/2022/scRNA/HNSCC/Result/Fig4/GOtargetPOSTNtoSPP1.xlsx")
OurKEGG <- openxlsx::read.xlsx("/work/brj/Collaboration/2022/scRNA/HNSCC/Result/Fig4/KEGGtargetPOSTNtoSPP1.xlsx")
GSE182227_GO <- openxlsx::read.xlsx("/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Data/GSE182227/NichenetrN/GOtargetPOSTNtoSPP1.xlsx")
GSE182227_KEGG <- openxlsx::read.xlsx("/work/brj/Collaboration/2022/scRNA/HNSCC/Revise/Data/GSE182227/NichenetrN/KEGGtargetPOSTNtoSPP1.xlsx")


GOresult <- GSE182227_GO %>% dplyr::filter(Description %in% OurGO$Description[1:10]) 
GOresults <- GOresult %>% dplyr::mutate(Ratio = round(Count/59,4)) %>% 
  dplyr::select(c(Cluster,Description,Ratio,p.adjust)) %>%
  dplyr::arrange(Cluster,p.adjust) #%>% dplyr::group_by(Cluster) %>% dplyr::top_n(-topN,p.adjust)  
p1 <- ggplot(GOresults,aes(x=Cluster,y=Description))+
  geom_point(aes(color = -log10(p.adjust+10^-27),size=Ratio))+
  labs(color = "-log10(Adjust P)",size="GeneRatio")+
  scale_colour_gradientn(limit=c(0,max(abs(-log10(GOresults$p.adjust+10^-27)))),
                         colors= rev(RColorBrewer::brewer.pal(11, "Spectral")[2:11]))+
  theme_bw()+
  scale_y_discrete(limit=rev(OurGO$Description[1:10]))
openxlsx::write.xlsx(p1$data,file=paste(OutPath,"/",folder,"/Fig4kup.xlsx",sep=""),rownames = T)


### Select pwy 
KEGGresult <- GSE188737_KEGG %>% dplyr::filter(Description %in% OurKEGG$Description[1:10]) 

KEGGresults <- KEGGresult %>% dplyr::mutate(Ratio = round(Count/49,4)) %>% 
  dplyr::select(c(Cluster,Description,Ratio,p.adjust)) %>%
  dplyr::arrange(Cluster,p.adjust) #%>% dplyr::group_by(Cluster) %>% dplyr::top_n(-topN,p.adjust)  
p2 <- ggplot(KEGGresults,aes(x=Cluster,y=Description))+
  geom_point(aes(color = -log10(p.adjust+10^-27),size=Ratio))+
  labs(color = "-log10(Adjust P)",size="GeneRatio")+
  scale_colour_gradientn(limit=c(0,max(abs(-log10(KEGGresults$p.adjust+10^-27)))),
                         colors= rev(RColorBrewer::brewer.pal(11, "Spectral")[2:11]))+
  theme_bw()+
  scale_y_discrete(limit=rev(OurKEGG$Description[1:10]))
openxlsx::write.xlsx(p2$data,file=paste(OutPath,"/",folder,"/Fig4kdn.xlsx",sep=""),rownames = T)




