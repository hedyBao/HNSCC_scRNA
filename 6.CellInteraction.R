# Cell-cell interaction
cellchatMyeFib <- readRDS('/work/brj/Collaboration/2022/scRNA/HNSCC/Data/26sample_definetype_MyeFib_cellchat_stageSplit.rds')
cellchatMyeFib <- cellchatMyeFib[c(6,7,2,1,4,3,5,8)]
NetWeightMyeFib <- lapply(c(1:8), function(sub){
  mat <- cellchatMyeFib[[sub]]@net$weight
  mat <- as.data.frame(mat)
  MyeFib <- grep("macro|macor|DC|Mast|myeloid|Mono|fibroblast|Fib",colnames(mat),value = T)
  mat2 <- mat[c("POSTN+ fibroblast","RSPO1+ fibroblast","SPP1+ macorphages","FOLR2+ macorphages"),MyeFib] 
  mat2$stage <- names(cellchatMyeFib)[sub]
  return(mat2)
})

NetWeightMyeFib <- bind_rows(NetWeightMyeFib)
NetWeightMyeFib_melt <- melt(NetWeightMyeFib)
NetWeightMyeFib_melt$cluster <- rep(c("POSTN+ fibroblast","RSPO1+ fibroblast","SPP1+ macrophages","FOLR2+ macrophages"),times = 152)# 152: length(NetWeightMyeFib_melt$variable)/4
NetWeightMyeFib_melt$cellchat <- paste(NetWeightMyeFib_melt$cluster,'to',NetWeightMyeFib_melt$variable,sep = ' ')

NetWeightMyeFib_meltF2 <- NetWeightMyeFib_melt %>% 
  filter(cluster %in% c('POSTN+ fibroblast')) %>% 
  filter(variable %in% c("SPP1+ macorphages","FOLR2+ macorphages")) %>%
  filter(stage %in% c("NT","Pre","E","A")) %>%
  mutate(cellchat2 = paste(cellchat,stage))
NetWeightMyeFib_meltF2$stage <- factor(NetWeightMyeFib_meltF2$stage,levels = c("NT","Pre","E","A"))
NetWeightMyeFib_meltF2$cellchat2 <- factor(NetWeightMyeFib_meltF2$cellchat2,levels = rev(unique(NetWeightMyeFib_meltF2$cellchat2)))
Yorder <- c('POSTN+ fibroblast to SPP1+ macorphages NT','POSTN+ fibroblast to SPP1+ macorphages Pre',
            "POSTN+ fibroblast to SPP1+ macorphages E","POSTN+ fibroblast to SPP1+ macorphages A",
            '',
            "POSTN+ fibroblast to FOLR2+ macorphages NT",'POSTN+ fibroblast to FOLR2+ macorphages Pre',
            "POSTN+ fibroblast to FOLR2+ macorphages E","POSTN+ fibroblast to FOLR2+ macorphages A")
NetWeightMyeFib_meltF1 <- NetWeightMyeFib_melt %>% 
  filter(cluster %in% c('SPP1+ macrophages')) %>% 
  filter(variable %in% c("POSTN+ fibroblast","RSPO1+ fibroblast")) %>%
  filter(stage %in% c("NT","Pre","E","A")) %>%
  mutate(cellchat2 = paste(cellchat,stage))



NetWeightMyeFib_meltF1$stage <- factor(NetWeightMyeFib_meltF1$stage,levels = c("NT","Pre","E","A"))
NetWeightMyeFib_meltF1$cellchat2 <- factor(NetWeightMyeFib_meltF1$cellchat2,levels = rev(unique(NetWeightMyeFib_meltF1$cellchat2)))
Yorder <- c('SPP1+ macrophages to POSTN+ fibroblast NT','SPP1+ macrophages to POSTN+ fibroblast Pre',
            "SPP1+ macrophages to POSTN+ fibroblast E","SPP1+ macrophages to POSTN+ fibroblast A",
            '',
            "SPP1+ macrophages to RSPO1+ fibroblast NT",'SPP1+ macrophages to RSPO1+ fibroblast Pre',
            "SPP1+ macrophages to RSPO1+ fibroblast E","SPP1+ macrophages to RSPO1+ fibroblast A")
HNSCC_Whole <- readRDS('/work/smy/Project/HNSCC/2.data/DeinfeTypes/HNSCC_26sampleAll_DefineTypes.rds.gz')
Idents(HNSCC_Whole) <- HNSCC_Whole$DefineTypes
HNSCC_Whole$celltype_RE <- HNSCC_Whole$DefineTypes

## add tumor in epi : metadata$Malignancy  == 'malignant'
metadata <- readRDS('/work/smy/Project/HNSCC_26sample/2.data/26sample_tumor_cell_metadata.rds')
malignant <- rownames(metadata[grep('malignant',metadata$Malignancy),])
HNSCC_Whole$celltype_RE[malignant] <- 'TumorCell'

HNSCC_TumorFib_Allstage <- subset(HNSCC_Whole,subset = celltype_RE %in% c("POSTN+ fibroblast","TumorCell"))
HNSCC_TumorFib <- subset(HNSCC_TumorFib_Allstage,subset = stage %in% c("NT","Pre","E","A"))
HNSCC_TumorFib$celltype_Re <- paste(HNSCC_TumorFib$celltype_RE,HNSCC_TumorFib$stage,sep = "_")
#HNSCC_TumorFib$celltype_aggregate <- ifelse(HNSCC_TumorFib$celltype_Re %in% c("TumorCell_E","TumorCell_NT","TumorCell_Pre"),"TumorCell_NTPreE",HNSCC_TumorFib$celltype_Re)

HNSCC_TumorFib$celltype_aggregate <- HNSCC_TumorFib$celltype_Re

seurat_obj <- HNSCC_TumorFib
celltype_id = "celltype_aggregate" # metadata column name of the cell type of interest
seurat_obj = SetIdent(seurat_obj, value = seurat_obj[[celltype_id]])

####Read in the NicheNet ligand-receptor network and ligand-target matrix
# ##the following of human origin
ligand_target_matrix <- readRDS("/work/brj/Collaboration/2022/scRNA/HNSCC/NichenetrData/ligand_target_matrix.rds")
lr_network <- readRDS("/work/brj/Collaboration/2022/scRNA/HNSCC/NichenetrData/lr_network.rds")
lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)

organism = "human" # user adaptation required on own dataset



####nichnet
#####1. Define the niches/microenvironments of interest
#! Important: your receiver cell type should consist of 1 cluster!
table(seurat_obj$celltype_aggregate)

niches = list(
  "Stage_A" = list(
    "sender" = "POSTN+ fibroblast_A",
    "receiver" = "TumorCell_A"),
  "Stage_E" = list(
    "sender" = c("POSTN+ fibroblast_E"),
    "receiver" = "TumorCell_E"),
  "Stage_Pre" = list(
    "sender" = c("POSTN+ fibroblast_Pre"),
    "receiver" = "TumorCell_Pre"),
  "Stage_NT" = list(
    "sender" = c("POSTN+ fibroblast_NT"),
    "receiver" = "TumorCell_NT")
) # user adaptation required on own dataset



#####2. Calculate differential expression between the niches # 得到差异性受体配体矩阵，受体--配体
assay_oi = "RNA" # other possibilities: RNA,...
DE_sender = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), niches = niches, type = "sender", assay_oi = assay_oi) # only ligands important for sender cell types
DE_receiver = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), niches = niches, type = "receiver", assay_oi = assay_oi) # only receptors now, later on: DE analysis to find targets
DE_sender = DE_sender %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
DE_receiver = DE_receiver %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
expression_pct = 0.10  ###Process DE results and filter
DE_sender_processed = process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")
DE_receiver_processed = process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")
specificity_score_LR_pairs = "min_lfc" ###Combine sender-receiver DE based on L-R pairs:
DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)


#####3. Optional: Calculate differential expression between the different spatial regions
include_spatial_info_sender = FALSE # if not spatial info to include: put this to false # user adaptation required on own dataset
include_spatial_info_receiver = FALSE # if spatial info to include: put this to true # user adaptation required on own dataset
#spatial_info = tibble(celltype_region_oi = "CAF_High", celltype_other_region = "myofibroblast_High", niche =  "pEMT_High_niche", celltype_type = "sender") # user adaptation required on own dataset
#specificity_score_spatial = "lfc"
# this is how this should be defined if you don't have spatial info
# mock spatial info
if(include_spatial_info_sender == FALSE & include_spatial_info_receiver == FALSE){
  spatial_info = tibble(celltype_region_oi = NA, celltype_other_region = NA) %>% mutate(niche =  niches %>% names() %>% head(1), celltype_type = "sender")
} 

if(include_spatial_info_sender == TRUE){
  sender_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "sender"))
  sender_spatial_DE_processed = process_spatial_de(DE_table = sender_spatial_DE, type = "sender", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
  
  # add a neutral spatial score for sender celltypes in which the spatial is not known / not of importance
  sender_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)
  
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
  
} else {
  # # add a neutral spatial score for all sender celltypes (for none of them, spatial is relevant in this case)
  sender_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))  
  
}
## [1] "Calculate Spatial DE between: CAF_High and myofibroblast_High"

if(include_spatial_info_receiver == TRUE){
  receiver_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "receiver"))
  receiver_spatial_DE_processed = process_spatial_de(DE_table = receiver_spatial_DE, type = "receiver", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
  
  # add a neutral spatial score for receiver celltypes in which the spatial is not known / not of importance
  receiver_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)
  
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  
} else {
  # # add a neutral spatial score for all receiver celltypes (for none of them, spatial is relevant in this case)
  receiver_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
}

#####4. Calculate ligand activities and infer active ligand-target links 配体-- target
lfc_cutoff = 0.15 # recommended for 10x as min_lfc cutoff. 
specificity_score_targets = "min_lfc"
DE_receiver_targets = calculate_niche_de_targets(seurat_obj = seurat_obj, niches = niches, lfc_cutoff = lfc_cutoff, expression_pct = expression_pct, assay_oi = assay_oi) 
DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches, expression_pct = expression_pct, specificity_score = specificity_score_targets)

background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
geneset_niche1 = DE_receiver_processed_targets %>% filter(receiver == niches[[1]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
geneset_niche2 = DE_receiver_processed_targets %>% filter(receiver == niches[[2]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
geneset_niche3 = DE_receiver_processed_targets %>% filter(receiver == niches[[3]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
geneset_niche4 = DE_receiver_processed_targets %>% filter(receiver == niches[[4]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()


# Good idea to check which genes will be left out of the ligand activity analysis (=when not present in the rownames of the ligand-target matrix).
# If many genes are left out, this might point to some issue in the gene naming (eg gene aliases and old gene symbols, bad human-mouse mapping)
geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))
geneset_niche3 %>% setdiff(rownames(ligand_target_matrix))
geneset_niche4 %>% setdiff(rownames(ligand_target_matrix))

print(length(geneset_niche1))##We recommend having between 20 and 1000 genes in the geneset of interest
print(length(geneset_niche2))
print(length(geneset_niche3))
print(length(geneset_niche4))

top_n_target = 250


niche_geneset_list = list(
  "Stage_A" = list(
    "receiver" = niches[[1]]$receiver,
    "geneset" = geneset_niche1,
    "background" = background),
  "Stage_E" = list(
    "receiver" = niches[[2]]$receiver,
    "geneset" = geneset_niche2 ,
    "background" = background),
  "Stage_Pre" = list(
    "receiver" = niches[[3]]$receiver,
    "geneset" = geneset_niche3 ,
    "background" = background),
  "Stage_NT" = list(
    "receiver" = niches[[4]]$receiver,
    "geneset" = geneset_niche4 ,
    "background" = background)
)
ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)
#[1] "Calculate Ligand activities for: TumorCell_A"
#[1] "Calculate Ligand activities for: SNCG+ Basal_R"

#####5. Calculate (scaled) expression of ligands, receptors and targets across cell types of interest (log expression values and expression fractions)
features_oi = union(lr_network$ligand, lr_network$receptor) %>% union(ligand_activities_targets$target) %>% setdiff(NA)
dotplot = suppressWarnings(Seurat::DotPlot(seurat_obj %>% subset(idents = niches %>% unlist() %>% unique()), features = features_oi, assay = assay_oi))
exprs_tbl = dotplot$data %>% as_tibble()
exprs_tbl = exprs_tbl %>% rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
            mutate(fraction = fraction/100) %>% as_tibble() %>% select(celltype, gene, expression, expression_scaled, fraction) %>% distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))
exprs_tbl_ligand = exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction) 
exprs_tbl_receptor = exprs_tbl %>% filter(gene %in% lr_network$receptor) %>% rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
exprs_tbl_target = exprs_tbl %>% filter(gene %in% ligand_activities_targets$target) %>% rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)

exprs_tbl_ligand = exprs_tbl_ligand %>%  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% mutate(ligand_fraction_adapted = ligand_fraction) %>% mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))
exprs_tbl_receptor = exprs_tbl_receptor %>% mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% mutate(receptor_fraction_adapted = receptor_fraction) %>% mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))
#####6. Expression fraction and receptor
exprs_sender_receiver = lr_network %>% 
  inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
  inner_join(exprs_tbl_receptor, by = c("receptor")) %>% inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))
ligand_scaled_receptor_expression_fraction_df = exprs_sender_receiver %>% group_by(ligand, receiver) %>% mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% mutate(ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% distinct() %>% ungroup() 

#####7. Prioritization of ligand-receptor and ligand-target links
prioritizing_weights = c("scaled_ligand_score" = 5,
                         "scaled_ligand_expression_scaled" = 1,
                         "ligand_fraction" = 1,
                         "scaled_ligand_score_spatial" = 2, 
                         "scaled_receptor_score" = 0.5,
                         "scaled_receptor_expression_scaled" = 0.5,
                         "receptor_fraction" = 1, 
                         "ligand_scaled_receptor_expression_fraction" = 1,
                         "scaled_receptor_score_spatial" = 0,
                         "scaled_activity" = 0,
                         "scaled_activity_normalized" = 1,
                         "bona_fide" = 1)
output = list(DE_sender_receiver = DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed,
              ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand,  exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target)
prioritization_tables = get_prioritization_tables(output, prioritizing_weights)
prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(receiver == niches[[2]]$receiver) %>% head(10)
prioritization_tables$prioritization_tbl_ligand_target %>% filter(receiver == niches[[2]]$receiver) %>% head(10)

####8. Visualization of the Differential NicheNet output
top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)
top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand, receptor) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)
ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, prioritization_score) %>% group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% ungroup() %>% distinct() %>% inner_join(top_ligand_niche_df) %>% filter(niche == top_niche) %>% group_by(niche) %>% top_n(20, prioritization_score) %>% ungroup() # get the top50 ligands per niche

outpath <- "/work/brj/Collaboration/2022/scRNA/HNSCC/Result/nichenetr/"
targetcell <- "POSTNTumorsplit"

receiver_oi = "TumorCell_A" 
sender_oi="POSTN+ fibroblast_A"
filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% filter(sender == sender_oi) %>% pull(ligand) %>% unique()
prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 
lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, prioritized_tbl_oi, prioritization_tables$prioritization_tbl_ligand_receptor, plot_legend = FALSE, heights = NULL, widths = NULL)

pdf(paste(outpath,targetcell,"/1.Top20ligand_plot.pdf",sep = ""),width = 10,height = 10)
lfc_plot
dev.off()


##Ligand expression, activity and target genes 配体表达，活性和靶基因
exprs_activity_target_plot = make_ligand_activity_target_exprs_plot(receiver_oi, prioritized_tbl_oi,  prioritization_tables$prioritization_tbl_ligand_receptor,  prioritization_tables$prioritization_tbl_ligand_target, output$exprs_tbl_ligand,  output$exprs_tbl_target, lfc_cutoff, ligand_target_matrix, plot_legend = T, heights = 1, widths = 1)

pdf(paste(outpath,targetcell,"/2.LigandExpressionActivityTargetgenes_plotLegend.pdf",sep = ""),width = 40,height = 10)
exprs_activity_target_plot$combined_plot
dev.off()

filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% filter(sender == sender_oi)%>% top_n(20, prioritization_score) %>% pull(ligand) %>% unique()
prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 
exprs_activity_target_plot = make_ligand_activity_target_exprs_plot(receiver_oi, prioritized_tbl_oi,  prioritization_tables$prioritization_tbl_ligand_receptor,  prioritization_tables$prioritization_tbl_ligand_target, output$exprs_tbl_ligand,  output$exprs_tbl_target, lfc_cutoff, ligand_target_matrix, plot_legend = FALSE, heights = NULL, widths = NULL)

pdf(paste(outpath,targetcell,"/3.LigandExpressionActivityTargetgenes_Filterplot.pdf",sep = ""),width = 40,height = 10)
exprs_activity_target_plot$combined_plot
dev.off()


####### L-R pairs 圈图# 展示四个阶段
# barplot(c(1:5),col = c("lavender","#F7F7F7","#CCCCCC", "#969696","#636363"))
receiver_oi <- c("TumorCell_NT","TumorCell_Pre","TumorCell_E","TumorCell_A")
sender_oi <-  c("POSTN+ fibroblast_NT","POSTN+ fibroblast_Pre","POSTN+ fibroblast_E","POSTN+ fibroblast_A")
filtered_ligands =  prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(receiver == receiver_oi) %>% filter(sender == sender_oi)%>% top_n(20, prioritization_score) %>% pull(ligand) %>% unique()
prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% dplyr::select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

prioritized_tbl_oi$sender <- factor(prioritized_tbl_oi$sender,levels = sender_oi)
prioritized_tbl_oi$receiver <- factor(prioritized_tbl_oi$receiver,levels = receiver_oi )


colors_sender = RColorBrewer::brewer.pal(n = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% length(), name = 'Set2') %>% magrittr::set_names(prioritized_tbl_oi$sender %>% unique() %>% sort())
colors_receiver = c("#E0E0E0","#BABABA","#878787","#4D4D4D")  %>% magrittr::set_names(prioritized_tbl_oi$receiver %>% unique() %>% sort())
circos_output = make_circos_lr(prioritized_tbl_oi, colors_sender, colors_receiver)

pdf(paste(outpath,targetcell,"/4.Circosplot.pdf",sep = ""),width = 7,height =6)
circos_output$p_circos
dev.off()
pdf(paste(outpath,targetcell,"/4.Circosplotlegend.pdf",sep = ""),width = 5,height = 5)
circos_output$p_legend
dev.off()
