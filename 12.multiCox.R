### 加上HPV 及margin 重新做multiCox 
rm(list=ls())
library(tidyverse)
library(tibble)
library(data.table)
library(survival)
library(survminer)
### 1. 5 tumor subcluster : Tumor sub 使用的分组是median
TumorSub_proportion <- read.delim('/server2_work/smy/Project/HNSCC_26sample/2.data/CIBERSORTx_Results_TumorSubcluster.txt')
TumorSub_proportion$Mixture <- gsub('[.]','-',TumorSub_proportion$Mixture)
TumorSub_proportion <- TumorSub_proportion[str_sub(TumorSub_proportion$Mixture,14,16) %in% c("01A"),]
TumorSub_proportion$Mixture <- str_sub(TumorSub_proportion$Mixture,1,12)
TumorSub <- TumorSub_proportion %>% dplyr::select(c(Mixture,cluster_0,cluster_1,cluster_2,cluster_3,cluster_4))

### PSOTN fib, RSPO 
GSVAscore <- readr::read_rds("/work/brj/Collaboration/scRNA/2022/HNSCC/Revise/Result/HNSC_POSTN_SPP1.rds.gz")
GSVAF <- GSVAscore$POSTN_SPP1_score[[1]] %>% t() %>% as.data.frame() 
GSVAF$Mixture <- substr(rownames(GSVAF),0,12)
HNSC_SignatureF <- GSVAF[substr(rownames(GSVAF),14,16) %in% c("01A"),]
HNSC_SignatureF$Mixture <- gsub("[.]","-",HNSC_SignatureF$Mixture)

HNSC_POSTN_SPP1 <- HNSC_SignatureF %>% dplyr::select(c(Mixture,POSTNTop50,SPP1Top50))


GSVAscore1 <- readr::read_rds("/work/brj/Collaboration/scRNA/2022/HNSCC/Revise/Result/HNSC_RSPO1_FOLR2.rds.gz")
GSVAF1 <- GSVAscore1$RSPO1_FOLR2_score[[1]] %>% t() %>% as.data.frame() 
GSVAF1$Mixture <- substr(rownames(GSVAF1),0,12)
HNSC_SignatureF1 <- GSVAF1[substr(rownames(GSVAF1),14,16) %in% c("01A"),]
HNSC_SignatureF1$Mixture <- gsub("[.]","-",HNSC_SignatureF1$Mixture)

HNSC_RSPO1_FOLR2 <- HNSC_SignatureF1 %>% dplyr::select(c(Mixture,RSPO1Top50,FOLR2Top50))

# SignatureD
GSVAscore <- readr::read_rds("/work/brj/Collaboration/scRNA/2022/HNSCC/Result/HNSC_ExprFPKM_SigantureDE.rds.gz")
GSVAD <- GSVAscore$CD8A_IFNg_TcellSingaling_score[[1]] %>% 
  t() %>% as.data.frame() %>%
  dplyr::filter(substr(rownames(.),14,16) %in% "01A") %>%
  dplyr::mutate(Mixture = substr(rownames(.),0,12)) %>%
  dplyr::select(c(Mixture,SignatureD))
GSVAD$Mixture <- gsub('[.]','-',GSVAD$Mixture)


HNSC_FPKM <- readr::read_rds("/server2_work/smy/Public/data/RNA-seq/TCGA/HNSCC_FPKM.rds")
HNSC_FPKM_F <- HNSC_FPKM[c("MYBL2","TFDP1"),]
HNSC_FPKM_F <- t(HNSC_FPKM_F) %>% as.data.frame()
HNSC_FPKM_F[,1:2] <- apply(HNSC_FPKM_F[,1:2],2,function(x)log2(x+1))
HNSC_FPKM_F$Mixture <- substr(rownames(HNSC_FPKM_F),0,12)
HNSC_FPKM_FF <- HNSC_FPKM_F[substr(rownames(HNSC_FPKM_F),14,16) %in% c("01A"),]
HNSC_FPKM_FF$Mixture <- str_sub(HNSC_FPKM_FF$Mixture,1,12)

HNSCC_proportion <- TumorSub %>% 
  inner_join(GSVAD) %>% 
  inner_join(HNSC_FPKM_FF) %>%
  inner_join(HNSC_POSTN_SPP1) %>%
  inner_join(HNSC_RSPO1_FOLR2) 

##clinical info old
survivaO <- read.delim("/work/brj/GEO/CRC/multiCox/TCGA_ClinicalData_20180420.txt")#,","gender"
survivaOld <- survivaO %>% dplyr::select("bcr_patient_barcode",'OS','OS.time',"age_at_initial_pathologic_diagnosis","gender","ajcc_pathologic_tumor_stage") %>% 
  dplyr::filter(bcr_patient_barcode %in% HNSCC_proportion$Mixture) %>% 
  dplyr::mutate(OS.time = as.numeric(as.numeric(as.character(.$OS.time))/30),
                age_at_initial_pathologic_diagnosis = as.numeric(.$age_at_initial_pathologic_diagnosis)) %>% 
  dplyr::filter(!is.na(OS.time))  %>% 
  set_colnames(c("Mixture","OS","OS.time","age","gender","stage"))

##clinical info new 
survivaN <- openxlsx::read.xlsx("/work/brj/Collaboration/scRNA/2022/HNSCC/Data/TCGAHNSCclinical.xlsx")
survival <- survivaN %>% dplyr::select("sampleID","margin_status") %>% #,"hpv_status_by_p16_testing"
  set_colnames(c("sampleID","margin")) %>% #,"HPV"
  dplyr::mutate(Mixture = substr(.$sampleID,0,12)) %>%
  dplyr::inner_join(survivaOld)%>%
  #dplyr::filter(!is.na(HPV)) %>%
  dplyr::filter(!is.na(margin)) %>%
  dplyr::select(-sampleID) %>%
  dplyr::select(c(Mixture,OS,OS.time),everything()) %>%
  dplyr::filter(!duplicated(.))

names(survival) <- c("sample","status","time","margin","age","gender","stage") # stage 有的是432 ,,"stage" ,"HPV"

survival$stage[grep("Discrepancy",survival$stage)] <- NA
survival$stage[grep("Not",survival$stage)] <- NA
survival$stage[grep("Unknown",survival$stage)] <- NA

### 1.按照Stage I、II、III、IV 的条件区分
survival$stage[grep("Stage I$|IA",survival$stage)]  <- "Stage I"
survival$stage[grep("Stage II$|IIA|IIB|IIC",survival$stage)]  <- "Stage II"
survival$stage[grep("III",survival$stage)]  <- "Stage III"
survival$stage[grep("IV",survival$stage)]  <- "Stage IV"
survival$stage <- factor(survival$stage,levels = c("Stage I","Stage II","Stage III","Stage IV"))

### 2.按照Stage I、II、III、IV 的条件区分
survival$stage[grep("Stage I$|IA|II$|IIA|IIB|IIC",survival$stage)]  <- "Early_Stage"
survival$stage[grep("III|IV",survival$stage)]  <- "Advanced_Stage"
survival$stage <- factor(survival$stage,levels = c("Early_Stage","Advanced_Stage"))

survival <- survival %>%
  dplyr::filter(!is.na(age)) %>% 
  dplyr::filter(!is.na(gender)) %>%
  dplyr::filter(!is.na(stage)) 


mulCox <- function(types){
  print(types)
  
  HNSC_proportion <- HNSCC_proportion %>% dplyr::select(c(Mixture,types)) %>% set_colnames(c("sample","Cluster"))
  surv_expr <- survival %>%
    inner_join(HNSC_proportion) %>% 
    column_to_rownames("sample") %>%
    dplyr::filter(!is.na(stage)) 
  which(is.na(surv_expr))
  surv_expr$status <- as.numeric(as.character(surv_expr$status))
  surv_expr$time <- as.numeric(as.character(surv_expr$time)) 
  surv_expr$age <- as.numeric(surv_expr$age)
  surv_expr$gender <- as.factor(surv_expr$gender)
  surv_expr$Group <- ifelse(surv_expr[,"Cluster"] > median(surv_expr[,"Cluster"]),paste0(types,'_High'),paste0(types,'_Low'))
  surv_expr$Group <- factor(surv_expr$Group,levels = c(paste0(types,'_Low'),paste0(types,'_High')))
  surv_expr$status <- as.numeric(as.character(surv_expr$status)) 
  surv_expr$time <- as.numeric(as.character(surv_expr$time)) 
  #surv_expr$HPV <- as.factor(surv_expr$HPV)
  surv_expr$margin <- as.factor(surv_expr$margin)
  res.cox <- coxph(Surv(time, status) ~ Group+age+gender+margin+stage, data = surv_expr,na.action = na.exclude)
  
  summary(res.cox)
  tty <- ggforest(
    model = res.cox,
    data = surv_expr,
    main = "Hazard ratio",
    #cpositions = c(0.00, 0.20, 0.35),
    fontsize = 0.6,
    noDigits = 2)
  return(tty)
}

types <- c(names(HNSCC_proportion)[-1])

plist <- lapply(types, mulCox)

OutPath <- "/work/brj/Collaboration/scRNA/2022/HNSCC/Revise2/Fig/mCox/"
pdf(paste(OutPath,"mCoxHLStage1-2_TumorSigantureNoHPV.pdf",sep=""),width = 24,height = 12)
pp <- cowplot::plot_grid(plotlist = plist,ncol = 4)
print(pp)
dev.off()


