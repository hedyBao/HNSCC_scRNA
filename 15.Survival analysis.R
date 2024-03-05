### survival 
rm(list=ls())
library(magrittr)
library(maxstat)
library(ggplot2)
library(dplyr)
library(DOSE)
library(clusterProfiler) #GSEA
library(fgsea)
library(magrittr)
library(tibble)
library(RColorBrewer)
options(stringsAsFactors = F)
# Clinical 
ClinicalData <- read.delim("/data/public/TCGA_data/TCGA_preliminary/TCGA_ClinicalData_20180420.txt")
ClinicalData$PFI.time <- as.numeric(as.character(ClinicalData$PFI.time))
ClinicalData$OS.time <- as.numeric(as.character(ClinicalData$OS.time))

OutPath <- "/work/brj/Collaboration/scRNA/2022/HNSCC/Revise/Fig/Survival/"

#signature score
GSVAscore <- readr::read_rds("/work/brj/Collaboration/scRNA/2022/HNSCC/Revise/Result/HNSC_Regulons.rds.gz")
names(GSVAscore)[2] <- "Score"
GSVAscore <- GSVAscore %>% dplyr::mutate(
  xcell = purrr::map(.x=Score,function(.x){ 
    try(
      {	
        #.x=Expr$expr[[1]] 
        result=data.frame(geneset=rownames(.x),.x)#%>%
      },silent = F)	
  })) 

colnames(GSVAscore)[1]='cancer_types'
Candidate <- "RegulonsTFDP1" #RegulonsMYBL2
#1.提取 目标基因
GSVAscore %>% dplyr::mutate(Expr=purrr::map(.x=xcell,function(.x){
  .x %>% as.data.frame() %>% dplyr::filter(geneset %in% Candidate) %>% tidyr::gather(key="barcode",value="Expression",-geneset)})) %>%
  dplyr::mutate(Expr = purrr::map(.x=Expr,function(.x){.x %>% 
      dplyr::mutate(types = data.frame(do.call(rbind,strsplit(barcode,"-")))$X1,
                    sample_class = data.frame(do.call(rbind,strsplit(barcode,"-")))$X2)})) -> Candidateexp # 3colnames(Gene,barcode())

#.x <- Candidateexp$Expr[[1]] %>% as.data.frame()
#spls <- substr(.x$barcode,14,16)


#2.目标基因分高低组 
Candidateexp <- Candidateexp %>% dplyr::mutate(Candidate_sur=purrr::map(.x=Expr,function(.x){
  .x <- .x %>% dplyr::filter(substr(barcode,14,16) %in% c("01A")) %>%
    dplyr::mutate(barcode = gsub("\\.","-",substr(barcode,1,12))) %>% 
    dplyr::inner_join(ClinicalData,by=c("barcode"="bcr_patient_barcode"))
  .x <- .x %>% dplyr::mutate(OS.time=as.numeric(OS.time)/30,PFI.time=as.numeric(PFI.time)/30) %>% dplyr::filter(!is.na(OS.time))
  .x <- .x %>% 
    dplyr::mutate(ExpressionGroup = ifelse(Expression > median(Expression),"High","Low"),
                  OS = as.numeric(as.character(OS)),
                  PFI = as.numeric(as.character(PFI)))
  aa_try <- try(aa <- maxstat.test(Surv(OS.time, OS) ~ Expression, data= .x,smethod="LogRank"),silent=TRUE)
  .x <- .x %>% dplyr::mutate(ExpressionGroupB = ifelse(Expression >aa_try$estimate,"High","Low"))
  return(.x)
}))

Candidateexp <- Candidateexp %>% dplyr::mutate(Candidate_HR=purrr::map(.x=Candidate_sur,function(.x){
  model1 <- coxph(survival::Surv(OS.time, OS) ~ Expression, data=.x,  na.action=na.exclude)
  HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=4)
  Coxp <- signif(as.numeric(summary(model1)$coefficients[1,c( "Pr(>|z|)")]),digit=4)
  HR_detail = summary(model1)
  CILow =  HR_detail$conf.int[,"lower .95"]
  CIHigh =  HR_detail$conf.int[,"upper .95"]
  model1 <- survival::survdiff(survival::Surv(OS.time, OS) ~ ExpressionGroup,data=.x, na.action=na.exclude)
  
  KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(.x$ExpressionGroup)))-1)
  
  model2 <- survival::survdiff(survival::Surv(OS.time, OS) ~ ExpressionGroupB,data=.x, na.action=na.exclude)
  
  KMP2 <- 1-pchisq(model2$chisq, df=length(levels(factor(.x$ExpressionGroupB)))-1)
  
  print(nrow(.x))
  data.frame(cancertypes=unique(.x$type),HR, Coxp,CILow,CIHigh,KMP,BestKMP=KMP2)
  data.frame(cancertypes=unique(.x$type),HR, Coxp,CILow,CIHigh,KMP)
}))

Candidateexp <- Candidateexp %>% dplyr::mutate(Candidate_plot=purrr::map(.x= Candidate_sur,function(.x){
  model1 <- survival::coxph(survival::Surv(OS.time, OS) ~ Expression, data=.x,  na.action=na.exclude)
  HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=4)
  Coxp <- signif(as.numeric(summary(model1)$coefficients[1,c( "Pr(>|z|)")]),digit=4)
  HR_detail = summary(model1)
  CI =  paste(signif(HR_detail$conf.int[,"lower .95"],2),'-',signif(HR_detail$conf.int[,"upper .95"],2),sep=" ") #Result[i1, c('CI_95%_for_HR')]
  fit<- survival::survfit(survival::Surv(OS.time, OS) ~ ExpressionGroup,data=.x)
  model1 <- survival::survdiff(survival::Surv(OS.time, OS) ~ ExpressionGroup,data=.x, na.action=na.exclude)
  
  KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(.x$ExpressionGroup)))-1)
  
  p <- survminer::ggsurvplot(fit,palette = c("red", "blue"),
                             # risk.table = TRUE, risk.table.col = "strata",
                             pval = F, title=paste(unique(.x$type),": p = ",signif(as.numeric(KMP),digits = 2),
                                                   "\nHR ",HR, " [95% CI:",CI,"]",sep=""),
                             fun = function(y) y*100,legend = c(0.2, 0.2), legend.title=F,font.legend=10,font.title=10,
                             legend.labs = c(paste("High (n=",table(.x$ExpressionGroup)[[1]],")",sep=""),
                                             paste("Low (n=",table(.x$ExpressionGroup)[[2]],")",sep="")))
  p$plot <- p$plot+ labs(x="Months") + theme(legend.title = element_blank(),legend.background = element_blank())
  return(p$plot)
}))

Candidateexp <- Candidateexp %>% dplyr::mutate(Candidate_plotBestCurv=purrr::map(.x=Candidate_sur,function(.x){
  model1 <- survival::coxph(survival::Surv(OS.time, OS) ~ Expression, data=.x,  na.action=na.exclude)
  HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=4)
  Coxp <- signif(as.numeric(summary(model1)$coefficients[1,c( "Pr(>|z|)")]),digit=4)
  HR_detail = summary(model1)
  CI =  paste(signif(HR_detail$conf.int[,"lower .95"],2),'-',signif(HR_detail$conf.int[,"upper .95"],2),sep=" ") #Result[i1, c('CI_95%_for_HR')]
  fit<- survival::survfit(survival::Surv(OS.time, OS) ~ ExpressionGroupB,data=.x)
  model1 <- survival::survdiff(survival::Surv(OS.time, OS) ~ ExpressionGroupB,data=.x, na.action=na.exclude)
  
  KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(.x$ExpressionGroupB)))-1)
  
  p <- survminer::ggsurvplot(fit,palette = c("red", "blue"),
                             # risk.table = TRUE, risk.table.col = "strata",
                             pval = F, title=paste(unique(.x$type),": p = ",signif(as.numeric(KMP),digits = 2),
                                                   "\nHR ",HR, " [95% CI:",CI,"]",sep=""),
                             fun = function(y) y*100,legend = c(0.2, 0.2), legend.title=F,font.legend=10,font.title=10,
                             legend.labs = c(paste("High (n=",table(.x$ExpressionGroupB)[[1]],")",sep=""),
                                             paste("Low (n=",table(.x$ExpressionGroupB)[[2]],")",sep="")))
  p$plot <- p$plot + labs(x="Months") + theme(legend.title = element_blank(),legend.background = element_blank())
  return(p$plot)
}))


pdf(paste(OutPath,Candidate,"_OS_Survival.pdf",sep=""),width = 5,height = 5)
pp <- cowplot::plot_grid(plotlist = Candidateexp$Candidate_plot,labels = Candidate)
print(pp)
dev.off()

pdf(paste(OutPath,Candidate,"BC_OS_Survival.pdf",sep=""),width =  5,height = 5)
pp <- cowplot::plot_grid(plotlist = Candidateexp$Candidate_plotBestCurv,labels= Candidate)
print(pp)
dev.off()















