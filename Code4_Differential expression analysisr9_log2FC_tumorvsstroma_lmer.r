
## STUDYADDRESS: /projects/bsi/fl/studies/s214070.Thompson.DSP2/re_analyze_FinXX950/data
## INVESTIGATOR:  Dr. Thompson
## STATISTICIAN:  Yaohua Ma

## aims : 1. find each subset data, 2, using linear mixed model adjusted by Finxx iD find log2FC
##########################################################   part 1   ##########################################
##------------------------------- compare tumor vs stroma---------CD45tumor vs CD68tumor ------- CD45tumor vs tumortumor ---------
##-------------------------------CD68tumor vs tumortumor--------CD45stroma vs CD68stroma ------- -------

library(ggplot2)
library(readr)
library(readxl)
library(writexl)
library(tidyverse)
library(dplyr)
getwd()
setwd("/projects/bsi/fl/studies/s214070.Thompson.DSP2")

dat  <- read_excel("re_analyze_FinXX950/data/r8_DSP_FinXXdata_normalized_by3housekeeping.xlsx")
dat <- dat %>% mutate (class = paste0(classEnr, Tumor_Stroma.x))
vars <- colnames(dat)
ct_v <- c("ID" ,                  "Scan_Name"  ,          "ROI_Name"  ,           "Tumor_Stroma.x" ,
          "Segment_Name"  ,       "Chemo_group",          "RFS_status"  ,         "alive_die"  ,
          "classEnr" ,            "AOI_surface_area" ,    "AOI_nuclei_count" ,    "FinXX_ID" ,  "class"   )
cn_v <- vars[!(vars %in% ct_v)]
dat[ct_v] <- lapply(dat[ct_v], factor)
dat[cn_v] <- lapply(dat[cn_v], as.numeric)

dat$CD86 <- dat$CD86 +1  ## avoid denominator=0 in lmer model  logY~ x
dat$CD40L <- dat$CD40L +1  ## avoid denominator=0 in  lmer  model  logY~ x
protein_name <- cn_v[1:42]
protein_name <- as.data.frame(protein_name)

d_all <- dat  ##    $Tumor_Stroma.x    each involve 2 group : tumor , stroma ,  stroma is baseline
d_CD45 <- dat[which(dat$classEnr == c("CD45enr")), ] ## each involve 2 group : tumor , stroma  ##  data for later
d_CD68 <- dat[which(dat$classEnr == c("CD68Enr")), ]  ##  data for later
d_Tumor <- dat[which(dat$classEnr == c("TumorEnr")), ]  ##  data for later

d_CD45Tumor_CD68Tumor <- dat[which(dat$class == c("CD45enrTumor", "CD68EnrTumor" )), ] ## each involve 2 class
levels(d_CD45Tumor_CD68Tumor$class) <- list( "CD68EnrTumor" = c("CD68EnrTumor"), "CD45enrTumor" = c("CD45enrTumor") )

d_CD45Tumor_TumorTumor <- dat[which(dat$class == c("CD45enrTumor",  "TumorEnrTumor")), ]
levels(d_CD45Tumor_TumorTumor$class) <- list( "TumorEnrTumor" = c("TumorEnrTumor"), "CD45enrTumor" = c("CD45enrTumor") )

d_CD68Tumor_TumorTumor <- dat[which(dat$class == c( "CD68EnrTumor", "TumorEnrTumor")), ]
levels(d_CD68Tumor_TumorTumor$class) <- list( "TumorEnrTumor" = c("TumorEnrTumor"), "CD68EnrTumor" = c("CD68EnrTumor") )

d_CD45Stroma_CD68Stroma <- dat[which(dat$class == c("CD45enrStroma", "CD68EnrStroma" )), ]
levels(d_CD45Stroma_CD68Stroma$class) <- list( "CD68EnrStroma" = c("CD68EnrStroma"), "CD45enrStroma" = c("CD45enrStroma") )

d_CD45Stroma_TumorStroma <- dat[which(dat$class == c("CD45enrStroma",  "TumorEnrStroma")), ]  ##  data for later
d_CD68Stroma_TumorStroma <- dat[which(dat$class == c( "CD68EnrStroma", "TumorEnrStroma")), ]  ##  data for later





    ################################ find log2FC by tumor vs stroma     #####################################################
library(arsenal)
library(MASS)
require(foreign)
library(lme4)
library(lmerTest)
options(lmerControl=list(check.nobs.vs.rankZ = "warning",
                         check.nobs.vs.nlev = "warning",
                         check.nobs.vs.nRE = "warning",
                         check.nlev.gtreq.5 = "warning",
                         check.nlev.gtr.1 = "warning"))

fun_lmer1 <- function( dat) {
  dat <- as.data.frame(dat)  ## data can not include NA
  dat <- dat[!is.na(dat$Tumor_Stroma.x), ]
  vars <- colnames(dat)
  ct_v <- c("ID" ,                  "Scan_Name"  ,          "ROI_Name"  ,           "Tumor_Stroma.x" ,
            "Segment_Name"  ,       "Chemo_group",          "RFS_status"  ,         "alive_die"  ,
            "classEnr" ,            "AOI_surface_area" ,    "AOI_nuclei_count" ,    "FinXX_ID" ,  "class"   )
  cn_v <- vars[!(vars %in% ct_v)]
  dat[ct_v] <- lapply(dat[ct_v], factor)
  dat[cn_v] <- lapply(dat[cn_v], as.numeric)

  protein_name <- cn_v[1:42]
  protein_name <- as.data.frame(protein_name)

  out<-NULL
  for (i in 1:42){
    dependent_var <- cn_v[i]
    independent_var <- "Tumor_Stroma.x"
    f <- as.formula(paste0("log2(",  dependent_var, ")", "~", independent_var,  "+ (1|FinXX_ID)"))
    mixed.lmer <- lmer(f, data = dat)    ### use   ----------------------- data ---------- carefully---------------------
    m1 <- summary(mixed.lmer)
    log2FC <- round(m1$coefficients[2,1], 3)
    ci95 <- confint(mixed.lmer)
    log2ci <- round(ci95[4, ], digits = 3)
    log2FC_Lower <- log2ci[1]
    log2FC_Upper <- log2ci[2]

    FC <- round(2^(log2FC), 3)
    FC_Lower <- round (2^(log2FC_Lower), digits = 3)
    FC_Upper <- round (2^(log2FC_Upper), digits = 3)
    # t_value <-  m1$coefficients[2,4]
    # n <- length(m1$residuals)
    # p_value <- 2*pt(-abs(t_value),df=n-1)
    # p_value <- round(p_value, digits = 3)
    p_value  <-  m1$coefficients[2,5]
    p_value <- round(p_value, digits = 3)
    FC_95CI <- paste0(FC,   "(", FC_Lower,", ", FC_Upper, ")")
    log2FC_95CI <- paste0(log2FC,   "(", log2FC_Lower, ", ", log2FC_Upper, ")")
    ci <- data.frame(  FC, FC_Lower, FC_Upper, log2FC, log2FC_Lower, log2FC_Upper, FC_95CI,  log2FC_95CI, p_value)
    out <- rbind(out, ci)
  }
  out<-as.data.frame(out)
  out<-data.frame(protein_name, out)
  colnames(out)<-c( "protein_name", "FC" , "FC_Lower", "FC_Upper", "log2FC",
                    "log2FC_Lower", "log2FC_Upper", "FC_95CI",  "log2FC_95CI", "p_value")
 return(out)
}

  table01 <- fun_lmer1(d_all)  ## log2FC tumor vs stroma

  out001 <- as.data.frame(table01)

c_29 <- c (  "CD25"	, "TGFB1",	"CD20",	"PD_L1",	"VISTA",	"CD56",	"CTLA4",	"ICOS",
             "Beta_2_microglobulin",	"Tim_3", "CD127",	"B7_H3",	"PD_L2",	"CD40",	"IDO1",	   "STING",	"CD11c",	"CD4",
             "CD8",    "GZMB",	 "CD3",	"CD68",	"Ki_67",	"HLA_DR",	"PanCk",	"CD45",	"CD44",	"Fibronectin",	"SMA")
out001_29 <- out001 %>% filter(protein_name %in% c_29)  ## "TGFB1 not in 756 DSP data . so get final 28 proteins"


fun_lmer2 <- function( dat) {
  dat <- as.data.frame(dat)
  dat <- dat[!is.na(dat$class), ]
  vars <- colnames(dat)
  ct_v <- c("ID" ,                  "Scan_Name"  ,          "ROI_Name"  ,           "Tumor_Stroma.x" ,
            "Segment_Name"  ,       "Chemo_group",          "RFS_status"  ,         "alive_die"  ,
            "classEnr" ,            "AOI_surface_area" ,    "AOI_nuclei_count" ,    "FinXX_ID" ,  "class"   )
  cn_v <- vars[!(vars %in% ct_v)]
  dat[ct_v] <- lapply(dat[ct_v], factor)
  dat[cn_v] <- lapply(dat[cn_v], as.numeric)

  protein_name <- cn_v[1:42]
  protein_name <- as.data.frame(protein_name)
  out<-NULL
  for (i in 1:42){
    dependent_var <- cn_v[i]
    independent_var <- "class"
    f <- as.formula(paste0("log2(",  dependent_var, ")", "~", independent_var,  "+ (1|FinXX_ID)"))
    mixed.lmer <- lmer(f, data = dat)
    m1 <- summary(mixed.lmer)
    log2FC <- round(m1$coefficients[2,1], 3)
    ci95 <- confint(mixed.lmer)
    log2ci <- round(ci95[4, ], digits = 3)
    log2FC_Lower <- log2ci[1]
    log2FC_Upper <- log2ci[2]

    FC <- round(2^(log2FC), 3)
    FC_Lower <- round (2^(log2FC_Lower), digits = 3)
    FC_Upper <- round (2^(log2FC_Upper), digits = 3)
    # t_value <-  m1$coefficients[2,4]
    # n <- length(m1$residuals)
    # p_value <- 2*pt(-abs(t_value),df=n-1)
    # p_value <- round(p_value, digits = 3)
    p_value  <-  m1$coefficients[2,5]
    p_value <- round(p_value, digits = 3)
    FC_95CI <- paste0(FC,   "(", FC_Lower,", ", FC_Upper, ")")
    log2FC_95CI <- paste0(log2FC,   "(", log2FC_Lower, ", ", log2FC_Upper, ")")
    ci <- data.frame(  FC, FC_Lower, FC_Upper, log2FC, log2FC_Lower, log2FC_Upper, FC_95CI,  log2FC_95CI, p_value)
    out <- rbind(out, ci)
  }
  out<-as.data.frame(out)
  out<-data.frame(protein_name, out)
  colnames(out)<-c( "protein_name", "FC" , "FC_Lower", "FC_Upper", "log2FC",
                    "log2FC_Lower", "log2FC_Upper", "FC_95CI",  "log2FC_95CI", "p_value")
  return(out)
}


table02 <- fun_lmer2(d_CD45Tumor_CD68Tumor)
out002 <- as.data.frame(table02)
out002_29 <- out002 %>% filter(protein_name %in% c_29)  ## "TGFB1 not in 756 DSP data . so get final 28 proteins"


table0 <- fun_lmer2(d_CD45Tumor_TumorTumor)
out003 <- as.data.frame(table0)
out003_29 <- out003 %>% filter(protein_name %in% c_29)  ##  so get final 29 proteins"

table0 <- fun_lmer2(d_CD68Tumor_TumorTumor)
out004 <- as.data.frame(table0)
out004_29 <- out004 %>% filter(protein_name %in% c_29)  ##  so get final 29 proteins"

table0 <- fun_lmer2(d_CD45Stroma_CD68Stroma )
out005 <- as.data.frame(table0)
out005_29 <- out005 %>% filter(protein_name %in% c_29)  ## ". so get final 29proteins"

##################################################################################################################part 4
## aims: 4. repeat 1,2,3, use 29 protein only  just select 29 protein from part 1,2,3,output,  done for part 4  result
readme <- data.frame(Sheet_Name=c( 'Sheet1' ,'Sheet2','Sheet3','Sheet4','Sheet5', 'Sheet6' ,'Sheet7','Sheet8','Sheet9','Sheet10'),
                     comparison =c('tumor_vs_stroma', 'CD45tumor_vs_CD68tumor', 'CD45tumor_vs_tumortumor',
                                   'CD68tumor_vs_tumortumor', 'CD45stroma_vs_CD68stroma',
                                   'tumor_vs_stroma29protein', 'CD45tumor_vs_CD68tumor29protein', 'CD45tumor_vs_tumortumor29protein',
                                   'CD68tumor_vs_tumortumor29protein', 'CD45stroma_vs_CD68stroma29protein'))

out_log2FC_all <-  list( out001, out002, out003, out004, out005,  out001_29, out002_29, out003_29 ,  out004_29, out005_29 , readme )
## write_xlsx (out_log2FC_all, "/projects/bsi/fl/studies/s214070.Thompson.DSP2/re_analyze_FinXX950/output/report1_FinXX_DSP950_log2FC_r9tumorvsstroma.xlsx")
########################################### good data done ----------------ready for forest plot------------------------------------





table(d_all$classEnr)
library(arsenal)
table000000 <-  summary(tableby( Tumor_Stroma.x ~ classEnr + Chemo_group +	RFS_status,
                 data = d_all,
                 numeric.stats = c("mean"),
                 numeric.test = "kwt",
                 cat.test = "chisq",
                 chisq.correct = F, digits = 2
), text = T, pfootnote = F)

out_freq <- kable(table000000)  ## kable is good for write2word, data.frame is good for write_xlsx
write2word (out_freq, "/projects/bsi/fl/studies/s214070.Thompson.DSP2/re_analyze_FinXX950/output/report1_FinXX_DSP950_freq.docx")
########################################### good data done ----------------ready for basic check ----------------------------------


# all(is.na(dummy))
# which(is.na(dummy))
# any(is.na(d_CD45Tumor_CD68Tumor))


