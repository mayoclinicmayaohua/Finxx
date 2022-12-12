
## STUDYADDRESS: /projects/bsi/fl/studies/s214070.Thompson.DSP2/re_analyze_FinXX950/data
## INVESTIGATOR:  Dr. Thompson
## STATISTICIAN:  Yaohua Ma



## aims: 1. find segment for each subset-data, 2. find mean by Finxx ID for each subdata,  3. using glm model find log2FC RFS yes vs no


## aims: 1. find segment for each subset-data, ------------------------------------------------------------------------------
library(ggplot2)
library(readxl)
library(writexl)
library(arsenal)

getwd()
setwd("/projects/bsi/fl/studies/s214070.Thompson.DSP2")

dat <- read_excel("re_analyze_FinXX950/data/r8_DSP_FinXXdata_normalized_by3housekeeping.xlsx")
 vars <- colnames(dat)
ct_v <- c("ID" ,       "Scan_Name" ,       "ROI_Name"  ,    "Tumor_Stroma.x" ,      "Segment_Name" ,  "Chemo_group" ,
          "RFS_status" ,    "alive_die"  ,   "classEnr"  ,     "AOI_surface_area",     "AOI_nuclei_count" ,    "FinXX_ID" )

dat[ct_v] <- lapply(dat[ct_v], factor)

cn_v <- vars[!(vars %in% ct_v)]  ## same as : cn_v <- vars[ vars %nin% ct_v]
cn_v29 <- c ("CD25"	, "TGFB1",	"CD20",	"PD_L1",	"VISTA",	"CD56",	"CTLA4",	"ICOS",
           "Beta_2_microglobulin",	"Tim_3",
           "CD127",	"B7_H3",	"PD_L2",	"CD40",	"IDO1",	"STING",	"CD11c",	"CD4",	"CD8",
           "GZMB",	"CD3",	"CD68",	"Ki_67",	"HLA_DR",	"PanCk",	"CD45",	"CD44",	"Fibronectin",	"SMA")


dat[cn_v] <- lapply(dat[cn_v], as.numeric)
# dat$CD86 <- dat$CD86 +1   ## avoid denominator=0 in  later model  logY~ x
# dat$CD40L <- dat$CD40L +1  ## avoid denominator=0 in  later model  logY~ x


protein_name <- cn_v
protein_name <- as.data.frame(protein_name)
# dat_T <- dat[which(dat$Chemo_group %in% c("TX_CEX")), ]  ###---------- find subset data ------------------
# dat_T_CEF <- dat[which(dat$Chemo_group %in% c("T_CEF")), ]

d_alltumor <- dat[which(dat$Tumor_Stroma.x  == c("Tumor")), ]
d_allstroma <- dat[which(dat$Tumor_Stroma.x  == c("Stroma")), ]

d_CD45tumor <- d_alltumor[which(d_alltumor$classEnr  == c("CD45enr")), ]
d_CD68tumor <- d_alltumor[which(d_alltumor$classEnr  == c("CD68Enr")), ]
d_tumortumor <- d_alltumor[which(d_alltumor$classEnr  == c("TumorEnr")), ]

d_CD45stroma <- d_allstroma[which(d_allstroma$classEnr  == c("CD45enr")), ]
d_CD68stroma <- d_allstroma[which(d_allstroma$classEnr  == c("CD68Enr")), ]

d_allsegment <- dat

## aims:  2. find mean by Finxx ID for each subdata, ------------------------------------------------

############################################################################
#below is find mean by clinic ID in tumor data, then  re-do it in stroma data
############################################################################
fun_mean <- function (dat) {
  library(arsenal)
  library(knitr)
  library(dplyr)
  dat <- as.data.frame(dat)
  vars <- colnames(dat)
  ct_v <- c("ID" ,       "Scan_Name" ,       "ROI_Name"  ,    "Tumor_Stroma.x" ,      "Segment_Name" ,  "Chemo_group" ,
            "RFS_status" ,    "alive_die"  ,   "classEnr"  ,     "AOI_surface_area",     "AOI_nuclei_count" ,    "FinXX_ID" )

  cn_v <- vars[!(vars %in% ct_v)]
  dat[ct_v] <- lapply(dat[ct_v], factor)
  dat[cn_v] <- lapply(dat[cn_v], as.numeric)

  n <-   cn_v[1:42]  ##  n can choose any columns in cn_v
  f <- as.formula (paste('FinXX_ID ~', paste(n[!n %in% 'FinXX_ID'], collapse = ' + ')))

  table_mean <- summary(tableby( f,
                                 data = dat,
                                 numeric.stats = c("mean"),
                                 numeric.test = "kwt",
                                 cat.test = "chisq",
                                 chisq.correct = F, digits = 2
  ), text = T, pfootnote = F)
  dat1 <- as.data.frame(table_mean)

  np1 <- dim(dat1)[1]
  np2 <- dim(dat1)[2] - 2
  even_indexes <- seq(2,np1,2)
  odd_indexes <- seq(1,np1,2)

  col1 <- dat1[,1]
  ProteinName <- data.frame(ProteinName = col1[odd_indexes])
  dat_exp <- dat1[even_indexes  , 2:np2]  ## col 313=total, col 314=pvalue , col 1= name ,   so take 2:312= 311 patients
  new <- cbind(ProteinName, dat_exp)
  clinic <-  colnames(new)[2:np2]
  clinic <- as.data.frame(clinic)

  newt <- t(new[, 2:np2])
  colnames(newt) <- ProteinName$ProteinName
  out <- cbind(clinic , newt)
  out <- as.data.frame(out)
   out1 <- out %>% separate(clinic, into=c("FinXX","clinic", "numb"), sep=" ") %>%
     arrange(clinic) %>% dplyr::select(-FinXX, -numb)
  return(out1)
}

table01 <- fun_mean(d_alltumor)
table02 <- fun_mean(d_allstroma)
table03 <- fun_mean(d_CD45tumor)
table04 <- fun_mean(d_CD68tumor)
table05 <- fun_mean(d_tumortumor)
table06 <- fun_mean(d_CD45stroma)
table07 <- fun_mean(d_CD68stroma)
table08 <- fun_mean(d_allsegment)



tableall00 <- list(table01,  table02, table03,table04, table05,table06, table07,table08)
##  write_xlsx (tableall00 , "/projects/bsi/fl/studies/s214070.Thompson.DSP2/re_analyze_FinXX950/output/data_FinXXDSP44_mean_byclinicID_r10.xlsx")

### -----------------data tidyverse done -----------------------------------good data done ---------- part 1 done ------------


 ####################################### match annotation and meanbycliniID, then find log2FC below ---------------------------
 d_44anno950 <- read_excel("re_analyze_FinXX950/data/r8_DSP_FinXXdata_normalized_by3housekeeping.xlsx")
 d_44anno950 <- d_44anno950[,1:12]
 d_44anno950 <- d_44anno950 %>% dplyr::select(FinXX_ID, everything())

 d_44annoclinic <- read_excel("/projects/bsi/fl/studies/s211340_thompson_DSP/DSPsurvivalRFS_Final44sample_45protein_950segment/data/r5_DSP_data_final_44_samples.xlsx")
 d_44annoclinic<- as.data.frame(d_44annoclinic)

 d_44anno <- d_44annoclinic[, c(1,33,34)]  ## can choose more columns :  tomor-stroma, chemoyesno, classenrich.....
  ## this data will be used 8 times below. re-analysis 950 data 44 samples

 d_alltumormean <- read_xlsx("re_analyze_FinXX950/output/data_FinXXDSP44_mean_byclinicID_r10.xlsx", sheet = "Sheet1")
 d_allstromamean <- read_xlsx("re_analyze_FinXX950/output/data_FinXXDSP44_mean_byclinicID_r10.xlsx", sheet = "Sheet2")

 d_CD45tumormean <- read_xlsx("re_analyze_FinXX950/output/data_FinXXDSP44_mean_byclinicID_r10.xlsx", sheet = "Sheet3")
 d_CD68tumormean <- read_xlsx("re_analyze_FinXX950/output/data_FinXXDSP44_mean_byclinicID_r10.xlsx", sheet = "Sheet4")
 d_tumortumormean <- read_xlsx("re_analyze_FinXX950/output/data_FinXXDSP44_mean_byclinicID_r10.xlsx", sheet = "Sheet5")

 d_CD45stromamean <- read_xlsx("re_analyze_FinXX950/output/data_FinXXDSP44_mean_byclinicID_r10.xlsx", sheet = "Sheet6")
 d_CD68stromamean <- read_xlsx("re_analyze_FinXX950/output/data_FinXXDSP44_mean_byclinicID_r10.xlsx", sheet = "Sheet7")

 d_allsegmentmean <- read_xlsx("re_analyze_FinXX950/output/data_FinXXDSP44_mean_byclinicID_r10.xlsx", sheet = "Sheet8")

 fun_merge <- function(dat1, dat2) {
   dat1_sort <-  dat1[order(dat1$PatientID),]
   dat2_sort <-  dat2 %>% rename(PatientID = clinic) %>% arrange(PatientID)
   d_RFS_merge <- merge(x=dat1_sort , y=dat2_sort, by="PatientID", all.y = TRUE)
   return( d_RFS_merge)
 }
    df1 <- fun_merge(d_44anno,  d_alltumormean) ## matched data = RFS+ meanby ID in tumor  patients
    df2 <- fun_merge(d_44anno,  d_allstromamean) ## matched data = RFS+ meanby ID
    df3 <- fun_merge(d_44anno,  d_CD45tumormean) ## matched data = RFS+ meanby ID
    df4 <- fun_merge(d_44anno,  d_CD68tumormean) ## matched data = RFS+ meanby ID
    df5 <- fun_merge(d_44anno,   d_tumortumormean) ## matched data = RFS+ meanby ID
    df6 <- fun_merge(d_44anno,  d_CD45stromamean) ## matched data = RFS+ meanby ID
    df7 <- fun_merge(d_44anno,   d_CD68stromamean) ## matched data = RFS+ meanby ID
    df8 <- fun_merge(d_44anno,   d_allsegmentmean) ## matched data = RFS+ meanby ID




## aims: 3. find fold change by RFSyes vs. no in tumor data , using generalized linear model   forest plot



library(arsenal)
library(tidyverse)
library(MASS)
require(foreign)

glm_function3 <- function( dat) {
   dat <- dat %>% mutate (RFS_status = Censor_RFS)    ### RFS_censor= RFS_status:  1 = recurrence free yes=RFS_yes ,  0=recurrence free- no
   dat <- dat %>% dplyr::select(-Censor_RFS) ## rfs 0=recurrence free yes = censor ,  rfs  1=recur free - no = death = event
  vars <- colnames(dat)
  ct_v <-  c("PatientID" ,  "RFS_status"   )

  # dat[1:6,60:66]    ## rfs_time = survival days,  RFS = survival years
  cn_v <- vars[!(vars %in% ct_v)]
  dat[ct_v] <- lapply(dat[ct_v], factor)
  dat[cn_v] <- lapply(dat[cn_v], as.numeric)
  dat$RFS_status <- factor(dat$RFS_status, levels =0:1, labels = c("RFS_no", "RFS_yes"))
  ## rfs : 0 = yes = free , 1=no = not free ,      ###   RFS_status 1=yes free  ,  0=no free
  dat <- as.data.frame(dat)
  dat <- dat[!is.na(dat$RFS_status), ]

  protein_name <- cn_v[2:43]
  protein_name <- as.data.frame(protein_name)
  out<-NULL
  for (i in 1:42){
    dependent_var <- cn_v[i+1]
    independent_var <- "RFS_status"
    f <- as.formula(paste0(dependent_var, "~", independent_var))
    m1 <- glm.nb(f, data = dat)
    ci1<-cbind(LogFC= coef(m1)[2], LogFCLower = confint(m1)[2,1], LogFCUpper = confint(m1)[2,2])
    log2FC <- ci1/log(2)  # or log2FC <-  log(ci2, 2)  = log2(ci2)
    ci2<-cbind(EstimateFC = exp(coef(m1)[2]), Lower = exp(confint(m1)[2,1]), Upper = exp(confint(m1)[2,2]))

    test <- summary(m1)
    pvalue <- test$coefficients[2, 4]
    ci <- cbind(ci2, log2FC, pvalue)
    out <- rbind(out, ci)
  }
  out<-as.data.frame(out)
  out<-data.frame(protein_name, out)

  ci95FC <-paste0(round(out$EstimateFC,digits=2), "(", round(out$Lower, digits=2),
                  ",", round(out$Upper, digits=2), ")")
  ci95logFC <-paste0(round(out$LogFC,digits=2), "(", round(out$LogFCLower, digits=2),
                     ",", round(out$LogFCUpper, digits=2), ")")
  p_value3digits <-  round(out$pvalue, digits=3)
  table01 <- data.frame(out, ci95FC, ci95logFC, p_value3digits)
  table01 <- table01 %>% dplyr::select(-pvalue)
  colnames(table01) <- c( "protein_name", "FC" , "FC_Lower", "FC_Upper", "log2FC",
                          "log2FC_Lower", "log2FC_Upper", "FC_95CI",  "log2FC_95CI", "p_value")
  return(table01)
}


table01 <- glm_function3(df1)
out001 <- as.data.frame(table01)
out001_29 <- out001 %>% filter(protein_name %in% cn_v29)  ## ". so get final 29 proteins"

table01 <- glm_function3(df2)
out002 <- as.data.frame(table01)
out002_29 <- out002 %>% filter(protein_name %in% cn_v29)  ## ". so get final 29 proteins"

table01 <- glm_function3(df3)
out003 <- as.data.frame(table01)
out003_29 <- out003 %>% filter(protein_name %in% cn_v29)  ## ". so get final 29 proteins"

table01 <- glm_function3(df4)
out004 <- as.data.frame(table01)
out004_29 <- out004 %>% filter(protein_name %in% cn_v29)  ## ". so get final 29 proteins"

table01 <- glm_function3(df5)
out005 <- as.data.frame(table01)
out005_29 <- out005 %>% filter(protein_name %in% cn_v29)  ## ". so get final 29 proteins"

table01 <- glm_function3(df6)
out006 <- as.data.frame(table01)
out006_29 <- out006 %>% filter(protein_name %in% cn_v29)  ## ". so get final 29 proteins"

table01 <- glm_function3(df7)
out007 <- as.data.frame(table01)
out007_29 <- out007 %>% filter(protein_name %in% cn_v29)  ## ". so get final 29 proteins"

table01 <- glm_function3(df8)
out008 <- as.data.frame(table01)
out008_29 <- out008 %>% filter(protein_name %in% cn_v29)  ## ". so get final 29 proteins"

readme <- data.frame(Sheet_Name=c( 'Sheet1' ,'Sheet2','Sheet3','Sheet4','Sheet5', 'Sheet6' ,'Sheet7','Sheet8','Sheet9',
                                   'Sheet10','Sheet11' ,'Sheet12','Sheet13','Sheet14','Sheet15', 'Sheet16'),
                     segmentName =c('tumor','stroma', 'CD45tumor', 'CD68tumor', 'tumortumor',
                                    'CD45stroma', 'CD68stroma', 'allsegments', 'tumor_29protein','stroma_29protein',
                                    'CD45tumor_29protein', 'CD68tumor_29protein', 'tumortumor_29protein',
                                    'CD45stroma_29protein', 'CD68stroma_29protein', 'allsegments_29protein' ),
                     comparison = c(rep("RFS_yes vs RFS_no", 16)))

out_log2FC_all <-  list( out001, out002, out003, out004, out005, out006, out007, out008 ,
                           out001_29, out002_29, out003_29 ,  out004_29, out005_29 ,  out006_29 ,  out007_29, out008_29 , readme )
write_xlsx (out_log2FC_all, "/projects/bsi/fl/studies/s214070.Thompson.DSP2/re_analyze_FinXX950/output/report2_FinXX_DSP950_log2FC_r10RFSyesvsno.xlsx")
########################################### good data done ----------------ready for forest plot------------------------------------

