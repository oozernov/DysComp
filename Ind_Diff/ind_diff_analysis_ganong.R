Packages <- c("dplyr", "stats", "psych", "ggplot2", "lme4", "lmerTest","Jmisc",
              "gridExtra","olsrr",'relaimpo','BayesFactor','MASS','psych','mice','VIM')
lapply(Packages, library, character.only = TRUE)
setwd("~/Dropbox (MIT)/GitHub/DysComp2/Ind_Diff")


###load and organize ABCD ###
abcd=read.csv("~/Dropbox (MIT)/Annals_SVR/ABCD.csv")
abcd_raw=read.csv('ABCD_raw.csv')
abcd=merge(abcd,abcd_raw)
#Clean data
names(abcd)[names(abcd)=="wrmt_lc_ss_2"] <- "LC"
names(abcd)[names(abcd)=="gort_comp_ss_2"] <- "RC"
names(abcd)[names(abcd)=="gort_ori_ss_2"] <- "ORI"
names(abcd)[names(abcd)=="wrmt_wa_ss_2"] <- "WA"
names(abcd)[names(abcd)=="wrmt_id_ss_2"] <- "WID"
names(abcd)[names(abcd)=="towre_sw_ss_2"] <- "SWE"
names(abcd)[names(abcd)=="towre_pde_ss_2"] <- "PDE"
names(abcd)[names(abcd)=="ran_letters_ss_2"] <- "RANL"
names(abcd)[names(abcd)=="ran_objects_ss_2"] <- "RANO"
names(abcd)[names(abcd)=="kbit_ss_2"] <- "IQ"
names(abcd)[names(abcd)=="ppvt_vocab_ss_2"] <- "Vocab"
names(abcd)[names(abcd)=="ctopp_elision_ss_2"] <- "Elision"
names(abcd)[names(abcd)=="ctopp_blending_ss_2"] <- "Blending"
names(abcd)[names(abcd)=="ctopp_nonword_ss_2"] <- "NonWord"
names(abcd)[names(abcd)=="wais_dsf_ss_2"] <- "DigitsFw"
names(abcd)[names(abcd)=="wais_dsb_ss_2"] <- "DigitsBck"

abcd$Sex<-as.factor(abcd$Sex)

## calculate dyslexia 
abcd$DD<-ifelse(abcd$WA < 90 & abcd$WID < 90 | abcd$WA < 90 & abcd$SWE < 90|
                  abcd$WA < 90 & abcd$PDE < 90 |
                  abcd$WID < 90 & abcd$SWE < 90 | abcd$WID <90 & abcd$PDE < 90 |
                  abcd$SWE < 90 & abcd$PDE < 90,1,
                ifelse (abcd$WA >= 90 & abcd$WID >= 90 & abcd$PDE >= 90 & abcd$SWE >= 90,0,NA))
#abcd$DD2=ifelse(abcd$subgroup==1 & abcd$DD==0,0,ifelse(abcd$subgroup==2 & abcd$DD==1,1,"NA"))

#convert to factors
abcd$Sex<-as.factor(abcd$Sex)
abcd$Group<-as.factor(abcd$Group)
abcd$DD<-as.factor(abcd$DD)

abcd2<-abcd%>%dplyr::select('ID','Age','IQ','Sex','DD','WID','WA','Elision','Blending','NonWord','PDE','SWE','Vocab',
                            'RANL','RANO','DigitsFw','DigitsBck','LC','RC','ORI')

#impute missing scores
set.seed(182)
init = mice(abcd2, maxit=0)
meth = init$method
predM = init$predictorMatrix
predM[, c("ID")]=0
imputed = mice(abcd2, method=meth, predictorMatrix=predM, m=50)
imputed <- complete(imputed)
#apply(imputed,2,pMiss)
abcd3<-imputed
abcd3$Study<-'abcd'
abcd3$a_c<-'a'
###load and organize READ ###
read=read.csv("~/Dropbox (MIT)/Annals_SVR/READ.csv",na.strings=c("NA","NaN", "9999","8888","999","7777"))

read$DD=ifelse(read$W3WAss_T4 < 90 & read$W3WIss_T4 < 90 | read$W3WAss_T4 < 90 & read$TSWEss_T4 < 90|
                 read$W3WAss_T4 < 90 & read$TPDEss_T4 < 90 |
                 read$W3WIss_T4 < 90 & read$TSWEss_T4 < 90 | read$W3WIss_T4 <90 & read$TPDEss_T4 < 90 |
                 read$TSWEss_T4 < 90 & read$TPDEss_T4 < 90,"1", 
               ifelse (read$W3WAss_T4 >= 90 & read$W3WIss_T4 >= 90 & read$W3WAss_T4 >= 90 & read$TSWEss_T4 >= 90,"0","NA"))

read$Sex=ifelse(read$Sex=="Female",1,2)


names(read)[names(read)=="W3WAss_T4"] <- "WA"
names(read)[names(read)=="W3WIss_T4"] <- "WID"
names(read)[names(read)=="RANLss_T4"] <- "RANL"
names(read)[names(read)=="RANOss_T4"] <- "RANO"
names(read)[names(read)=="CELF5ss_T4"] <- "LC"
names(read)[names(read)=="GORTACOMPss_T4"] <- "RC"
names(read)[names(read)=="GORTORIss_T4"] <- "ORI"
names(read)[names(read)=="TSWEss_T4"] <- "SWE"
names(read)[names(read)=="TPDEss_T4"] <- "PDE"
names(read)[names(read)=="KBITss_T4"] <- "IQ"
names(read)[names(read)=="PPVTss_T4"] <- "Vocab"
names(read)[names(read)=="CTELss_T4"] <- "Elision"
names(read)[names(read)=="CTBWss_T4"] <- "Blending"
names(read)[names(read)=="CTNRss_T4"] <- "NonWord"
names(read)[names(read)=="CTMDss_T4"] <- "MD"
names(read)[names(read)=="WISCDFss_T4"] <- "DigitsFw"
names(read)[names(read)=="WISCDBss_T4"] <- "DigitsBck"

read$T4mos<-as.integer(read$T4mos)
read$Age<-read$T4mos/12

read2<-read%>%dplyr::select('ID','IQ','Sex','Age','DD','WID','WA','Elision','Blending','NonWord','PDE','SWE','Vocab',
                            'RANL','RANO','DigitsFw','DigitsBck','LC','RC','ORI')


#impute missing scores
set.seed(182)
init = mice(read2, maxit=0)
meth = init$method
predM = init$predictorMatrix
predM[, c("ID")]=0
imputed = mice(read2, method=meth, predictorMatrix=predM, m=50)
imputed <- complete(imputed)
#apply(imputed,2,pMiss)
read3<-imputed
read3$Study<-'read'
read3$a_c<-'c'

### READER
reader=read.csv("~/Dropbox (MIT)/Annals_SVR/READER.csv",na.strings=c("invalid"))
reader <- reader %>% filter(ID != "READER_124")
reader <- reader %>%dplyr::filter(reader$ID!='READER_inel')

reader$DD=ifelse(reader$wrmt_wa_std_behav1 < 90 & reader$wrmt_id_std_behav1 < 90 | reader$wrmt_wa_std_behav1 < 90 & reader$towre_sight_std_behav1 < 90|
                   reader$wrmt_wa_std_behav1 < 90 & reader$towre_phon_std_behav1 < 90 |
                   reader$wrmt_id_std_behav1 < 90 & reader$towre_sight_std_behav1 < 90 | reader$wrmt_id_std_behav1 <90 & reader$towre_phon_std_behav1 < 90 |
                   reader$towre_sight_std_behav1 < 90 & reader$towre_phon_std_behav1 < 90,1,
                 ifelse (reader$wrmt_wa_std_behav1 >= 90 & reader$wrmt_id_std_behav1 >= 90 & reader$towre_phon_std_behav1 >= 90 & reader$towre_sight_std_behav1 >= 90,0,NA))

names(reader)[names(reader)=="celf_usp_std_behav1"] <- "LC"
names(reader)[names(reader)=="gort_comp_std_behav1"] <- "RC"
names(reader)[names(reader)=="gort_ori_std_behav1"] <- "ORI"
names(reader)[names(reader)=="wrmt_wa_std_behav1"] <- "WA"
names(reader)[names(reader)=="wrmt_id_std_behav1"] <- "WID"
names(reader)[names(reader)=="ran_let_std_behav1"] <- "RANL"
names(reader)[names(reader)=="ran_obj_std_behav1"] <- "RANO"
names(reader)[names(reader)=="wisc_matrix_std_behav1"]<-"IQ"
names(reader)[names(reader)=="reader$calculated_age_behav_1"]<-"Age"
names(reader)[names(reader)=="towre_sight_std_behav1"] <- "SWE"
names(reader)[names(reader)=="towre_phon_std_behav1"] <- "PDE"
names(reader)[names(reader)=="ef_vocab_std_behav1"] <- "Vocab"
names(reader)[names(reader)=="ctopp_elis_std_behav1"] <- "Elision"
names(reader)[names(reader)=="ctopp_bw_std_behav1"] <- "Blending"
names(reader)[names(reader)=="ctopp_nwr_std_behav1"] <- "NonWord"
names(reader)[names(reader)=="wisc_digitfwd_std_behav1"] <- "DigitsFw"
names(reader)[names(reader)=="wisc_digitbkwd_std_behav1"] <- "DigitsBck"
reader$DD<-as.factor(reader$DD)
reader$Sex<-as.factor(reader$Sex)

reader2<-reader%>%dplyr::select('ID','IQ','Sex','Age','DD','WID','WA','Elision','Blending','NonWord','PDE','SWE','Vocab',
                            'RANL','RANO','DigitsFw','DigitsBck','LC','RC','ORI')


Missing_patterns <-md.pattern(reader2)
Missing_patterns
#impute missing scores
set.seed(182)
init = mice(reader2, maxit=0)
meth = init$method
predM = init$predictorMatrix
predM[, c("ID")]=0
imputed = mice(reader2, method=meth, predictorMatrix=predM, m=50)
imputed <- complete(imputed)
#apply(imputed,2,pMiss)
reader3<-imputed
reader3$Study<-'reader'
reader3$a_c<-'c'

d_all<-rbind(read3, reader3)
d_all<-rbind(d_all,abcd3)
write.csv(d_all,"all_beh.csv")
