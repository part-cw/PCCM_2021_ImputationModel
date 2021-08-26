VPS$age_PALS_cat <- NULL
VPS$age_PALS_cat[VPS$age.in.months>=0 & VPS$age.in.months<1] <- "0 days to <1 mo"
VPS$age_PALS_cat[VPS$age.in.months>=1 & VPS$age.in.months<3] <- "1 mo to <3 mo"
VPS$age_PALS_cat[VPS$age.in.months>=3 & VPS$age.in.months<12] <- "3 mo to <1 yr"
VPS$age_PALS_cat[VPS$age.in.months>=12 & VPS$age.in.months<2*12] <- "1 yr to <2 yrs"
VPS$age_PALS_cat[VPS$age.in.months>=2*12 & VPS$age.in.months<10*12] <- "2 yrs to <10 yrs"
VPS$age_PALS_cat[VPS$age.in.months>=10*12 & VPS$age.in.months<18*12] <- "10 yrs to <18 yrs"
VPS$age_PALS_cat <- factor(VPS$age_PALS_cat, levels=c("0 days to <1 mo","1 mo to <3 mo","3 mo to <1 yr","1 yr to <2 yrs",
                                                      "2 yrs to <10 yrs","10 yrs to <18 yrs"),ordered=TRUE)
## HR:  
Tachycardia <- (VPS$high.heart.rate..bpm.>205 & VPS$age_PALS_cat=="0 days to <1 mo")|
               (VPS$high.heart.rate..bpm.>205 & VPS$age_PALS_cat=="1 mo to <3 mo")|
               (VPS$high.heart.rate..bpm.>190 & VPS$age_PALS_cat=="3 mo to <1 yr")|
               (VPS$high.heart.rate..bpm.>190 & VPS$age_PALS_cat=="1 yr to <2 yrs")|
               (VPS$high.heart.rate..bpm.>140 & VPS$age_PALS_cat=="2 yrs to <10 yrs")|
               (VPS$high.heart.rate..bpm.>100 & VPS$age_PALS_cat=="10 yrs to <18 yrs")
Tachycardia[is.na(VPS$high.heart.rate..bpm.)] <- NA
VPS$PALS_abnormal_HR <- Tachycardia

# SBP
SysBloodPress <- 
    (VPS$low.systolic.blood.pressure..mmhg.<60 & VPS$age_PALS_cat=="0 days to <1 mo")|
    (VPS$low.systolic.blood.pressure..mmhg.<70 & VPS$age_PALS_cat=="1 mo to <3 mo")| 
    (VPS$low.systolic.blood.pressure..mmhg.<70 & VPS$age_PALS_cat=="3 mo to <1 yr")| 
    (VPS$low.systolic.blood.pressure..mmhg.<(70+floor(VPS$age.in.months)/12*2) & VPS$age_PALS_cat=="1 yr to <2 yrs")|
    (VPS$low.systolic.blood.pressure..mmhg.<(70+floor(VPS$age.in.months)/12*2) & VPS$age_PALS_cat=="2 yrs to <10 yrs")|
    (VPS$low.systolic.blood.pressure..mmhg.<90 & VPS$age_PALS_cat=="10 yrs to <18 yrs")
SysBloodPress[is.na(VPS$low.systolic.blood.pressure..mmhg.)]<-NA
VPS$PALS_abnormal_SBP <- SysBloodPress

# RR [Infant= <1year, Toddler= 1-3yrs, Preschooler 4-5, School Age = 6-12, Adolescent = 12-18]
RespRate <- 
  (VPS$high.respiratory.rate..bpm.>60 & VPS$age.in.months<12)|
  (VPS$high.respiratory.rate..bpm.>40 & VPS$age.in.months>=1*12 & VPS$age.in.months<4*12) | 
  (VPS$high.respiratory.rate..bpm.>34 & VPS$age.in.months>=4*12 & VPS$age.in.months<6*12) |
  (VPS$high.respiratory.rate..bpm.>30 & VPS$age.in.months>=6*12 & VPS$age.in.months<12*12) |
  (VPS$high.respiratory.rate..bpm.>16 & VPS$age.in.months>=12*12 & VPS$age.in.months<18*12)
RespRate[is.na(VPS$high.respiratory.rate..bpm.)]<-NA
VPS$PALS_abnormal_RR <- RespRate

# GCS as per Sepsis-3
#VPS$PALS_abnormal_Mentation <- VPS$worst.glascow.coma.score <11
VPS$PALS_abnormal_Mentation <- VPS$worst.glascow.coma.score <=13

# Cleanup
rm(Tachycardia,SysBloodPress,RespRate)
