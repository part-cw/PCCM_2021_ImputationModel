##Age (mo) <1, 1–11, 12–23, 24–59, 60–143, ≥ 144
VPS$age_qPELOD2_cat <- NULL
VPS$age_qPELOD2_cat[VPS$age.in.months>=0 & VPS$age.in.months<1] <- "0 days to <1 mo"
VPS$age_qPELOD2_cat[VPS$age.in.months>=1 & VPS$age.in.months<12] <- "1 mo to <1 yr"
VPS$age_qPELOD2_cat[VPS$age.in.months>=12 & VPS$age.in.months<2*12] <- "1 yr to <2 yrs"
VPS$age_qPELOD2_cat[VPS$age.in.months>=2*12 & VPS$age.in.months<5*12] <- "2 yrs to <5 yrs"
VPS$age_qPELOD2_cat[VPS$age.in.months>=5*12 & VPS$age.in.months<12*12] <- "5 yrs to <12 yrs"
VPS$age_qPELOD2_cat[VPS$age.in.months>=12*12 & VPS$age.in.months<18*12] <- "12 yrs to <18 yrs"
VPS$age_qPELOD2_cat <- factor(VPS$age_qPELOD2_cat, levels=c("0 days to <1 mo","1 mo to <1 yr","1 yr to <2 yrs","2 yrs to <5 yrs",
                                                            "5 yrs to <12 yrs","12 yrs to <18 yrs"),ordered=TRUE)
## HR: Table 1 Cutoff 207, 215, 203, 191, 176, 167
Tachycardia <- 
  (VPS$high.heart.rate..bpm.<207 & VPS$age_qPELOD2_cat=="0 days to <1 mo")|
  (VPS$high.heart.rate..bpm.<215 & VPS$age_qPELOD2_cat=="1 mo to <1 yr")| 
  (VPS$high.heart.rate..bpm.<203 & VPS$age_qPELOD2_cat=="1 yr to <2 yrs")| 
  (VPS$high.heart.rate..bpm.<191 & VPS$age_qPELOD2_cat=="2 yrs to <5 yrs")|
  (VPS$high.heart.rate..bpm.<176 & VPS$age_qPELOD2_cat=="5 yrs to <12 yrs")|
  (VPS$high.heart.rate..bpm.<167 & VPS$age_qPELOD2_cat=="12 yrs to <18 yrs")
Tachycardia[is.na(VPS$high.heart.rate..bpm.)] <- NA
VPS$qPELOD2_abnormal_HR <- Tachycardia

# TABLE 2. Systolic arterial pressure Cutoff 5: 63 79  87  90  95  98
#                                     Cutoff 1: 22, 38, 47, 49, 54, 58
SysBloodPress_Cutoff5 <- 
    (VPS$low.systolic.blood.pressure..mmhg.<63 & VPS$age_qPELOD2_cat=="0 days to <1 mo")|
    (VPS$low.systolic.blood.pressure..mmhg.<79 & VPS$age_qPELOD2_cat=="1 mo to <1 yr")| 
    (VPS$low.systolic.blood.pressure..mmhg.<87 & VPS$age_qPELOD2_cat=="1 yr to <2 yrs")| 
    (VPS$low.systolic.blood.pressure..mmhg.<90 & VPS$age_qPELOD2_cat=="2 yrs to <5 yrs")|
    (VPS$low.systolic.blood.pressure..mmhg.<95 & VPS$age_qPELOD2_cat=="5 yrs to <12 yrs")|
    (VPS$low.systolic.blood.pressure..mmhg.<98 & VPS$age_qPELOD2_cat=="12 yrs to <18 yrs")
SysBloodPress_Cutoff5[is.na(VPS$low.systolic.blood.pressure..mmhg.)]<-NA
SysBloodPress_Cutoff1 <- 
  (VPS$low.systolic.blood.pressure..mmhg.<22 & VPS$age_qPELOD2_cat=="0 days to <1 mo")|
  (VPS$low.systolic.blood.pressure..mmhg.<38 & VPS$age_qPELOD2_cat=="1 mo to <1 yr")| 
  (VPS$low.systolic.blood.pressure..mmhg.<47 & VPS$age_qPELOD2_cat=="1 yr to <2 yrs")| 
  (VPS$low.systolic.blood.pressure..mmhg.<49 & VPS$age_qPELOD2_cat=="2 yrs to <5 yrs")|
  (VPS$low.systolic.blood.pressure..mmhg.<54 & VPS$age_qPELOD2_cat=="5 yrs to <12 yrs")|
  (VPS$low.systolic.blood.pressure..mmhg.<58 & VPS$age_qPELOD2_cat=="12 yrs to <18 yrs")
SysBloodPress_Cutoff1[is.na(VPS$low.systolic.blood.pressure..mmhg.)]<-NA
VPS$qPELOD2_abnormal_SBP <- SysBloodPress_Cutoff5

# TABLE 3: Glasgow Coma Score; 3-4, 5-10
VPS$qPELOD2_abnormal_Mentation <- VPS$worst.glascow.coma.score <11

# Cleanup
rm(Tachycardia,SysBloodPress_Cutoff5,SysBloodPress_Cutoff1)
