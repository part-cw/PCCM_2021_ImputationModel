####################
#--- calculate SIRS criteria based on Goldstein 2005, PMID: 15636651
# NOTE: They screwed up their SYS thresholds, which got corrected in a letter to the editor: PMID 16003219
####################
VPS$age_cat <- NULL
VPS$age_cat[VPS$age.in.months>=0 & VPS$age.in.months<0.25] <- "0 days to <1 wk"
VPS$age_cat[VPS$age.in.months>=0.25 & VPS$age.in.months<1] <- "1 wk to <1 mo"
VPS$age_cat[VPS$age.in.months>=1 & VPS$age.in.months<2*12] <- "1 mo to <2 yrs"
VPS$age_cat[VPS$age.in.months>=2*12 & VPS$age.in.months<6*12] <- "2 yrs to <6 yrs" 
VPS$age_cat[VPS$age.in.months>=6*12 & VPS$age.in.months<13*12] <- "6 yrs to <13 yrs"
VPS$age_cat[VPS$age.in.months>=13*12 & VPS$age.in.months<18*12] <- "13 yrs to <18 yrs"
VPS$age_cat <- factor(VPS$age_cat, levels=c("0 days to <1 wk","1 wk to <1 mo","1 mo to <2 yrs","2 yrs to <6 yrs",
                                            "6 yrs to <13 yrs","13 yrs to <18 yrs"),ordered=TRUE)

## TEMP:
CoreTemp <- VPS$low.temperature..c.<36 | VPS$high.temperature..c. > 38.5
CoreTemp[is.na(VPS$low.temperature..c.)|is.na(VPS$high.temperature..c.)] <- NA
VPS$abnormal_Temp <- CoreTemp

## HR:  
Tachycardia <- 
  (VPS$high.heart.rate..bpm.>180 & VPS$age_cat=="0 days to <1 wk")|
  (VPS$high.heart.rate..bpm.>180 & VPS$age_cat=="1 wk to <1 mo") | 
  (VPS$high.heart.rate..bpm.>180 & VPS$age_cat=="1 mo to <2 yrs") |
  (VPS$high.heart.rate..bpm.>140 & VPS$age_cat=="2 yrs to <6 yrs") |
  (VPS$high.heart.rate..bpm.>130 & VPS$age_cat=="6 yrs to <13 yrs") |
  (VPS$high.heart.rate..bpm.>110 & VPS$age_cat=="13 yrs to <18 yrs")
Tachycardia[is.na(VPS$high.heart.rate..bpm.)] <- NA

Bradycardia <- 
  (VPS$low.heart.rate..bpm.<100 & VPS$age_cat=="0 days to <1 wk")|
  (VPS$low.heart.rate..bpm.<100 & VPS$age_cat=="1 wk to <1 mo") |
  (VPS$low.heart.rate..bpm.<90 & VPS$age_cat=="1 mo to <2 yrs")
Bradycardia[is.na(VPS$low.heart.rate..bpm.)] <- NA
VPS$abnormal_HR = Tachycardia|Bradycardia

# RR
RespRate <- 
  (VPS$high.respiratory.rate..bpm.>50 & VPS$age_cat=="0 days to <1 wk")|
  (VPS$high.respiratory.rate..bpm.>40 & VPS$age_cat=="1 wk to <1 mo") | 
  (VPS$high.respiratory.rate..bpm.>34 & VPS$age_cat=="1 mo to <2 yrs") |
  (VPS$high.respiratory.rate..bpm.>22 & VPS$age_cat=="2 yrs to <6 yrs") |
  (VPS$high.respiratory.rate..bpm.>18 & VPS$age_cat=="6 yrs to <13 yrs") |
  (VPS$high.respiratory.rate..bpm.>13 & VPS$age_cat=="13 yrs to <18 yrs")
RespRate[is.na(VPS$high.respiratory.rate..bpm.)]<-NA
VPS$abnormal_RR <- RespRate

# Leukos
Leukocytes <- 
  (VPS$high.white.blood.cell.count..10.9..l.>34 & VPS$age_cat=="0 days to <1 wk") |
  (VPS$high.white.blood.cell.count..10.9..l.>19.5 & VPS$age_cat=="1 wk to <1 mo") |
  (VPS$low.white.blood.cell.count..10.9..l.<5 & VPS$age_cat=="1 wk to <1 mo") |
  (VPS$high.white.blood.cell.count..10.9..l.>17.5 & VPS$age_cat=="1 mo to <2 yrs") |
  (VPS$low.white.blood.cell.count..10.9..l.<5 & VPS$age_cat=="1 mo to <2 yrs") |
  (VPS$high.white.blood.cell.count..10.9..l.>15.5 & VPS$age_cat=="2 yrs to <6 yrs") |
  (VPS$low.white.blood.cell.count..10.9..l.<6 & VPS$age_cat=="2 yrs to <6 yrs") |
  (VPS$high.white.blood.cell.count..10.9..l.>13.5 & VPS$age_cat=="6 yrs to <13 yrs") |
  (VPS$low.white.blood.cell.count..10.9..l.<4.5 & VPS$age_cat=="6 yrs to <13 yrs") |
  (VPS$high.white.blood.cell.count..10.9..l.>11 & VPS$age_cat=="13 yrs to <18 yrs") |
  (VPS$low.white.blood.cell.count..10.9..l.<4.5 & VPS$age_cat=="13 yrs to <18 yrs") 
Leukocytes[is.na(VPS$high.white.blood.cell.count..10.9..l.)] <- NA
Leukocytes[is.na(VPS$low.white.blood.cell.count..10.9..l.)] <- NA
VPS$abnormal_Leukos <- Leukocytes

# SBP - Updated to fix the issue with the 50th instead of 5th centile!
# SysBloodPress <- 
#   (VPS$low.systolic.blood.pressure..mmhg.<65 & VPS$age_cat=="0 days to <1 wk")|
#   (VPS$low.systolic.blood.pressure..mmhg.<75 & VPS$age_cat=="1 wk to <1 mo") | 
#   (VPS$low.systolic.blood.pressure..mmhg.<100 & VPS$age_cat=="1 mo to <2 yrs") |
#   (VPS$low.systolic.blood.pressure..mmhg.<94 & VPS$age_cat=="2 yrs to <6 yrs") |
#   (VPS$low.systolic.blood.pressure..mmhg.<105 & VPS$age_cat=="6 yrs to <13 yrs") |
#   (VPS$low.systolic.blood.pressure..mmhg.<117 & VPS$age_cat=="13 yrs to <18 yrs")
SysBloodPress <- 
  (VPS$low.systolic.blood.pressure..mmhg.<59 & VPS$age_cat=="0 days to <1 wk")|
  (VPS$low.systolic.blood.pressure..mmhg.<79 & VPS$age_cat=="1 wk to <1 mo") | 
  (VPS$low.systolic.blood.pressure..mmhg.<75 & VPS$age_cat=="1 mo to <2 yrs") |
  (VPS$low.systolic.blood.pressure..mmhg.<74 & VPS$age_cat=="2 yrs to <6 yrs") |
  (VPS$low.systolic.blood.pressure..mmhg.<83 & VPS$age_cat=="6 yrs to <13 yrs") |
  (VPS$low.systolic.blood.pressure..mmhg.<90 & VPS$age_cat=="13 yrs to <18 yrs")
SysBloodPress[is.na(VPS$low.systolic.blood.pressure..mmhg.)]<-NA
VPS$abnormal_SBP <- SysBloodPress

# GCS is missing in 66,066 children!
VPS$abnormal_Mentation <- VPS$worst.glascow.coma.score <=11 # (per Goldstein)

# Order the age categories
VPS$age_cat <- factor(VPS$age_cat, levels=c("0 days to <1 wk","1 wk to <1 mo",
                      "1 mo to <2 yrs","2 yrs to <6 yrs","6 yrs to <13 yrs","13 yrs to <18 yrs"),ordered=TRUE)

# Cleanup
rm(Leukocytes,SysBloodPress,RespRate,Tachycardia,Bradycardia,CoreTemp)