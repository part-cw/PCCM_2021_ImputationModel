

##########################################################################################################

# Pre-processing VPS dataset

VPS_preprocessing <- function(the_VPS_df){
  the_VPS_df$worst.glascow.coma.score.unknown <- is.na(the_VPS_df$worst.glascow.coma.score)
  the_VPS_df$low.white.blood.cell.count.unknown <- is.na(the_VPS_df$low.white.blood.cell.count..10.9..l.)
  the_VPS_df$high.white.blood.cell.count.unknown <- is.na(the_VPS_df$high.white.blood.cell.count..10.9..l.)
  
  missing_low_not_missing_high <- the_VPS_df$low.white.blood.cell.count.unknown & (!the_VPS_df$high.white.blood.cell.count.unknown)
  missing_high_not_missing_low <- the_VPS_df$high.white.blood.cell.count.unknown & (!the_VPS_df$low.white.blood.cell.count.unknown)
  
  ## Missing low white blood cell count replaced with recorded high white blood cell count
  the_VPS_df[missing_low_not_missing_high,"low.white.blood.cell.count..10.9..l."] <- the_VPS_df[missing_low_not_missing_high,"high.white.blood.cell.count..10.9..l."]
  ## Missing high white blood cell count replaced with recorded low white blood cell count
  the_VPS_df[missing_high_not_missing_low,"high.white.blood.cell.count..10.9..l."] <- the_VPS_df[missing_high_not_missing_low,"low.white.blood.cell.count..10.9..l."]
  the_VPS_df$low.white.blood.cell.count.unknown <- is.na(the_VPS_df$low.white.blood.cell.count..10.9..l.)
  the_VPS_df$high.white.blood.cell.count.unknown <- is.na(the_VPS_df$high.white.blood.cell.count..10.9..l.)
  return(the_VPS_df)
}



##########################################################################################################



##########################################################################################################
# USED FOR PLOTTING

plot_theme <- function(a_ggplot, the_angle=0, x_angle=0, y_angle){
  a_ggplot <- a_ggplot + theme(axis.text = element_text(size = 9, angle = the_angle), 
                               axis.text.x = element_text(size=9, angle = x_angle),
                               axis.text.y = element_text(size =9, angle = y_angle),
                               axis.title = element_text(face = "bold", size = 12),
                               legend.text = element_text(size = 12),
                               legend.title = element_text(size = 12, face = "bold"),
                               strip.text = element_text(face = "bold", size=12),  
                               plot.title = element_text(face="bold", size = 14, hjust =0.5), 
                               plot.caption = element_text(size=10, face = "bold", hjust = 0))
  return(a_ggplot)
}

##########################################################################################################


##########################################################################################################


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

##########################################################################################################
