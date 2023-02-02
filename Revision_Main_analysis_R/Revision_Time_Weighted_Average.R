require(skimr)
require(tidyverse)
require(mgcv)
require(gridExtra)
require(grid)
require(caTools)
#require(zoo)

#Set the working directory
#Saving the figures in the next directory named "Figures"
input_dir <- '/Users/takahirokiritoshi/OneDrive/Hypotension_Project/Analysis/Revision/'
output_dir <- '/Users/takahirokiritoshi/OneDrive/Hypotension_Project/Analysis/Revision/Figures/'
df.TWA<- read.csv(paste0(input_dir, 'Sepsis_TWA_revision.csv'))
df.TWA %>% skim

binary_variables <- c('male', 'aspirin', 'diuretics', 'ace_inhibitors', 'ARBs', 'beta_blockers',
                      'ca_channel_blockers', 'PastHistory_Hypertension', 'PastHistory_Diabetes',
                      'PastHistory_COPD', 'PastHistory_Congestive_Heart_Failure', 'PastHistory_Peripheral_Vascular_Disease',
                      'PastHistory_Valve_disease', 'PastHistory_Pulmonary_Embolism', 'PastHistory_Neuromuscular_Disease',
                      'PastHistory_Hypothyroidism', 'PastHistory_Liver_Disease', 'PastHistory_AIDS',
                      'PastHistory_Cancer_Tumor', 'PastHistory_Arthritis_Vasculitis', 'PastHistory_Coagulopathy',
                      'PastHistory_Anemia', 'PastHistory_Home_Oxygen', 'PastHistory_Organ_Transplant',
                      'PastHistory_Chronic_Kidney_Disease', 'PastHistory_Myocardial_Infarction',
                      'PastHistory_Stroke', 'PastHistory_Coronary_Artery_Disease',
                      'On_vent')

#generalized additive model analyses and predicted probabilities
get.df.analysis.TWA <- function(df){
  df <- df %>% mutate_if(is.character, as.factor)
  df$ethnicity <- relevel(df$ethnicity, ref='Caucasian')
  df$Admission_Type <- relevel(df$Admission_Type, ref='Emergency_Department')
  df$Admit_Year <- relevel(df$Admit_Year, ref='2013-2016')
  df$Hospital_Bed_Size <- relevel(df$Hospital_Bed_Size, ref='>500')
  df$Hb <- relevel(df$Hb, ref='>=11')
  df$Alb <- relevel(df$Alb, ref='>=3')
  df$WBCx1000 <- relevel(df$WBCx1000, ref='4<=WBCx1000<12')
  df$BUN <- relevel(df$BUN, ref='<=30')
  df$Lac <- relevel(df$Lac, ref='<2')
  df_analysis <- df[complete.cases(df),]
  
  return(df_analysis)
}

df.TWA.analysis <- get.df.analysis.TWA(df.TWA)

get.predict.df.TWA <- function(df, component, percentile.range){
  
  components <- c('TWA_MAP', 'TWA_SBP', 'TWA_DBP', 'TWA_PP')
  ind <- which(components==component)
  max.components <- c(20, 25, 25, 35)
  max.component <- max.components[[ind]]
  min.component <- 0
  l <- 2*(max.component-min.component)+1
  
  grid.range <- seq(min.component, max.component, by=0.5)
  
  predict.df <- data.frame(male=rep(1, l), aspirin=rep(0, l), diuretics=rep(0, l), ace_inhibitors=rep(0, l),
                           ARBs=rep(0, l), beta_blockers=rep(0, l), ca_channel_blockers=rep(0, l),
                           PastHistory_Hypertension=rep(0, l), PastHistory_Diabetes=rep(0, l),
                           PastHistory_COPD=rep(0, l), PastHistory_Congestive_Heart_Failure=rep(0, l),
                           PastHistory_Peripheral_Vascular_Disease=rep(0, l), PastHistory_Valve_disease=rep(0, l),
                           PastHistory_Pulmonary_Embolism=rep(0, l), PastHistory_Neuromuscular_Disease=rep(0, l),
                           PastHistory_Hypothyroidism=rep(0, l), PastHistory_Liver_Disease=rep(0, l),
                           PastHistory_AIDS=rep(0, l), PastHistory_Cancer_Tumor=rep(0, l),
                           PastHistory_Arthritis_Vasculitis=rep(0, l), PastHistory_Coagulopathy=rep(0, l),
                           PastHistory_Anemia=rep(0, l), PastHistory_Home_Oxygen=rep(0, l),
                           PastHistory_Organ_Transplant=rep(0, l), PastHistory_Chronic_Kidney_Disease=rep(0, l),
                           PastHistory_Myocardial_Infarction=rep(0, l), PastHistory_Stroke=rep(0, l),
                           PastHistory_Coronary_Artery_Disease=rep(0, l), On_vent=rep(0, l),
                           ethnicity=rep('Caucasian', l), Admission_Type=rep('Emergency_Department', l),
                           Admit_Year=rep('2013-2016', l), Hospital_Bed_Size=rep('>500', l),
                           Hb=rep('>=11', l), Alb=rep('>=3', l), WBCx1000=rep('4<=WBCx1000<12', l),
                           BUN=rep('<=30', l), Lac=rep('<2', l), age=rep(56, l), BMI=rep(29.3, l),
                           apachescore=rep(60.0, l))
  
  predict.df[, component] <- grid.range
  
  covariates.factor <- c(binary_variables, 'ethnicity', 'Admission_Type', 'Admit_Year', 'Hospital_Bed_Size', 'Hb',
                         'Alb', 'WBCx1000', 'BUN', 'Lac')
  
  covformula <- paste(covariates.factor, collapse='+')
  gamformula <- as.formula(paste0('Death ~ ', 's(', component, ')+', covformula, '+s(age)+s(BMI)+s(apachescore)'))
  
  gam.model <- gam(gamformula, data = df, family=binomial(link='logit'), method='REML')
  
  preds <- predict(gam.model, type='link', newdata=predict.df, se.fit=TRUE)
  fit <- preds$fit
  fit.low95 <- fit - 1.96*preds$se.fit
  fit.up95 <- fit + 1.96*preds$se.fit
  
  # Convert Log odds to probability
  predict.df$est <- 100*(exp(fit)/(1 + exp(fit)))
  predict.df$low95 <- 100*(exp(fit.low95)/(1 + exp(fit.low95)))
  predict.df$up95 <- 100*(exp(fit.up95)/(1 + exp(fit.up95)))
  
  return(list(gam.model, predict.df, min.component, max.component))
}

#Draw four spline curves with histogram
get.spline.hist.TWA <- function(df, predict.df, component, min.component, max.component){
  
  exposure <- str_remove(component, 'TWA_')
  exposures <- c('MAP', 'SBP', 'DBP', 'PP')
  thresholds <- c(69, 100, 60, 57)
  ind <- which(exposures == exposure)
  threshold <- thresholds[[ind]]
  
  start <- min.component
  end <- max.component
    
  upper <- 20
  lower <- 0
  count_col <- paste0(exposure, '_count')
  df[, count_col] <- 5*round(df[, component]/5, 1)
  df_hist <- dplyr::count_(df, count_col, sort = FALSE)
  predict.df$ID <- as.factor(predict.df[, component])
  df_hist$ID <- as.factor(df_hist[, count_col])
  df.combined <- merge(x = predict.df, y = df_hist, by = c('ID'), all.x = TRUE)
  scale_val <- upper/max(na.omit(df.combined[, 'n']))
  df.combined[, 'scaled_counts'] <- df.combined[, 'n']*scale_val
    
  spline <- ggplot() +
    xlab(paste('Time-Weighted Average of', exposure, '<', threshold, '(mmHg)', sep=' ')) +
    scale_x_continuous(breaks = seq(start, end, 5), limits = c(start, end)) +
    scale_y_continuous(breaks = seq(lower, upper, 10), limits = c(lower, upper)) +
    geom_bar(data = df.combined,
             aes_string(x = component, y='scaled_counts'),
             width = 0.5,
             position = position_nudge(x = 0.25),
             stat = 'identity',
             color = 'black',
             size = 0.02,
             alpha = 0.15) +
    geom_line(data = df.combined, aes_string(x = component, y = 'est'), color='red') +
    geom_ribbon(data = df.combined,
                aes_string(x = component,
                           ymin = 'low95',
                           ymax = 'up95',
                           fill = 'up95' > 'low95'),
                alpha = 0.5) +
    scale_fill_manual(values="lightpink") +
    theme_classic() +
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 9),
          plot.margin = unit(c(9, 3, 3, 3), 'mm')) +
    guides(fill = 'none')
  return(spline)
}

get.spline.hist.from.df.TWA <- function(df, component, percentile.range){
  list.component <- get.predict.df.TWA(df, component, percentile.range)
  gam.model <- list.component[[1]]
  predict.df <- list.component[[2]]
  min.component <- list.component[[3]]
  max.component <- list.component[[4]]
  
  fit <- predict(gam.model, type='link', newdata=df, se.fit=FALSE)
  
  # Convert Log odds to probability
  df$est <- 100*(exp(fit)/(1 + exp(fit)))
  auc <- round(as.numeric(caTools::colAUC(df$est, df$Death)), 3)
  
  #spline curve for risk of outcome using mean values for covariates
  spline <- get.spline.hist.TWA(df, predict.df, component, min.component, max.component)
  
  return(list(spline, auc, gam.model))
}

get.four.plots.hist.TWA <- function(df, percentile.range){
  
  mean <- get.spline.hist.from.df.TWA(df, 'TWA_MAP', percentile.range)
  systolic <- get.spline.hist.from.df.TWA(df, 'TWA_SBP', percentile.range)
  diastolic <- get.spline.hist.from.df.TWA(df, 'TWA_DBP', percentile.range)
  pulse <- get.spline.hist.from.df.TWA(df, 'TWA_PP', percentile.range)
  
  spline.mean <- mean[[1]]
  spline.systolic <- systolic[[1]]
  spline.diastolic <- diastolic[[1]]
  spline.pulse <- pulse[[1]]
  
  df.auc <- data.frame(component = c('MAP', 'SBP', 'DBP', 'PP'),
                       AUC = c(mean[[2]],systolic[[2]], diastolic[[2]], pulse[[2]])) %>% t
  
  models <- list(mean[[3]], systolic[[3]], diastolic[[3]], pulse[[3]])
  
  all.spline <- grid.arrange(
    spline.mean,
    spline.systolic,
    spline.diastolic,
    spline.pulse,
    nrow = 2,
    top=textGrob('Predicted Probability of ICU Mortality in Overall Sepsis Cohort (%)',
                 gp=gpar(fontsize=14))
  )
  
  return(list(all.spline, df.auc, models))
}

results <- get.four.plots.hist.TWA(df.TWA.analysis, 0.001)
df.auc <- results[[2]]
df.auc

gammodels <- results[[3]]

components <- c('MAP', 'SBP', 'DBP', 'PP')
Chi.sq.values <- as.numeric(c())
for (component in components){
  ind <- which(components == component)
  model <- gammodels[[ind]]
  mod.sum <- summary(model)
  Chi.sq.values[[ind]] <- mod.sum$s.table[1,3]
}

all.spline.hist.death.TWA <- results[[1]]
all.spline.hist.death.TWA

ggsave("./Figures/Fig_4.tiff", all.spline.hist.death.TWA, device='tiff', height = 12, width = 18, units='cm', dpi=300)

ggsave("./Figures/Sepsis_TWA_all.png", all.spline.hist.death.TWA, height = 12, width = 18, units='cm', dpi=300)

