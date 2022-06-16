require(skimr)
require(tidyverse)
require(tableone)
require(mgcv)
require(gridExtra)
require(grid)
require(zoo)

# Please set the working directory to run the code below
# The code saves the figures in the descendant directory named "Figures"
current_working_directory <- ''
setwd(current_working_directory)
df.Death <- read.csv('Sepsis_Death.csv')

df.Death %>% skim

# Create Table 1
df.tableone <- data.frame(df.Death)

binary_variables <- c('male', 'aspirin', 'diuretics', 'ace_inhibitors', 'ARBs', 'beta_blockers',
                      'ca_channel_blockers', 'PastHistory_Hypertension', 'PastHistory_Diabetes',
                      'PastHistory_COPD', 'PastHistory_Congestive_Heart_Failure', 'PastHistory_Peripheral_Vascular_Disease',
                      'PastHistory_Valve_disease', 'PastHistory_Pulmonary_Embolism', 'PastHistory_Neuromuscular_Disease',
                      'PastHistory_Hypothyroidism', 'PastHistory_Liver_Disease', 'PastHistory_AIDS',
                      'PastHistory_Cancer_Tumor', 'PastHistory_Arthritis_Vasculitis', 'PastHistory_Coagulopathy',
                      'PastHistory_Anemia', 'PastHistory_Home_Oxygen', 'PastHistory_Organ_Transplant',
                      'On_vent')

for (variable in binary_variables){
  df.tableone[, variable] <- df.tableone[, variable] %>% as.factor
}

covariates <- c(binary_variables, 'ethnicity', 'Admission_Type', 'Admit_Year', 'Hospital_Bed_Size', 'Hb',
                'Alb', 'WBCx1000', 'BUN', 'Lac', 'age', 'BMI', 'apachescore', 'NEE_mcg_perkg_perminx1000')
tableonevariables <- c(covariates, 'Systolic_below_120', 'Diastolic_below_120', 'Mean_below_120', 'Pulsepressure_below_120')

tab1 <- CreateTableOne(vars = tableonevariables, strata = "Death", data = df.tableone)
tab1Mat <- print(tab1, smd=TRUE)
write.csv(tab1Mat, file = "Sepsis_Table_1.csv")

get.df.analysis <- function(df){
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
  df$Systolic_below_120 <- as.integer(df$Systolic_below_120)
  df$Diastolic_below_120 <- as.integer(df$Diastolic_below_120)
  df$Mean_below_120 <- as.integer(df$Mean_below_120)
  df$Pulsepressure_below_120 <- as.integer(df$Pulsepressure_below_120)
  
  df_analysis <- df[complete.cases(df),]
  
  return(df_analysis)
}
df.Death.analysis <- get.df.analysis(df.Death)


# Generalized additive model analyses and predicted probabilities
get.predict.df <- function(df, component, outcome, percentile.range, n_bin, n_ma){
  percentile.list <- quantile(df[, component], c(percentile.range, 1-percentile.range))
  min.component <- round(percentile.list[[1]], 0)
  max.component <- round(percentile.list[[2]], 0)
  l <- max.component-min.component+1
  
  grid.range <- seq(min.component, max.component, length=l)
  
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
                           PastHistory_Organ_Transplant=rep(0, l), On_vent=rep(0, l),
                           ethnicity=rep('Caucasian', l), Admission_Type=rep('Emergency_Department', l),
                           Admit_Year=rep('2013-2016', l), Hospital_Bed_Size=rep('>500', l),
                           Hb=rep('>=11', l), Alb=rep('>=3', l), WBCx1000=rep('4<=WBCx1000<12', l),
                           BUN=rep('<=30', l), Lac=rep('<2', l), age=rep(56, l), BMI=rep(29.3, l),
                           apachescore=rep(60.0, l))
  
  predict.df[, component] <- grid.range
  
  covariates.factor <- c(binary_variables, 'ethnicity', 'Admission_Type', 'Admit_Year', 'Hospital_Bed_Size', 'Hb',
                         'Alb', 'WBCx1000', 'BUN', 'Lac')
  
  covformula <- paste(covariates.factor, collapse='+')
  gamformula <- as.formula(paste(outcome, '~', 's(', component, ')+', covformula, '+s(age)+s(BMI)+s(apachescore)', sep=''))
  
  gam.model <- gam(gamformula, data = df, family=binomial(link='logit'), method='REML')
  
  preds <- predict(gam.model, type='link', newdata=predict.df, se.fit=TRUE)
  fit <- preds$fit
  fit.low95 <- fit - 1.96*preds$se.fit
  fit.up95 <- fit + 1.96*preds$se.fit
  
  # Convert Log odds to probability
  predict.df$est <- 100*(exp(fit)/(1 + exp(fit)))
  predict.df$low95 <- 100*(exp(fit.low95)/(1 + exp(fit.low95)))
  predict.df$up95 <- 100*(exp(fit.up95)/(1 + exp(fit.up95)))
  
  #crude risk calculation using moving average
  df.ordered <- df[order(df[, component]), c(component, outcome)]
  n <- nrow(df.ordered)
  bin_size <- trunc(n/n_bin)
  r <- n%%bin_size
  
  m <- matrix(, nrow = n_bin, ncol = 2)
  s <- r%/%2
  for (i in 1:(n_bin)){
    index <- rep((s+1+(i-1)*bin_size):(s+i*bin_size))
    m[i, 1] <- mean(df.ordered[index, component])
    m[i, 2] <- 100*mean(df.ordered[index, outcome])
  }
  
  df.crude <- as.data.frame(m)
  colnames(df.crude) <- c(c('exposure', 'risk'))
  
  df.crude <- df.crude %>%
    mutate(ma = zoo::rollapply(df.crude$risk, width = 3L, FUN = mean, align = 'center', partial = TRUE))
  
  return(list(gam.model, predict.df, df.crude, min.component, max.component))
}

# Draw four spline curves with histogram
get.spline.hist <- function(df, predict.df, df.crude, component, exposure, min.component, max.component){
  
  if (component == 'Pulsepressure_below_120'){
    start <- floor(min.component/10)*10
    end <- ceiling(max.component/10)*10
    
    upper <- ceiling(max(df.crude$ma, predict.df$up95)/5)*5
    lower <- floor(min(df.crude$ma, predict.df$low95)/5)*5
    
    df_hist <- dplyr::count_(df, component, sort = FALSE)
    df.combined <- merge(x = predict.df, y = df_hist, by = c(component), all.x = TRUE)
    scale_val <- upper/max(df.combined[, 'n'])
    df.combined[, 'scaled_counts'] <- df.combined[, 'n']*scale_val
    
    spline <- ggplot() +
      xlab(paste('Lowest', exposure, '120-minutes', sep=' ')) +
      scale_x_continuous(breaks = seq(start, end, 10), limits = c(start, end)) +
      scale_y_continuous(breaks = seq(lower, upper, 5), limits = c(lower, upper)) +
      geom_bar(data = df.combined,
               aes_string(x = component, y='scaled_counts'),
               width = 1,
               stat='identity',
               color = 'black',
               size = 0.02,
               alpha = 0.15) +
      geom_line(data = df.combined, aes_string(x = component, y = 'est'), color='blue') +
      geom_ribbon(data = df.combined,
                  aes_string(x = component,
                             ymin = 'low95',
                             ymax = 'up95',
                             fill = 'up95' > 'low95'),
                  alpha = 0.5) +
      scale_fill_manual(values="skyblue") +
      theme_classic() +
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12)) +
      guides(fill = 'none') +
      geom_line(data = df.crude, aes(x=exposure, y=ma))
    
  } else if (component == 'Systolic_below_120'){
    
    start <- 50
    end <- 150
    
    upper <- ceiling(max(df.crude$ma, predict.df$up95)/10)*10
    lower <- floor(min(df.crude$ma, predict.df$low95)/10)*10
    
    df_hist <- dplyr::count_(df, component, sort = FALSE)
    df.combined <- merge(x = predict.df, y = df_hist, by = component, all.x = TRUE)
    scale_val <- upper/max(df.combined[, 'n'])
    df.combined[, 'scaled_counts'] <- df.combined[, 'n']*scale_val
    
    spline <- ggplot() +
      xlab(paste('Lowest', exposure, '120-minutes', sep=' ')) +
      scale_x_continuous(breaks = seq(start, end, 20), limits = c(start, end)) +
      scale_y_continuous(breaks = seq(lower, upper, 10), limits = c(lower, upper)) +
      geom_bar(data = df.combined,
               aes_string(x = component, y='scaled_counts'),
               width = 1,
               stat='identity',
               color = 'black',
               size = 0.02,
               alpha = 0.15) +
      geom_line(data = df.combined, aes_string(x = component, y = 'est'), color='blue') +
      geom_ribbon(data = df.combined,
                  aes_string(x = component,
                             ymin = 'low95',
                             ymax = 'up95',
                             fill = 'up95' > 'low95'),
                  alpha = 0.5) +
      scale_fill_manual(values="skyblue") +
      theme_classic() +
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            plot.margin = unit(c(9, 3, 3, 3), 'mm')) +
      guides(fill = 'none') +
      geom_line(data = df.crude, aes(x=exposure, y=ma))
    
  } else {
    
    start <- floor(min.component/10)*10
    end <- ceiling(max.component/10)*10
    
    upper <- ceiling(max(df.crude$ma, predict.df$up95)/10)*10
    lower <- floor(min(df.crude$ma, predict.df$low95)/10)*10
    
    df_hist <- dplyr::count_(df, component, sort = FALSE)
    df.combined <- merge(x = predict.df, y = df_hist, by = component, all.x = TRUE)
    scale_val <- upper/max(df.combined[, 'n'])
    df.combined[, 'scaled_counts'] <- df.combined[, 'n']*scale_val
    
    spline <- ggplot() +
      xlab(paste('Lowest', exposure, '120-minutes', sep=' ')) +
      scale_x_continuous(breaks = seq(start, end, 10), limits = c(start, end)) +
      scale_y_continuous(breaks = seq(lower, upper, 10), limits = c(lower, upper)) +
      geom_bar(data = df.combined,
               aes_string(x = component, y='scaled_counts'),
               width = 1,
               stat='identity',
               color = 'black',
               size = 0.02,
               alpha = 0.15) +
      geom_line(data = df.combined, aes_string(x = component, y = 'est'), color='blue') +
      geom_ribbon(data = df.combined,
                  aes_string(x = component,
                             ymin = 'low95',
                             ymax = 'up95',
                             fill = 'up95' > 'low95'),
                  alpha = 0.5) +
      scale_fill_manual(values="skyblue") +
      theme_classic() +
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            plot.margin = unit(c(9, 3, 3, 3), 'mm')) +
      guides(fill = 'none') +
      geom_line(data = df.crude, aes(x=exposure, y=ma))
  }
  
  return(spline)
}

get.spline.hist.from.df <- function(df, component, exposure, outcome, percentile.range, n_bin, n_ma){
  list.component <- get.predict.df(df, component, outcome, percentile.range, n_bin, n_ma)
  gam.model <- list.component[[1]]
  predict.df <- list.component[[2]]
  df.crude <- list.component[[3]]
  min.component <- list.component[[4]]
  max.component <- list.component[[5]]
  
  #spline curve for risk of outcome using mean values for covariates
  spline <- get.spline.hist(df, predict.df, df.crude, component, exposure, min.component, max.component)
  
  return(spline)
}

# Combine all plots to one figure
get.four.plots.hist <- function(df, outcome, outcome_description, percentile.range, n_bin=100, n_ma=5){
  
  spline.mean <- get.spline.hist.from.df(df, 'Mean_below_120', 'MAP', outcome, percentile.range, n_bin, n_ma)
  spline.systolic <- get.spline.hist.from.df(df, 'Systolic_below_120', 'SBP', outcome, percentile.range, n_bin, n_ma)
  spline.diastolic <- get.spline.hist.from.df(df, 'Diastolic_below_120', 'DBP', outcome, percentile.range, n_bin, n_ma)
  spline.pulse <- get.spline.hist.from.df(df, 'Pulsepressure_below_120', 'PP', outcome, percentile.range, n_bin, n_ma)
  
  all.spline <- grid.arrange(
    spline.mean,
    spline.systolic,
    spline.diastolic,
    spline.pulse,
    nrow = 2,
    top=textGrob(paste('Predicted Probability of ', outcome_description, (' (%)'), sep=''),
                 gp=gpar(fontsize=18))
  )
  
  return(all.spline)
}

all.spline.hist.death <- get.four.plots.hist(df.Death.analysis, 'Death', 'Mortality', 0.001)
all.spline.hist.death
ggsave("./Figures/Sepsis_spline_hist_death.png", all.spline.hist.death, height = 8, width = 12, dpi=600)

# AKI
df.AKI <- read.csv('Sepsis_AKI_composite.csv')
sum(df.AKI$AKI)
mean(df.AKI$AKI)
df.AKI.analysis <- get.df.analysis(df.AKI)

all.spline.hist.AKI <- get.four.plots.hist(df.AKI.analysis, 'AKI_composite', 'AKI or Death', 0.001)
all.spline.hist.AKI
ggsave("./Figures/Sepsis_spline_hist_AKI.png", all.spline.hist.AKI, height = 8 , width = 12, dpi=600)

# Myocardial Injury
df.MI <- read.csv('Sepsis_MI_composite.csv')
sum(df.MI$MI)
mean(df.MI$MI)
df.MI.analysis <- get.df.analysis(df.MI)

all.spline.hist.MI <- get.four.plots.hist(df.MI.analysis, 'MI_composite', 'Myocardial Injury or Death', 0.001)
all.spline.hist.MI
ggsave("./Figures/Sepsis_spline_hist_MI.png", all.spline.hist.MI, height = 8 , width = 12, dpi=600)
