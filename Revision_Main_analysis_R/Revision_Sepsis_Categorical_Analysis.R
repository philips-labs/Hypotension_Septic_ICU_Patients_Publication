require(tidyverse)

#Set the working directory
setwd('/Users/takahirokiritoshi/OneDrive/Hypotension_Project/Analysis/Revision')
df.Death <- read.csv('Sepsis_Death_revision.csv')
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
covariates.binary <- c('male', 'aspirin', 'diuretics', 'ace_inhibitors', 'ARBs', 'beta_blockers',
                       'ca_channel_blockers', 'PastHistory_Hypertension', 'PastHistory_Diabetes',
                       'PastHistory_COPD', 'PastHistory_Congestive_Heart_Failure', 'PastHistory_Peripheral_Vascular_Disease',
                       'PastHistory_Valve_disease', 'PastHistory_Pulmonary_Embolism', 'PastHistory_Neuromuscular_Disease',
                       'PastHistory_Hypothyroidism', 'PastHistory_Liver_Disease', 'PastHistory_AIDS',
                       'PastHistory_Cancer_Tumor', 'PastHistory_Arthritis_Vasculitis', 'PastHistory_Coagulopathy',
                       'PastHistory_Anemia', 'PastHistory_Home_Oxygen', 'PastHistory_Organ_Transplant',
                       'PastHistory_Chronic_Kidney_Disease', 'PastHistory_Myocardial_Infarction',
                       'PastHistory_Stroke', 'PastHistory_Coronary_Artery_Disease',
                       'On_vent')
covariates.factor <- c(covariates.binary, 'ethnicity', 'Admission_Type', 'Admit_Year', 'Hospital_Bed_Size', 'Hb',
                       'Alb', 'WBCx1000', 'BUN', 'Lac')
covariates.all <- c(covariates.factor, 'age', 'BMI', 'apachescore')

get.categorical.df <- function(df, component, threshold_ref, threshold_exp){
  
  var_comp <- as.name(component)
  
  if (component=='Pulsepressure_below_120'){
    
    if (threshold_exp>=50){
      
      df.categorical <- df %>% filter(((!!var_comp) >= threshold_ref & (!!var_comp) < threshold_ref + 10) | (!!var_comp) >= threshold_exp)
      df.categorical <- mutate(df.categorical, exposure = if_else((!!var_comp) >= threshold_exp, true = 1, false = 0))
      
    } else {
      
      df.categorical <- df %>% filter(((!!var_comp) >= threshold_ref & (!!var_comp) < threshold_ref + 10) | (!!var_comp) < threshold_exp)
      df.categorical <- mutate(df.categorical, exposure = if_else((!!var_comp) < threshold_exp, true = 1, false = 0))
      
    }
    
  } else {
    
    df.categorical <- df %>% filter((!!var_comp) >= threshold_ref | (!!var_comp) < threshold_exp)
    df.categorical <- mutate(df.categorical, exposure = if_else((!!var_comp) >= threshold_ref, true = 0, false = 1))
    
  }
  
  return(df.categorical)
  
}

get.OR <- function(df, component, threshold_ref, threshold_exp, outcome){
  
  df.categorical <- get.categorical.df(df, component, threshold_ref, threshold_exp)
  var_outcome <- as.name(outcome)
  
  n_ref <- nrow(df.categorical %>% filter(exposure==0))
  event_ref <- nrow(df.categorical %>% filter(exposure==0 & (!!var_outcome)==1))
  rate_ref <- round(100*event_ref/n_ref, 1)
  
  n_exp <- nrow(df.categorical %>% filter(exposure==1))
  event_exp <- nrow(df.categorical %>% filter(exposure==1 & (!!var_outcome)==1))
  rate_exp <- round(100*event_exp/n_exp, 1)
  
  crude <- c(n_ref, paste(event_ref, ' (', rate_ref, '%)', sep = ''),
             n_exp, paste(event_exp, ' (', rate_exp, '%)', sep = ''))
  
  form <- as.formula(paste(outcome, '~', 'exposure+', paste(covariates.all, collapse = '+'), sep=''))
  
  model <- glm(form, data = df.categorical, family = 'binomial')
  
  est <- coef(summary(model))[2, 1]
  se <- coef(summary(model))[2, 2]
  pval <- coef(summary(model))[2, 4]
  
  OR <- c(round(exp(est), 2), paste(round(exp(est - 1.96*se), 2), 'to', round(exp(est + 1.96*se), 2), sep=' '),
          round(pval, 3))
  
  return(list(crude, OR))
  #return(list(df.categorical, crude))
}

get.Mean.results <- function(df){
  m <- matrix(, nrow = 5, ncol = 5)
  thresholds <- seq(75, 45, by = -10)
  for (threshold in thresholds){
    i <- which(thresholds == threshold)
    if (i==1){
      results <- get.OR(df.Death.analysis, 'Mean_below_120', 75, threshold, 'Death')
      m[1, 1:2] <- results[[1]][1:2]
      m[i+1, 1:2] <- results[[1]][3:4]
      m[i+1, 3:5] <- results[[2]]
    } else {
      results <- get.OR(df.Death.analysis, 'Mean_below_120', 75, threshold, 'Death')
      m[i+1, 1:2] <- results[[1]][3:4]
      m[i+1, 3:5] <- results[[2]]
    }
  }
  df.results <- as.data.frame(m, row.names = c('Mean BP >=75mmHg', 'Mean BP <75mmHg', 'Mean BP <65mmHg', 'Mean BP <55mmHg', 'Mean BP <45mmHg'))
  colnames(df.results) <- c('n', 'Death (%)', 'Adjusted_OR', '95% CI', 'p-value')
  return(df.results)
}
get.Systolic.results <- function(df){
  m <- matrix(, nrow = 5, ncol = 5)
  thresholds <- seq(110, 80, by = -10)
  for (threshold in thresholds){
    i <- which(thresholds == threshold)
    if (i==1){
      results <- get.OR(df.Death.analysis, 'Systolic_below_120', 110, threshold, 'Death')
      m[1, 1:2] <- results[[1]][1:2]
      m[i+1, 1:2] <- results[[1]][3:4]
      m[i+1, 3:5] <- results[[2]]
    } else {
      results <- get.OR(df.Death.analysis, 'Systolic_below_120', 110, threshold, 'Death')
      m[i+1, 1:2] <- results[[1]][3:4]
      m[i+1, 3:5] <- results[[2]]
    }
  }
  df.results <- as.data.frame(m, row.names = c('sBP >=110mmHg', 'sBP <110mmHg', 'sBP <100mmHg', 'sBP <90mmHg', 'sBP <80mmHg'))
  colnames(df.results) <- c('n', 'Death (%)', 'Adjusted_OR', '95% CI', 'p-value')
  return(df.results)
}
get.Diastolic.results <- function(df){
  m <- matrix(, nrow = 5, ncol = 5)
  thresholds <- seq(60, 30, by = -10)
  for (threshold in thresholds){
    i <- which(thresholds == threshold)
    if (i==1){
      results <- get.OR(df.Death.analysis, 'Diastolic_below_120', 60, threshold, 'Death')
      m[1, 1:2] <- results[[1]][1:2]
      m[i+1, 1:2] <- results[[1]][3:4]
      m[i+1, 3:5] <- results[[2]]
    } else {
      results <- get.OR(df.Death.analysis, 'Diastolic_below_120', 60, threshold, 'Death')
      m[i+1, 1:2] <- results[[1]][3:4]
      m[i+1, 3:5] <- results[[2]]
    }
  }
  df.results <- as.data.frame(m, row.names = c('dBP >=60mmHg', 'dBP <60mmHg', 'dBP <50mmHg', 'dBP <40mmHg', 'dBP <30mmHg'))
  colnames(df.results) <- c('n', 'Death (%)', 'Adjusted_OR', '95% CI', 'p-value')
  return(df.results)
}
get.Pulse.results <- function(df){
  m <- matrix(, nrow = 5, ncol = 5)
  thresholds <- seq(50, 20, by = -10)
  for (threshold in thresholds){
    i <- which(thresholds == threshold)
    if (i==1){
      results <- get.OR(df.Death.analysis, 'Pulsepressure_below_120', 40, threshold, 'Death')
      m[1, 1:2] <- results[[1]][1:2]
      m[i+1, 1:2] <- results[[1]][3:4]
      m[i+1, 3:5] <- results[[2]]
    } else {
      results <- get.OR(df.Death.analysis, 'Pulsepressure_below_120', 40, threshold, 'Death')
      m[i+1, 1:2] <- results[[1]][3:4]
      m[i+1, 3:5] <- results[[2]]
    }
  }
  df.results <- as.data.frame(m, row.names = c('40 <= PP < 50mmHg', 'PP >=50mmHg', 'PP <40mmHg', 'PP <30mmHg', 'PP <20mmHg'))
  colnames(df.results) <- c('n', 'Death (%)', 'Adjusted_OR', '95% CI', 'p-value')
  return(df.results)
}

df.results.MBP <- get.Mean.results(df.Death.analysis)
df.results.sBP <- get.Systolic.results(df.Death.analysis)
df.results.dBP <- get.Diastolic.results(df.Death.analysis)
df.results.PP <- get.Pulse.results(df.Death.analysis)

results.table <- rbind(df.results.MBP, df.results.sBP, df.results.dBP, df.results.PP)
write.csv(results.table, 'Sepsis_Table_2.csv')

df.subgroup <- df.Death.analysis %>% 
  mutate(age_cat = case_when(age<45 ~ '<45',
                             (age>=45 & age<65) ~ '45 <=age <65',
                             age>=65 ~ '>65'),
         vasopressor = case_when(NEE_mcg_perkg_perminx1000==0 ~ 'none',
                                 (NEE_mcg_perkg_perminx1000>0 & NEE_mcg_perkg_perminx1000<=50) ~ '0 <Average Rate (mcg/kg/min) <=0.05',
                                 (NEE_mcg_perkg_perminx1000>50 & NEE_mcg_perkg_perminx1000<=100) ~ '0.05 <Average Rate (mcg/kg/min) <=0.1',
                                 NEE_mcg_perkg_perminx1000>100 ~ 'Average Rate (mcg/kg/min) >0.1'),
         exposure = if_else(Mean_below_120>=65, 0, 1))
df.subgroup$age_cat <- factor(df.subgroup$age_cat, levels = c('<45', '45 <=age <65', '>65'))
df.subgroup$vasopressor <- factor(df.subgroup$vasopressor, levels = c('none', '0 <Average Rate (mcg/kg/min) <=0.05', '0.05 <Average Rate (mcg/kg/min) <=0.1', 'Average Rate (mcg/kg/min) >0.1'))
df.subgroup$Admission_Type <- factor(df.subgroup$Admission_Type, levels = c('Emergency_Department', 'Other_Ward', 'Elective', 'Other_Hospital', 'Other/Unknown'))

#potential.covariates <- c('age', 'On_vent', 'Admission_Type', 'PastHistory_Cancer_Tumor', 'Lac', 'apachescore', 'Hb', 'Alb', 'WBCx1000')
potential.covariates <- c('age', 'On_vent', 'Admission_Type', 'PastHistory_Cancer_Tumor', 'Lac', 'apachescore')


correlation.check <- df.subgroup[,c('age', 'On_vent', 'PastHistory_Cancer_Tumor', 'apachescore')]
colnames(correlation.check) <- c('Age', 'Ventilator', 'Cancer/Tumor', 'APACHE IVa')
cormat <- round(cor(correlation.check), 2)

library(reshape2)
melted_cormat <- melt(cormat)
cormat.graph <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank())
cormat.graph

ggsave("./Figures/Correlation_Matrix.png", cormat.graph, dpi=600)

library(car)
vif.model <- glm(Death ~ age + On_vent + Admission_Type + PastHistory_Cancer_Tumor + Lac + apachescore,
                 data=df.subgroup, family = binomial())
vif.result <- tibble::rownames_to_column(as.data.frame(vif(vif.model)), 'Variable')
vif.result
vif.graph <- ggplot(vif.result, aes(y=GVIF, x=Variable)) +
  geom_col(fill='steelblue') +
  geom_hline(yintercept=5, linetype='dashed', color = 'blue') +
  coord_flip() +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5)) +
  ggtitle('Generalized Variance Inflation Factor') +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 18, hjust = 0.5))

vif.graph
ggsave("./Figures/VIF.png", vif.graph, height = 4 , width = 12, dpi=600)

get.subgroup.OR <- function(df, outcome, subgroup, potential.covariates){
  var_sub <- as.name(subgroup)
  var_outcome <- as.name(outcome)
  if (subgroup == 'age_cat') {
    covariates <- potential.covariates[-which(potential.covariates=='age')]
    levs <- levels(df[,subgroup])
  } else if (subgroup == 'vasopressor') {
    covariates <- potential.covariates
    levs <- levels(df[,subgroup])
  } else if (subgroup == 'Admission_Type'){
    covariates <- potential.covariates[-which(potential.covariates==subgroup)]
    levs <- levels(df[,subgroup])
    levs <- levs[-which(levs=='Other/Unknown')]
  } else {
    covariates <- potential.covariates[-which(potential.covariates==subgroup)]
    levs <- c(0, 1)
  }
  sub.formula <- as.formula(paste(outcome, '~exposure+', paste(covariates, collapse = '+'), sep=''))
  n_subgroups <- length(levs)
  
  subgroup_names <- vector(, length = n_subgroups)
  m <- matrix(, nrow = n_subgroups, ncol = 7)
  
  for (lev in levs){
    ind <- which(levs == lev)
    
    subgroup_names[ind] <- paste(subgroup, lev, sep = ' ')
    
    df.subanalysis <- df %>% filter((!!var_sub)==lev)
    exp <- df.subanalysis %>% filter(exposure == 1)
    exp_event <- nrow(exp %>% filter((!!var_outcome)==1))
    con <- df.subanalysis %>% filter(exposure == 0)
    con_event <- nrow(con %>% filter((!!var_outcome)==1))
    
    m[ind, 1] <- nrow(exp)
    m[ind, 2] <- paste(exp_event, ' (', round(100*exp_event/nrow(exp), 1), '%)', sep = '')
    m[ind, 3] <- nrow(con)
    m[ind, 4] <- paste(con_event, ' (', round(100*con_event/nrow(con), 1), '%)', sep = '')
    
    model <- glm(sub.formula, data = df.subanalysis, family = 'binomial')
    
    est <- coef(summary(model))[2, 1]
    se <- coef(summary(model))[2, 2]
    pval <- coef(summary(model))[2, 4]
    
    m[ind, 5] <- round(exp(est), 2)
    m[ind, 6] <- paste(round(exp(est-1.96*se), 2), ' to ', round(exp(est+1.96*se), 2), sep='')
    m[ind, 7] <- round(pval, 3)
  }
  
  df.results <- as.data.frame(m, row.names = subgroup_names)
  colnames(df.results) <- c('Mean BP <=65mmHg (n)', 'Death (%)', 'Mean BP >65mmHg (n)', 'Death (%)', 'Adjusted_OR', '95% CI', 'p-value')
  
  return(df.results)
}
results.subgroup.age <- get.subgroup.OR(df.subgroup, 'Death', 'age_cat', potential.covariates)
results.subgroup.vasopressor <- get.subgroup.OR(df.subgroup, 'Death', 'vasopressor', potential.covariates)
results.subgroup.vent <- get.subgroup.OR(df.subgroup, 'Death', 'On_vent', potential.covariates)
results.subgroup.cancer <- get.subgroup.OR(df.subgroup, 'Death', 'PastHistory_Cancer_Tumor', potential.covariates)
results.subgroup.admtype <- get.subgroup.OR(df.subgroup, 'Death', 'Admission_Type', potential.covariates)

overall.formula <- as.formula(paste('Death', '~exposure+', paste(potential.covariates, collapse = '+'), sep=''))
overall.model <- glm(overall.formula, data = df.subgroup, family = 'binomial')
overall.est <- coef(summary(overall.model))[2, 1]
overall.se <- coef(summary(overall.model))[2, 2]
overall.pval <- coef(summary(overall.model))[2, 4]
n_exposed <- nrow(df.subgroup %>% filter(exposure==1))
n_outcome_exposed <- nrow(df.subgroup %>% filter(exposure==1 & Death==1))
n_unexposed <- nrow(df.subgroup %>% filter(exposure==0))
n_outcome_unexposed <- nrow(df.subgroup %>% filter(exposure==0 & Death==1))
overall.vector <- c(n_exposed, paste(n_outcome_exposed, ' (', round(100*n_outcome_exposed/n_exposed, 1), '%)', sep=''),
                    n_unexposed, paste(n_outcome_unexposed, ' (', round(100*n_outcome_unexposed/n_unexposed, 1), '%)', sep=''),
                    round(exp(overall.est), 2), paste(round(exp(overall.est-1.96*overall.se), 2), ' to ', round(exp(overall.est+1.96*overall.se), 2), sep=''), round(overall.pval, 3))
results.subgroup.overall <- as.data.frame(t(overall.vector), row.names = 'Overall')
colnames(results.subgroup.overall) <- c('Mean BP <=65mmHg (n)', 'Death (%)', 'Mean BP >65mmHg (n)', 'Death (%)', 'Adjusted_OR', '95% CI', 'p-value')

subresults.table <- rbind(results.subgroup.age, results.subgroup.vasopressor, results.subgroup.vent,
                          results.subgroup.cancer, results.subgroup.admtype, results.subgroup.overall)
write.csv(subresults.table, 'Sepsis_Table_3.csv')