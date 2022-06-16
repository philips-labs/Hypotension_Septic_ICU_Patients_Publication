library(tidyverse)
library(magrittr)
library(splines)
library(gridExtra)
library(grid)

# Please set the working directory to run the code below
input_dir <- ''
output_dir <- paste0(input_dir, 'Figures/')

get.df.analysis <- function(df){
  df <- df %>% dplyr::mutate_if(is.character, as.factor)
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

binary_variables <- c('male', 'aspirin', 'diuretics', 'ace_inhibitors', 'ARBs',
                      'beta_blockers', 'ca_channel_blockers',
                      'PastHistory_Hypertension', 'PastHistory_Diabetes',
                      'PastHistory_COPD', 'PastHistory_Congestive_Heart_Failure',
                      'PastHistory_Peripheral_Vascular_Disease',
                      'PastHistory_Valve_disease',
                      'PastHistory_Pulmonary_Embolism',
                      'PastHistory_Neuromuscular_Disease',
                      'PastHistory_Hypothyroidism', 'PastHistory_Liver_Disease',
                      'PastHistory_AIDS', 'PastHistory_Cancer_Tumor',
                      'PastHistory_Arthritis_Vasculitis',
                      'PastHistory_Coagulopathy', 'PastHistory_Anemia',
                      'PastHistory_Home_Oxygen', 'PastHistory_Organ_Transplant',
                      'On_vent')


covariates.factor <- c(binary_variables, 'ethnicity', 'Admission_Type',
                       'Admit_Year', 'Hospital_Bed_Size', 'Hb',
                       'Alb', 'WBCx1000', 'BUN', 'Lac')

covformula <- paste(covariates.factor, collapse='+')

# Read csv file for each outcome
read_in_data <- function(outcome){
  df_init <- read.csv(file.path(input_dir, paste0('Sepsis_', outcome, '.csv')))
  return(df_init)
}

# Conduct threshold regression
get_thresh_model <- function(component, outcome, covformula, seed=123,
                             subsample_frac=1){
  df_init <- read_in_data(outcome)
  df_analysis <- get.df.analysis(df_init)
  if (subsample_frac < 1){
    set.seed(seed)
    df_analysis <- dplyr::sample_frac(df_analysis, size=subsample_frac)
  }
  cat('Start', outcome, 'and', component, '\n')
  
  formula.1 <- as.formula(paste(outcome, '~', covformula,
                                '+ns(age)+ns(BMI)+ns(apachescore)', sep=''))
  formula.2 <- as.formula(paste0('~', component))
  
  start <- Sys.time()
  thresh.model <- chngpt::chngptm(formula.1=formula.1, formula.2=formula.2,
                                  df_analysis, type='M10', family='binomial', var.type='none')
  end <- Sys.time()
  et <- difftime(end, start, units='hours')
  if (et < 1){
    et <- paste(round(difftime(end, start, units='min'), 3), 'minutes')
  } else {
    et <- paste(round(et, 3), 'hours')
  }
  cat('End', outcome, 'and', component, 'in', et, '\n')
  return(thresh.model)
}

outcomes <- c('MI_composite', 'AKI_composite', 'Death')
components <- c('Mean_below_120', 'Systolic_below_120',
                'Diastolic_below_120', 'Pulsepressure_below_120')
combn_df <- expand.grid(outcome=outcomes, component=components,
                        stringsAsFactors=FALSE)
run_parallel <- function(i, combn_df, covformula, subsample_frac=1){
  outcome <- combn_df$outcome[i]
  component <- combn_df$component[i]
  mdl <- get_thresh_model(component=component, outcome=outcome,
                          covformula=covformula, subsample_frac=subsample_frac)
  return(list(model=mdl, outcome=outcome, component=component))
}

thresh_mdls_init <- parallel::mclapply(seq_len(nrow(combn_df)), run_parallel,
                                        combn_df=combn_df, covformula=covformula,
 				       subsample_frac=1, mc.cores=6)
thresh_mdls <- list()
for (outcome in outcomes){
  thresh_mdls[[outcome]] <- list()
  for (component in components){
    i <- lapply(thresh_mdls_init,
                function(e)e$outcome == outcome && e$component == component) %>%
           unlist() %>%
      which()
    thresh_mdls[[outcome]][[component]] <- thresh_mdls_init[[i]]
  }
}

cat('Saving models...')
save(thresh_mdls, file=file.path(output_dir, 'thresh_mdls.RData'))
cat('Done.\n')

load(file.path(output_dir, 'thresh_mdls.RData'))

cat('Getting threshold results...')
thresh_results_df <- dplyr::tibble(
  outcome=as.character(c()), component=as.character(c()),
  changepoint=as.numeric(c()), effect_size=as.numeric(c()))
for (outcome in outcomes){
  for (component in components){
    mdl <- thresh_mdls[[outcome]][[component]]$model
    df <- dplyr::tibble(outcome=outcome, component=component,
                        changepoint=mdl$chngpt,
                        effect_size=tail(coef(mdl), 1))
    thresh_results_df <- rbind(thresh_results_df, df)
  }
}

readr::write_tsv(thresh_results_df,
                 file.path(output_dir, 'threshold_analysis_results.txt'))
cat('Done.\n')

get.predict.df <- function(df, component, outcome, percentile.range){
  percentile.list <- quantile(df[, component],
                              c(percentile.range, 1-percentile.range))
  min.component <- round(percentile.list[[1]], 0)
  max.component <- round(percentile.list[[2]], 0)
  l <- max.component-min.component+1
  
  grid.range <- seq(min.component, max.component, length=l)
  
  predict.df <- data.frame(
    male=rep(1, l), aspirin=rep(0, l), diuretics=rep(0, l),
    ace_inhibitors=rep(0, l), ARBs=rep(0, l), beta_blockers=rep(0, l),
    ca_channel_blockers=rep(0, l), PastHistory_Hypertension=rep(0, l),
    PastHistory_Diabetes=rep(0, l),
    PastHistory_COPD=rep(0, l), PastHistory_Congestive_Heart_Failure=rep(0, l),
    PastHistory_Peripheral_Vascular_Disease=rep(0, l),
    PastHistory_Valve_disease=rep(0, l),
    PastHistory_Pulmonary_Embolism=rep(0, l),
    PastHistory_Neuromuscular_Disease=rep(0, l),
    PastHistory_Hypothyroidism=rep(0, l), PastHistory_Liver_Disease=rep(0, l),
    PastHistory_AIDS=rep(0, l), PastHistory_Cancer_Tumor=rep(0, l),
    PastHistory_Arthritis_Vasculitis=rep(0, l),
    PastHistory_Coagulopathy=rep(0, l), PastHistory_Anemia=rep(0, l),
    PastHistory_Home_Oxygen=rep(0, l), PastHistory_Organ_Transplant=rep(0, l),
    On_vent=rep(0, l), ethnicity=rep('Caucasian', l),
    Admission_Type=rep('Emergency_Department', l),
    Admit_Year=rep('2013-2016', l), Hospital_Bed_Size=rep('>500', l),
    Hb=rep('>=11', l), Alb=rep('>=3', l), WBCx1000=rep('4<=WBCx1000<12', l),
    BUN=rep('<=30', l), Lac=rep('<2', l), age=rep(56, l), BMI=rep(29.8, l),
    apachescore=rep(60.9, l))
  
  predict.df[, component] <- grid.range
  
  preds <- predict(thresh_mdls[[outcome]][[component]]$model, type='link',
                   newdata=predict.df, se.fit=TRUE)
  fit <- preds$fit
  
  # Convert Log odds to probability
  bp_mean <- mean(df[[component]])
  bp_sd <- sd(df[[component]])
  prediction_df <- dplyr::tibble(
    bp=grid.range,
    bp_standardized=(grid.range - bp_mean) / bp_sd,
    predicted_risk=100*(exp(fit)/(1 + exp(fit))),
    outcome=outcome,
    component=component)
  
  return(prediction_df)
}

cat('Plotting results...')

plt_df_list <- list()
for (outcome in outcomes){
  df_init <- read_in_data(outcome)
  df_analysis <- get.df.analysis(df_init)
  for (component in components){
    plt_df_list[[length(plt_df_list) + 1]] <- get.predict.df(
      df=df_analysis, component=component, outcome=outcome,
      percentile.range=0.001)
  }
}

plt_df <- do.call(rbind, plt_df_list) %>%
  dplyr::mutate(
    outcome=ifelse(outcome == 'Death', outcome, paste0(stringr::str_replace(outcome, '_', ' ('), ')')),
    outcome=factor(outcome,
                   levels=c('Death', 'AKI (composite)', 'MI (composite)')),
    component=factor(component,
                     levels=c('Mean_below_120', 'Systolic_below_120',
                              'Diastolic_below_120', 'Pulsepressure_below_120'),
                     labels=c('MAP', 'SBP', 'DBP', 'PP')))

threshold_plots <- list()
plot_outcomes = c('Death', 'AKI (composite)', 'MI (composite)')
for (plot_outcome in plot_outcomes){
  ind <- which(plot_outcomes == plot_outcome)
  subset_plt_df <- dplyr::filter(plt_df, outcome==plot_outcome)
  if (ind != 3) {
    threshold_plots[[ind]] <- ggplot(data=subset_plt_df, aes(x=bp_standardized, y=predicted_risk, color=component)) +
      geom_line() +
      scale_y_continuous(breaks = seq(0, 12.5, 2.5), limits = c(0, 12.5)) +
      scale_color_brewer(name=NULL, palette = 'Set1') +
      labs(title = paste0(plot_outcome, '\n')) +
      theme_classic() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 8, hjust=0.5),
            legend.position = "none")
  } else {
    threshold_plots[[ind]] <- ggplot(data=subset_plt_df, aes(x=bp_standardized, y=predicted_risk, color=component)) +
      geom_line() +
      scale_y_continuous(breaks = seq(0, 12.5, 2.5), limits = c(0, 12.5)) +
      scale_color_brewer(name=NULL, palette = 'Set1') +
      labs(title = 'Myocardial injury\n(composite)') +
      theme_classic() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 8, hjust=0.5),
            legend.text = element_text(size = 8))
  }
}

# Save the figure for threshold logistic regression
png(file.path(output_dir, paste0('threshold_final.png')), height=4, width=7.2,
    units='in', res=600)

threshold_plot <- grid.arrange(
  threshold_plots[[1]],
  threshold_plots[[2]],
  threshold_plots[[3]],
  nrow = 1,
  widths = c(1, 1, 1.43),
  left=textGrob('Predicted Probability (%)',
                gp=gpar(fontsize=10),
                rot=90),
  bottom=textGrob('Blood pressure component (standardized)',
                  gp=gpar(fontsize=10))
)
dev.off()

cat('Done.\n')