library(tidyverse)
library(magrittr)
library(splines)
library(gridExtra)
library(grid)
library(skimr)

input_dir <- '/Users/takahirokiritoshi/OneDrive/Hypotension_Project/Analysis/Revision/'
output_dir <- '/Users/takahirokiritoshi/OneDrive/Hypotension_Project/Analysis/Revision/Figures/'

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


covariates.factor <- c(binary_variables, 'ethnicity', 'Admission_Type',
                       'Admit_Year', 'Hospital_Bed_Size', 'Hb',
                       'Alb', 'WBCx1000', 'BUN', 'Lac')

covformula <- paste(covariates.factor, collapse='+')

components <- c('Mean_below_120', 'Systolic_below_120',
                'Diastolic_below_120', 'Pulsepressure_below_120')

# Read csv file for each outcome
df <- read.csv(paste0(input_dir, 'Sepsis_Death_revision.csv'))
df.analysis <- get.df.analysis(df)

df.shock <- data.frame(df.analysis)
df.shock <- df.shock %>% 
  mutate(septic_shock = if_else(((Lac=='2<=Lac<5' | Lac=='>=5') & NEE_mcg_perkg_perminx1000>0), 1, 0))
df.shock <- df.shock %>% filter(septic_shock==1)
df.shock %>% skim

covariates.factor.shock <- c('male', 'On_vent', 'ethnicity', 'Hb', 'Alb', 'WBCx1000', 'BUN')
covformula.shock<- paste(covariates.factor.shock, collapse='+')

# Conduct threshold regression
get_thresh_model <- function(component, outcome, df.analysis, covformula, seed=123,
                             subsample_frac=1){
  if (subsample_frac < 1){
    set.seed(seed)
    df.analysis <- dplyr::sample_frac(df.analysis, size=subsample_frac)
  }
  cat('Start', outcome, 'and', component, '\n')
  
  formula.1 <- as.formula(paste(outcome, '~', covformula,
                                '+ns(age)+ns(BMI)+ns(apachescore)', sep=''))
  formula.2 <- as.formula(paste0('~', component))
  
  start <- Sys.time()
  thresh.model <- chngpt::chngptm(formula.1=formula.1, formula.2=formula.2,
                                  df.analysis, type='M10', family='binomial', var.type='none')
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


run_parallel <- function(i, components, df.analysis, covformula, subsample_frac=1){
  component <- components[i]
  mdl <- get_thresh_model(component=component, outcome='Death', df.analysis,
                          covformula=covformula, subsample_frac=subsample_frac)
  return(list(model=mdl, component=component))
}

thresh_mdls_init_shock <- parallel::mclapply(seq_len(4), run_parallel, df.analysis=df.shock,
                                       components=components, covformula=covformula.shock,
                                       subsample_frac=1, mc.cores=4)
thresh_mdls_shock <- list()
for (component in components){
  i <- lapply(thresh_mdls_init_shock,
              function(e)e$component == component) %>%
    unlist() %>%
    which()
  thresh_mdls_shock[[component]] <- thresh_mdls_init_shock[[i]]
}

cat('Saving models...')
save(thresh_mdls_shock, file=file.path(output_dir, 'thresh_mdls_shock.RData'))
cat('Done.\n')

load(file.path(output_dir, 'thresh_mdls_shock.RData'))

cat('Getting threshold results...')
thresh_results_df_shock <- dplyr::tibble(
  component=as.character(c()),
  changepoint=as.numeric(c()), effect_size=as.numeric(c()))
for (component in components){
  mdl_shock <- thresh_mdls_shock[[component]]$model
  df <- dplyr::tibble(component=component,
                      changepoint=mdl_shock$chngpt,
                      effect_size=tail(coef(mdl_shock), 1))
  thresh_results_df_shock <- rbind(thresh_results_df_shock, df)
}


readr::write_tsv(thresh_results_df_shock,
                 file.path(output_dir, 'threshold_analysis_results_shock.txt'))
cat('Done.\n')

get.predict.df <- function(df, component, thresh_mdls, percentile.range){
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
    PastHistory_Chronic_Kidney_Disease=rep(0, l), PastHistory_Myocardial_Infarction=rep(0, l),
    PastHistory_Stroke=rep(0, l), PastHistory_Coronary_Artery_Disease=rep(0, l),
    On_vent=rep(0, l), ethnicity=rep('Caucasian', l),
    Admission_Type=rep('Emergency_Department', l),
    Admit_Year=rep('2013-2016', l), Hospital_Bed_Size=rep('>500', l),
    Hb=rep('>=11', l), Alb=rep('>=3', l), WBCx1000=rep('4<=WBCx1000<12', l),
    BUN=rep('<=30', l), Lac=rep('<2', l), age=rep(56, l), BMI=rep(29.8, l),
    apachescore=rep(60.9, l))
  
  predict.df[, component] <- grid.range
  
  preds <- predict(thresh_mdls[[component]]$model, type='link',
                   newdata=predict.df, se.fit=TRUE)
  fit <- preds$fit
  
  # Convert Log odds to probability
  bp_mean <- mean(df[[component]])
  bp_sd <- sd(df[[component]])
  prediction_df <- dplyr::tibble(
    bp=grid.range,
    bp_standardized=(grid.range - bp_mean) / bp_sd,
    predicted_risk=100*(exp(fit)/(1 + exp(fit))),
    component=component)
  
  return(prediction_df)
}

cat('Plotting results...')
plt_df_list.shock <- list()
for (component in components){
  plt_df_list.shock[[length(plt_df_list.shock) + 1]] <- get.predict.df(
    df=df.shock, component=component, thresh_mdls_shock,
    percentile.range=0.001)
}

plt_df.shock <- do.call(rbind, plt_df_list.shock) %>%
  dplyr::mutate(
    component=factor(component,
                     levels=c('Pulsepressure_below_120', 'Diastolic_below_120',
                              'Systolic_below_120', 'Mean_below_120'),
                     labels=c('PP', 'DBP', 'SBP', 'MAP')))

#threshold_plots <- list()

threshold_plots.shock <- ggplot(data=plt_df.shock, aes(x=bp_standardized, y=predicted_risk,
                                                 color=component, linetype=component)) +
  geom_line(size=0.5) +
  scale_color_manual(name='BP Components', values=c('#BC3C29EF', '#0072B5FF', '#E18727FF', '#20854EFF')) +
  scale_linetype_manual(name='BP Components', values=c('dashed', 'twodash', 'longdash', 'solid')) +
  scale_y_continuous(breaks = seq(0, 60, 10), limits = c(0, 60)) +
  labs(title = 'Septic Shock Patients\n(n=4,211)') +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(size = 8, hjust=0.5),
        legend.title = element_text(size=7),
        legend.text = element_text(size=7),
        legend.position = c(0.6, 0.7),
        legend.key.height = unit(0.2, 'cm'),
        legend.key.width = unit(1.68, 'cm'))+
  guides(color = guide_legend(reverse=TRUE), linetype = guide_legend(reverse = TRUE))

tiff(file.path(output_dir, paste0('Fig_3.b.tiff')), height=12, width=9,
     units='cm', res=300)

threshold_plots.shock

dev.off()

cat('Done.\n')  


















# Entire Data
#library(caret)
#train.index <- createDataPartition(df.analysis$Death, p = .3, list = FALSE)
#df.analysis <- df.analysis[train.index,]

thresh_mdls_init <- parallel::mclapply(seq_len(4), run_parallel, df.analysis=df.analysis,
                                       components=components, covformula=covformula.shock,
                                       subsample_frac=1, mc.cores=4)
thresh_mdls <- list()
for (component in components){
  i <- lapply(thresh_mdls_init,
              function(e)e$component == component) %>%
    unlist() %>%
    which()
  thresh_mdls[[component]] <- thresh_mdls_init[[i]]
}

cat('Saving models...')
save(thresh_mdls, file=file.path(output_dir, 'thresh_mdls_local.RData'))
cat('Done.\n')

load(file.path(output_dir, 'thresh_mdls_local.RData'))

cat('Getting threshold results...')
thresh_results_df <- dplyr::tibble(
  component=as.character(c()),
  changepoint=as.numeric(c()), effect_size=as.numeric(c()))
for (component in components){
  mdl <- thresh_mdls[[component]]$model
  df <- dplyr::tibble(component=component,
                      changepoint=mdl$chngpt,
                      effect_size=tail(coef(mdl), 1))
  thresh_results_df <- rbind(thresh_results_df, df)
}

readr::write_tsv(thresh_results_df,
                 file.path(output_dir, 'threshold_analysis_results.txt'))
cat('Done.\n')

cat('Plotting results...')
plt_df_list<- list()
for (component in components){
  plt_df_list[[length(plt_df_list) + 1]] <- get.predict.df(
    df=df.analysis, component=component, thresh_mdls,
    percentile.range=0.001)
}

plt_df <- do.call(rbind, plt_df_list) %>%
  dplyr::mutate(
    component=factor(component,
                     levels=c('Pulsepressure_below_120', 'Diastolic_below_120',
                              'Systolic_below_120', 'Mean_below_120'),
                     labels=c('PP', 'DBP', 'SBP', 'MAP')))

#threshold_plots <- list()

threshold_plot <- ggplot(data=plt_df, aes(x=bp_standardized, y=predicted_risk,
                                                       color=component, linetype=component)) +
  geom_line(size=0.5) +
  scale_color_manual(name='BP Components', values=c('#BC3C29EF', '#0072B5FF', '#E18727FF', '#20854EFF')) +
  scale_linetype_manual(name='BP Components', values=c('dashed', 'twodash', 'longdash', 'solid')) +
  scale_y_continuous(breaks = seq(0, 60, 10), limits = c(0, 60)) +
  labs(title = 'Overall Sepsis Cohort\n(n=77,328)') +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(size = 8, hjust=0.5),
        legend.title = element_text(size=7),
        legend.text = element_text(size=7),
        legend.position = c(0.6, 0.7),
        legend.key.height = unit(0.2, 'cm'),
        legend.key.width = unit(1.68, 'cm')) +
  guides(color = guide_legend(reverse=TRUE), linetype = guide_legend(reverse = TRUE))

tiff(file.path(output_dir, paste0('Fig_3_a.tiff')), height=12, width=9,
     units='cm', res=300)

threshold_plot

dev.off()

library(ggpmisc)
library(patchwork)

table.all <- data.frame('Component' = c('MAP', 'SBP', 'DBP', 'PP'),
                        'col_1' = c(62-2*10, 94-2*14, 48-2*9, 35-2*12),
                        'col_2' = c(62, 94, 48, 35),
                        'col_3' = c(62+2*10, 94+2*14, 48+2*9, 35+2*12),
                        'col_4' = c(62+4*10, 94+4*14, 48+4*9, 35+4*12))
table.all <- table.all %>% 
  select('Component', 'col_1', 'col_2', 'col_3', 'col_4') %>%  
  setNames(c('Component', '-2SD', '0(mean)', '+2SD', '+4SD'))

table.shock <- data.frame('Component' = c('MAP', 'SBP', 'DBP', 'PP'),
                        'col_1' = c(56-2.5*8, 84-2.5*10, 44-2.5*8, 28-2.5*10),
                        'col_2' = c(56, 84, 44, 28),
                        'col_3' = c(56+2.5*8, 84+2.5*10, 44+2.5*8, 28+2.5*10))
table.shock <- table.shock %>% 
  select('Component', 'col_1', 'col_2', 'col_3') %>% 
  setNames(c('Component', '-2.5SD', '0(mean)', '+2.5SD'))


#ggp_table_all <- ggplot() +
#  theme_void() +
#  annotate(geom = 'table',
#           x=1,
#           y=1,
#           size=1,
#           label = list(table.all))
  
#ggp_table_shock <- ggplot() +
#  theme_void() +
#  annotate(geom = 'table',
#           x=1,
#           y=1,
#           size=1,
#           label = list(table.shock))

tt <- ttheme_minimal(
  core = list(fg_params=list(fontsize = 6), padding=unit(c(1.2, 1.6), 'mm')),
  colhead = list(fg_params=list(fontsize = 6, fontface = 1), padding=unit(c(1.2, 1.6), 'mm')),
  rowhead = list(fg_params=list(fontsize = 6)))

tbl1 <- tableGrob(table.all, rows=NULL, theme=tt)

separators <- list(segmentsGrob(x0 = unit(0, 'npc'),
                                y0 = unit(1, 'npc'),
                                x1 = unit(1, 'npc'),
                                y1 = unit(1, 'npc'),
                                gp = gpar(lwd=1)),
                   segmentsGrob(x0 = unit(0, 'npc'),
                                y0 = unit(1, 'npc'),
                                x1 = unit(1, 'npc'),
                                y1 = unit(1, 'npc'),
                                gp = gpar(lwd=1)),
                   segmentsGrob(x0 = unit(0, 'npc'),
                                y0 = unit(0, 'npc'),
                                x1 = unit(1, 'npc'),
                                y1 = unit(0, 'npc'),
                                gp = gpar(lwd=1)))

tbl1 <- gtable::gtable_add_grob(tbl1,
                                grobs=separators,
                                t=c(1,2,5), b=c(1,2,5), l=1, r=5)
tbl2 <- tableGrob(table.shock, rows=NULL, theme=tt)
tbl2 <- gtable::gtable_add_grob(tbl2,
                                grobs=separators,
                                t=c(1,2,5), b=c(1,2,5), l=1, r=4)


tbl <- grid.arrange(tbl1, tbl2, ncol=2, widths=c(1.2, 1))

threshold_plots <- grid.arrange(
  threshold_plot,
  threshold_plots.shock,
  nrow = 1,
  widths = c(1, 1.02),
  left=textGrob('Predicted Probability of ICU Mortality (%)',
                gp=gpar(fontsize=9),
                rot=90),
  bottom=textGrob('Blood pressure component (standardized)',
                  gp=gpar(fontsize=9))
)

tiff(file.path(output_dir, paste0('Fig_3.tiff')), height=11, width=9,
     units='cm', res=300)

thrresh_plot_table <- grid.arrange(threshold_plots, tbl,
                                   nrow=2,
                                   heights=c(4.2,1))

dev.off()

ggsave(paste0(output_dir, 'Threshold_analysis.png'), thrresh_plot_table, height = 10, width = 9, units='cm', dpi=300)

threshold_mmHg_plots <- list()
comps <- c('MAP', 'SBP', 'DBP', 'PP')
line_colors <- c('#BC3C29EF', '#0072B5FF', '#E18727FF', '#20854EFF')

for (component in components){
  ind <- which(components==component)
  df <- plt_df %>% filter(component==comps[[ind]])
  threshold_mmHg_plots[[ind]] <- ggplot() +
    geom_line(data=df, aes(x=bp, y=predicted_risk), color=line_colors[[ind]]) +
    xlab(paste('Lowest', comps[[ind]], '120-minutes (mmHg)', sep=' ')) +
    scale_y_continuous(breaks = seq(0, 40, 10), limits = c(0, 40)) +
    theme_classic() +
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 9),
          plot.margin = unit(c(9, 3, 3, 3), 'mm')) +
    guides(fill = 'none')
}

tiff(file.path(output_dir, paste0('./Supplementary/Additional_File_6.tiff')), height=12, width=18,
     units='cm', res=300)

threshold_mmHg_plot <- grid.arrange(
  threshold_mmHg_plots[[1]],
  threshold_mmHg_plots[[2]],
  threshold_mmHg_plots[[3]],
  threshold_mmHg_plots[[4]],
  nrow = 2,
  top=textGrob('Predicted ICU Mortality by Threshold Logistic Regression (%)',
               gp=gpar(fontsize=14)))

dev.off()
cat('Done.\n') 