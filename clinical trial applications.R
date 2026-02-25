## This R code allows to replicate the two simulation studies presented in the paper 
## "Using permutation-based tests to evaluate an adaptive molecular treatment
##  algorithm in randomized precision oncology trials" by Marelli et al.

## Install and upload packages and functions

# install.packages("readxl")
library(readxl)
# install.packages("dplyr")
library(dplyr)
# install.packages("survival")
library(survival)
# install.packages("survminer")
library(survminer)
# install.packages("lubridate")
library(lubridate)
# install.packages("Minirand")
library(Minirand)
#install.packages("furrr")
library(furrr)
# install.packages("tidyverse")
library(tidyverse)

## User-defined functions

## First randomization of the treatment labels and permutation test statistic calculation:
## the minimization procedure is re-run and the stratified log-rank test statistic is re-computed.
## Arguments:
## * 'data' is the dataset excluding patients treated with the discontinued treatments at both stages. 
## It must include the treatment label ('armAB' equal to "A" or "B"), the time of randomization ('date_rand'),
## the indicator variable of the event ('event'), and the ## time the event happens ('time_event');
## * 'cov_list' is the list of variables names to considere as covariates;
## * all the other entries follow the 'Minirand' function from the package 'Minirand'.
simulate_one_rand_EXT <- function(data, cov_list, ntrt, ratio, covwt) {
  
  ## controls
  if (is.null(data)) 
    stop("Missing dataset")
  if (!"armAB" %in% names(data)) 
    stop("The dataset must contains a column 'armAB'")
  if (!all(data$armAB %in% c("A", "B")))
      stop("Column 'armAB' must contain only values 'A' or 'B'")
  if (!"date_rand" %in% names(data)) 
    stop("The dataset must contains a column 'date_rand'")
  if (!"event" %in% names(data)) 
    stop("The dataset must contains a column 'event'")
  if (!"time_event" %in% names(data)) 
    stop("The dataset must contains a column 'time_event'")
  
  N <- nrow(data)
  
  covmat <- as.matrix(data[, cov_list])
  
  res <- rep(NA, N)
  res[1] <- sample(1:ntrt, 1)
  
  for(j in 2:N) {
    res[j] <- Minirand(
      covmat = covmat,
      j = j,
      covwt = covwt,
      ratio = ratio,
      ntrt = ntrt,
      trtseq = c(1,0),
      method = "Range",
      result = res,
      p = 0.9
    )
  }
  
  data_temp <- data
  data_temp$trt_rand <- NULL
  data_temp <- data_temp %>%
    arrange(date_rand)
  data_temp$trt_rand <- res
  data_temp$armAB_rand <- ifelse(data_temp$trt_rand == 1,"A","B")
  data_temp$armAB_rand <- as.factor(data_temp$armAB_rand)
  data_temp$armAB_rand <- relevel(data_temp$armAB_rand, ref = "B")
  
  data_temp <- data_temp %>%
    arrange(date_rand)
  PFS_surv_obj_temp <- Surv(time = data_temp$time_event, event = data_temp$event)
  
  strata_vars <- colnames(covmat)
  strata_terms <- paste0("strata(`", strata_vars, "`)", collapse = " + ")
  formula_str <- paste("PFS_surv_obj_temp ~ armAB_rand +", strata_terms)
  formula_cox <- as.formula(formula_str)
  
  S_rand_s <- survdiff(formula_cox, data = data_temp)$chisq
  
  
  tibble::tibble(S_rand_strat = S_rand_s) 
}

## First randomization of the treatment labels and permutation test statistic calculation:
## the minimization procedure is re-run on stage I data and the Z-score for the log hazard ratio 
## is compared to an interim futility threshold. If the trial continues to stage II, 
## the stratified log-rank test statistic is computed.
## Arguments:
## * 'data' is the final dataset (stage I + stage II). It must include the enrollment period ('stage'), 
## the time of randomization ('date_rand'), the indicator variable of the event ('event') 
## and the time the event happens ('time_event');
## * 'target' is the name of the molecular alteration on which to perform the interim futility analysis;
## * 'cov_list' is the list of variables names to considere as covariates;
## * 'crit_val' is the interim futility threshold;
## * 'L' is the number of valid randomization needed;
## * 'max_attempts' is the maximum number of attempts to run to find the 'L' required valid randomizations;
## * all the other entries follow the 'Minirand' function from the package 'Minirand'.
simulate_one_rand_IA <- function(data, target, cov_list, crit_val=0, L=15000, 
                                 max_attempts=100000, ntrt, ratio, covwt) {
  
  ## controls
  if (is.null(data)) 
    stop("Missing dataset")
  if (!"armAB" %in% names(data))
    stop("The dataset must contains a column 'armAB'")
  if (!all(data$armAB %in% c("A", "B")))
    stop("Column 'armAB' must contain only values 'A' or 'B'")
  if (!"date_rand" %in% names(data)) 
    stop("The dataset must contains a column 'date_rand'")
  if (!"event" %in% names(data)) 
    stop("The dataset must contains a column 'event'")
  if (!"time_event" %in% names(data)) 
    stop("The dataset must contains a column 'time_event'")
  if (is.null(target)) 
    stop("Missing molecular alteration")
  if (is.null(cov_list)) 
    stop("Missing list of covariates")
  if (length(crit_val)!=1 || !is.numeric(crit_val) || is.nan(crit_val)) 
    stop("'crit_val' must be numeric or infinit")
  if (length(L)!=1 || !is.numeric(L) || L <= 0 || L != round(L)) 
    stop("'L' must be a positive integer")
  if (length(max_attempts)!=1 || !is.numeric(max_attempts) || is.nan(max_attempts)) 
    stop("'max_attempts' must be numeric or infinit")
  
  
  # STAGE I
  
  data_stageI <- data[data$stage==1,]
  
  N1 <- nrow(data_stageI)
  
  covmat1 <- as.matrix(data_stageI[, cov_list])
  
  S_rand_s <- rep(0,L)
  I_L_s <- rep(0,L)
  l <- 1  # valid randomization counter
  attempts <- 0 # number of attempts counter
  
  while(l <= L){
    attempts <- attempts + 1
    
    if(attempts > max_attempts){
      warning(paste("Only ", l-1, " randomizations obtained out of ", L, " randomizations"))
      S_rand_s <- S_rand_s[1:(l-1)]
      I_L_s <- I_L_s[1:(l-1)]
      break
    }
    
    res1 <- rep(NA, N1)
    res1[1] <- sample(1:ntrt, 1)
    
    for(j in 2:N1) {
      res1[j] <- Minirand(
        covmat = covmat1,
        j = j,
        covwt = covwt,
        ratio = ratio,
        ntrt = ntrt,
        trtseq = c(1,0),
        method = "Range",
        result = res1,
        p = 0.9
      )
    }
    
    data_stageI$trt_rand <- res1
    data_stageI$armAB_rand <- ifelse(data_stageI$trt_rand == 1, "A", "B")
    data_stageI$armAB_rand <- as.factor(data_stageI$armAB_rand)
    data_stageI$armAB_rand <- relevel(data_stageI$armAB_rand, ref = "B")
    
    data_interim <- data_stageI[data_stageI$TARGET==target,]
    if (nrow(data_interim) == 0)
      stop(paste("No rows found for TARGET =", target))
    
    data_interim <- data_interim %>%
      arrange(date_rand)
    PFS_surv_obj_interim <- Surv(time = data_interim$time_event, event = data_interim$event)
    
    strata_vars1 <- colnames(covmat1)
    strata_terms1 <- paste0("strata(`", strata_vars1, "`)", collapse = " + ")
    formula_str1 <- paste("PFS_surv_obj_interim ~ armAB_rand +", strata_terms1)
    formula_cox1 <- as.formula(formula_str1)
    
    cox_interim <- coxph(formula_cox1, data = data_interim)
    
    logHR_interim <- coef(cox_interim)        
    SE_logHR_interim <- sqrt(diag(vcov(cox_interim))) 
    
    Z_score_interim <- logHR_interim / SE_logHR_interim
    
    if(Z_score_interim<crit_val){
      
      data_stageII <- data[data$stage==2,]
      
      N2 <- nrow(data_stageII)
      
      covmat2 <- as.matrix(data_stageII[, cov_list])
      
      res2 <- rep(NA, N2)
      res2[1] <- sample(1:ntrt, 1)
      
      for(j in 2:N2) {
        res2[j] <- Minirand(
          covmat = covmat2,
          j = j,
          covwt = covwt,
          ratio = ratio,
          ntrt = ntrt,
          trtseq = c(1,0),
          method = "Range",
          result = res2,
          p = 0.9
        )
      }
      
      data_stageII$trt_rand <- res2
      data_stageII$armAB_rand <- ifelse(data_stageII$trt_rand == 1, "A", "B")
      
      data_rand <- rbind(data_stageI,data_stageII)
      data_rand$armAB_rand <- as.factor(data_rand$armAB_rand)
      data_rand$armAB_rand <- relevel(data_rand$armAB_rand, ref = "B")
      
      data_rand <- data_rand %>%
        arrange(date_rand)
      PFS_surv_obj_rand <- Surv(time = data_rand$time_event, event = data_rand$event)
      
      strata_vars2 <- colnames(covmat2)
      strata_terms2 <- paste0("strata(`", strata_vars2, "`)", collapse = " + ")
      formula_str2 <- paste("PFS_surv_obj_rand ~ armAB_rand +", strata_terms2)
      formula_cox2 <- as.formula(formula_str2)
      
      S_rand_s[l] <- survdiff(formula_cox2, data = data_rand)$chisq
      
      l <- l+1
    }
  }
  
  tibble::tibble(S_rand_strat=S_rand_s,total_attempts = attempts[1])
}


################################################################################


## Load the synthetic dataset generated from the SAFIR02-BREAST/PI3K dataset.
## It includes the following variables:
## * "targeted_therapy": therapy tailored to the characteristics of the patient;
## * "armAB": treatment arm ('A' or 'B');
## * "breast_molec_alt": breast molecular alteration;
## * "chemo_line": line of chemotherapy at the time of randomization ;      
## * "disease_status": disease status at randomization;
## * "escat_class": escat classification;
## * "date_rand": date of randomization;
## * "PFS": follow-up;              
## * "PFS_event": event (1=event, 0=censored).      

syn_BREAST <- read_excel("C:/Users/C_MARELLI/OneDrive - gustaveroussy.fr/PhD PROJECT/papers/axis 1/github code/syn_BREAST.xlsx")
head(syn_BREAST)
nrow(syn_BREAST)

syn_BREAST$armAB <- as.factor(syn_BREAST$armAB)
syn_BREAST$armAB <- relevel(syn_BREAST$armAB, ref = "B")
syn_BREAST$event <- as.numeric(syn_BREAST$PFS_event)
syn_BREAST$time_event <- as.numeric(syn_BREAST$PFS)

## ESCAT I/II population
db_escat <- syn_BREAST[syn_BREAST$escat_class %in% c("I","II"),]
nrow(db_escat)

# EXTERNAL MODIFICATIONS:
# drop of: AZD8931 (and AZD4547, AZD2014, VANDETANIB)
db_escat %>%
  filter(armAB == "A") %>%
  with(table(targeted_therapy))

# Stratified Cox models with AZD8931, AZD4547, AZD2014, and VANDETANIB
table(db_escat$armAB,db_escat$event)
db_escat <- db_escat %>%
  arrange(date_rand)
PFS_surv_obj_escat <- Surv(time = db_escat$time_event, event = db_escat$event)
summary(coxph(PFS_surv_obj_escat ~ armAB + strata(breast_molec_alt) + 
                strata(chemo_line) + strata(disease_status), data = db_escat))
# Stratified Cox models without AZD8931, AZD4547, AZD2014, and VANDETANIB
db_ext_escat <- db_escat[!(db_escat$targeted_therapy %in% c("AZD8931","AZD4547","AZD2014","VANDETANIB")
                           & db_escat$armAB=='A'),]
table(db_ext_escat$armAB,db_ext_escat$event)
db_ext_escat <- db_ext_escat %>%
  arrange(date_rand)
PFS_surv_obj_ext_escat  <- Surv(time = db_ext_escat$time_event, event = db_ext_escat$event)
summary(coxph(PFS_surv_obj_ext_escat ~ armAB + strata(breast_molec_alt) + 
                strata(chemo_line) + strata(disease_status), data = db_ext_escat))

# Kaplan-Meyer curve (with AZD8931, AZD4547, AZD2014, and VANDETANIB)
PFS_km_fit_escat_arm <- survfit(PFS_surv_obj_escat ~ armAB, data = db_escat, conf.int = 0.90)
PFS_plot_escat <- ggsurvplot(PFS_km_fit_escat_arm,
                             risk.table = TRUE,
                             censored = TRUE,
                             risk.table.col = "strata",
                             palette = c("blue", "red"),
                             xlab = "PFS in months",
                             xlim = c(0, 60),
                             break.time.by = 6,
                             ylab = "Survival probability",
                             legend.title = "",
                             legend.labs = c("Targeted therapy","Standard of care"))
print(PFS_plot_escat)

# Permutation-based test

# Stratified log-rank test
S_obs_strat_escat <- survdiff(PFS_surv_obj_ext_escat ~ armAB + strata(breast_molec_alt) + 
                                strata(chemo_line) + strata(disease_status), data = db_ext_escat)$chisq
pvalue_obs_strat_escat <- survdiff(PFS_surv_obj_ext_escat ~ armAB + strata(breast_molec_alt) + 
                                     strata(chemo_line) + strata(disease_status), data = db_ext_escat)$pvalue

# Parallelization
future::availableCores() # check the number of available cores 
plan(multisession, workers = future::availableCores() - 2)

L <- 15000
cov_list <- c("breast_molec_alt","chemo_line","disease_status")
ntrt <- length(unique(db_ext_escat$armAB))
ratio <- c(2, 1) # treated:control
covwt <- rep(1, length(cov_list))
rslt1 <- furrr::future_map(1:L, ~{
  simulate_one_rand_EXT(data = db_ext_escat, 
                        cov_list = cov_list, 
                        ntrt = ntrt, 
                        ratio = ratio, 
                        covwt = covwt
                        ) %>%
    mutate(sim=.x)
}, .progress = TRUE, .options = furrr_options(seed = 1234))
rslt1_df <- list_rbind(rslt1) 

print(paste("Observed p-value: ",pvalue_obs_strat_escat))
pvalue_rand_strat_escat <- 1/L*sum(ifelse(abs(rslt1_df$S_rand_strat) >= abs(S_obs_strat_escat),1,0))
print(paste("Permutation p-value: ",pvalue_rand_strat_escat))


## OVERALL population

# EXTERNAL MODIFICATIONS:
# drop of: AZD8931, AZD4547, AZD2014, and VANDETANIB
syn_BREAST %>%
  filter(armAB == "A") %>%
  with(table(targeted_therapy))

# Stratified Cox models with AZD8931, AZD4547, AZD2014, and VANDETANIB
table(syn_BREAST$armAB,syn_BREAST$event)
syn_BREAST <- syn_BREAST %>%
  arrange(date_rand)
PFS_surv_obj_BREAST <- Surv(time = syn_BREAST$time_event, event = syn_BREAST$event)
summary(coxph(PFS_surv_obj_BREAST ~ armAB + strata(breast_molec_alt) + 
                strata(chemo_line) + strata(disease_status), data = syn_BREAST))
# Stratified Cox models without AZD8931, AZD4547, AZD2014, and VANDETANIB
db_ext_BREAST <- syn_BREAST[!(syn_BREAST$targeted_therapy %in% c("AZD8931","AZD4547","AZD2014","VANDETANIB") 
                              & syn_BREAST$armAB=='A'),]
table(db_ext_BREAST$armAB,db_ext_BREAST$event)
db_ext_BREAST <- db_ext_BREAST %>%
  arrange(date_rand)
PFS_surv_obj_ext_BREAST  <- Surv(time = db_ext_BREAST$time_event, event = db_ext_BREAST$event)
summary(coxph(PFS_surv_obj_ext_BREAST ~ armAB + strata(breast_molec_alt) + 
                strata(chemo_line) + strata(disease_status), data = db_ext_BREAST))

# Kaplan-Meyer curve (with AZD8931, AZD4547, AZD2014, and VANDETANIB)
PFS_km_fit_BREAST_arm <- survfit(PFS_surv_obj_BREAST ~ armAB, data = syn_BREAST, conf.int = 0.95)
PFS_plot_BREAST <- ggsurvplot(PFS_km_fit_BREAST_arm,
                              risk.table = TRUE,
                              censored = TRUE,
                              risk.table.col = "strata",
                              palette = c("blue", "red"),
                              xlab = "PFS in months",
                              xlim = c(0, 60),
                              break.time.by = 6,
                              ylab = "Survival probability",
                              legend.title = "",
                              legend.labs = c("Targeted therapy","Standard of care"))
print(PFS_plot_BREAST)

# Permutation test

# stratified Cox models: stratified log-rank test
S_obs_strat_BREAST <- survdiff(PFS_surv_obj_ext_BREAST ~ armAB + strata(breast_molec_alt) + 
                                 strata(chemo_line) + strata(disease_status), data = db_ext_BREAST)$chisq
pvalue_obs_strat_BREAST <- survdiff(PFS_surv_obj_ext_BREAST ~ armAB + strata(breast_molec_alt) + 
                                      strata(chemo_line) + strata(disease_status), data = db_ext_BREAST)$pvalue

# Parallelization
# L <- 15000
cov_list <- c("breast_molec_alt","chemo_line","disease_status")
ntrt <- length(unique(db_ext_BREAST$armAB))
ratio <- c(2, 1) # treated:control
covwt <- rep(1, length(cov_list))
rslt2 <- furrr::future_map(1:L, ~{
  simulate_one_rand_EXT(data = db_ext_BREAST, 
                        cov_list = cov_list, 
                        ntrt = ntrt, 
                        ratio = ratio, 
                        covwt = covwt
                        ) %>%
    mutate(sim=.x)
}, .progress = TRUE, .options = furrr_options(seed = 1234))
rslt2_df <- list_rbind(rslt2) 

print(paste("Observed p-value: ",pvalue_obs_strat_BREAST))
pvalue_rand_strat_BREAST <- 1/L*sum(ifelse(abs(rslt2_df$S_rand_strat) >= abs(S_obs_strat_BREAST),1,0))
print(paste("Permutation p-value: ",pvalue_rand_strat_BREAST))


################################################################################


## Load the synthetic dataset generated from the SAFIR02-LUNG dataset.
## It includes the following variables:
## * "TARGET": name of the biomarker;
## * "treatment": therapy tailored to the characteristics of the patient;
## * "ARMCD_2cl": treatment arm ('Arm A: Targeted therapy' or 'Arm B: Standard maintenance therapy');
## * "histo_2cl": histologic classification;     
## * "Disease_status": disease status;
## * "Smoking_status": smoking status; 
## * "molecular_cat": molecular category;
## * "date_rand": date of randomization;
## * "PFS": follow-up;              
## * "PFS_event": event (1=event, 0=censored).   

syn_LUNG <- read_excel("C:/Users/C_MARELLI/OneDrive - gustaveroussy.fr/PhD PROJECT/papers/axis 1/github code/syn_LUNG.xlsx")
head(syn_LUNG)
nrow(syn_LUNG)

syn_LUNG$armAB <- as.factor(ifelse(syn_LUNG$ARMCD_2cl == "Arm B: Standard maintenance therapy","B","A"))
syn_LUNG$armAB <- relevel(syn_LUNG$armAB, ref = "B")
syn_LUNG$PFS_event <- as.numeric(syn_LUNG$PFS_event)


# INTERIM FUTILITY ANALYSIS:
# an interim futility analysis is performed after having observed 60 events 
# in the most frequent subgroup (KRAS) (irrespective of the treatment arm)
table(syn_LUNG$TARGET)
syn_LUNG %>%
  filter(TARGET == "KRAS") %>%
  with(table(PFS_event))
syn_LUNG %>%
  filter(TARGET == "KRAS") %>%
  with(table(armAB))

# date of event
syn_LUNG$date_rand <- as.Date(syn_LUNG$date_rand, format = "%d-%m-%Y", tz = "UTC")
syn_LUNG$PFS_date <- syn_LUNG$PFS + syn_LUNG$date_rand

# definition of the stage on the date of the 60th observed events in KRAS subgroup
syn_LUNG_KRAS_T <- syn_LUNG[syn_LUNG$TARGET %in% c("KRAS"),]
events <- syn_LUNG_KRAS_T %>%
  filter(PFS_event == 1) %>%
  arrange(PFS_date)
date_60 <- events$PFS_date[60]
date_60

syn_LUNG$stage <- rep(NA,nrow(syn_LUNG))
syn_LUNG$event <- rep(NA,nrow(syn_LUNG))
syn_LUNG$time_event <- rep(NA,nrow(syn_LUNG))
for(i in 1:nrow(syn_LUNG)){
  if(syn_LUNG$date_rand[i]<=date_60 & syn_LUNG$PFS_date[i]>date_60){
    syn_LUNG$stage[i] <- 1
    syn_LUNG$event[i] <- 0
    syn_LUNG$time_event[i] <- date_60
  }
  else if(syn_LUNG$date_rand[i]<=date_60 & syn_LUNG$PFS_date[i]<=date_60){
    syn_LUNG$stage[i] <- 1
    syn_LUNG$event[i] <- syn_LUNG$PFS_event[i]
    syn_LUNG$time_event[i] <- syn_LUNG$PFS[i]
  }
  else if(syn_LUNG$date_rand[i]>date_60){ # & syn_LUNG$PFS_date[i]>date_60
    syn_LUNG$stage[i] <- 2
    syn_LUNG$event[i] <- syn_LUNG$PFS_event[i]
    syn_LUNG$time_event[i] <- syn_LUNG$PFS[i]
  }
}
table(syn_LUNG$stage)
syn_LUNG %>%
  filter(TARGET == "KRAS",event == 1) %>%
  with(table(stage))

# Stratified Cox model with SELUMETINIB
syn_LUNG <- syn_LUNG %>%
  arrange(date_rand)
PFS_surv_obj_LUNG  <- Surv(time = syn_LUNG$time_event, event = syn_LUNG$event)
summary(coxph(PFS_surv_obj_LUNG ~ armAB + strata(histo_2cl) + strata(Disease_status) + 
                strata(Smoking_status) + strata(molecular_cat), data = syn_LUNG))
# Stratified Cox model without SELUMETINIB
syn_LUNG_noSELUMETINIB <- syn_LUNG[!(syn_LUNG$treatment %in% c("SELUMETINIB")),]
syn_LUNG_noSELUMETINIB <- syn_LUNG_noSELUMETINIB %>%
  arrange(date_rand)
PFS_surv_obj_noSELUMETINIB  <- Surv(time = syn_LUNG_noSELUMETINIB$time_event, 
                                    event = syn_LUNG_noSELUMETINIB$event)
summary(coxph(PFS_surv_obj_noSELUMETINIB ~ armAB + strata(histo_2cl) + strata(Disease_status) + 
                strata(Smoking_status) + strata(molecular_cat), data = syn_LUNG_noSELUMETINIB))

# Kaplan-Meyer curve (with SELUMETINIB)
PFS_km_fit_LUNG_arm <- survfit(PFS_surv_obj_LUNG ~ armAB, data = syn_LUNG)
PFS_plot_LUNG <- ggsurvplot(PFS_km_fit_LUNG_arm,
                            risk.table = TRUE,
                            censored = TRUE,
                            risk.table.col = "strata",
                            palette = c("blue", "red"),
                            xlab = "PFS in months",
                            # xlim = c(0, 12),
                            break.time.by = 3,
                            ylab = "Survival probability",
                            legend.title = "",
                            legend.labs = c("Targeted therapy","Standard of care"))
print(PFS_plot_LUNG)

# Permutation test

# stratified Cox models: stratified log-rank test
S_obs_strat_LUNG <- survdiff(PFS_surv_obj_LUNG ~ armAB + strata(histo_2cl) + strata(Disease_status) + 
                               strata(Smoking_status) + strata(molecular_cat), data = syn_LUNG)$chisq
pvalue_obs_strat_LUNG <- survdiff(PFS_surv_obj_LUNG ~ armAB + strata(histo_2cl) + strata(Disease_status) + 
                                    strata(Smoking_status) + strata(molecular_cat), data = syn_LUNG)$pvalue

# Parallelization
future::availableCores() # check the number of available cores 
plan(multisession, workers = future::availableCores() - 2)

target <- "KRAS"
cov_list <- c("histo_2cl","Disease_status","Smoking_status","molecular_cat")
ntrt <- length(unique(syn_LUNG$armAB))
ratio <- c(2, 1) # treated:control
covwt <- rep(1, length(cov_list))
crit_val <- 1.96
L <- 15000
max_attempts <- 100000
rslt3 <- furrr::future_map(1:1, ~{
  simulate_one_rand_IA(data=syn_LUNG,
                       target = target,
                       crit_val=crit_val, 
                       L=L, 
                       max_attempts=max_attempts,
                       cov_list = cov_list,
                       ntrt = ntrt,
                       ratio = ratio,
                       covwt = covwt
                       ) %>%
    mutate(sim=.x)
  
}, .progress = TRUE, .options = furrr_options(seed = 1234))
rslt3_df <- list_rbind(rslt3)

print(paste("Observed p-value: ",pvalue_obs_strat_LUNG))
pvalue_rand_strat_LUNG <- (1/rslt3_df$total_attempts)*
  sum(ifelse(abs(rslt3_df$S_rand_strat) >= abs(S_obs_strat_LUNG),1,0))
print(paste("Permutation p-value: ",pvalue_rand_strat_LUNG))
print(paste("Total number of attemps: ",rslt3_df$total_attempts))





