## This R code allows to replicate the two clinical trial applications presented in the paper:
## "Using permutation-based tests to evaluate an adaptive molecular treatment
## algorithm in randomized precision oncology trials" by Marelli et al.
## Before running the code, please download the synthetic datasets syn_BREAST.xlsx and syn_LUNG.xlsx.

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

## External event:
## First randomization of the treatment labels and randomization test statistic calculation:
## the minimization procedure is re-run and the unadjusted or adjusted score test
## statistic is re-computed.
## Arguments:
## * 'data' is the full dataset (including patients treated with the discontinued
## treatments at stage I);
## * 'keepID' is the list of patients ID retained;
## * 'cov_list' is the list of variables names to consider as minimization factors;
## * 'test' is the type of test to be performed (unadjusted or adjusted);
## * all the other entries follow the 'Minirand' function from the package 'Minirand'.
simulate_one_rand_EXT <- function(data, keepID, cov_list, ntrt, ratio, trtseq, prob_trt, covwt, p, test = "unadj") {
  
  data <- data %>%
    arrange(date_rand)
  
  N <- nrow(data)
  
  covmat <- as.matrix(data[, cov_list])
  
  res <- rep(NA, N)
  res[1] <- sample(trtseq, prob=prob_trt, size=1)
  
  for(j in 2:N) {
    res[j] <- Minirand(
      covmat = covmat,
      j = j,
      covwt = covwt,
      ratio = ratio,
      ntrt = ntrt,
      trtseq = trtseq,
      method = "Range",
      result = res,
      p = p
    )
  }
  
  data_temp <- data
  data_temp$trt_rand <- res
  data_analysis <- data_temp %>% filter(ID %in% keepID)
  data_temp$armAB_rand <- ifelse(data_temp$trt_rand == 1,"A","B")
  data_temp$armAB_rand <- as.factor(data_temp$armAB_rand)
  data_temp$armAB_rand <- relevel(data_temp$armAB_rand, ref = "B")
  data_temp$armAB_01_rand <- ifelse(data_temp$armAB_rand == "A",1,0)
  
  data_temp <- data_temp %>%
    arrange(date_rand)
  
  if(test=="unadj"){
    S_rand <- survdiff(Surv(time_event,event) ~ armAB_01_rand, data = data_temp)$chisq
  }
  else if(test=="adj"){
    data_temp$NEW_breast_molec_alt <- ifelse(data_temp$breast_molec_alt=="C","A",data_temp$breast_molec_alt)
    data_temp$NEW_breast_molec_alt_f <- factor(data_temp$NEW_breast_molec_alt)
    # data_temp$breast_molec_alt_f      <- factor(data_temp$breast_molec_alt)
    data_temp$chemo_line_f            <- factor(data_temp$chemo_line)
    data_temp$disease_status_f        <- factor(data_temp$disease_status)
    data_temp$NEW_breast_molec_alt_f_drop <- droplevels(data_temp$NEW_breast_molec_alt_f)
    
    data_temp <- data_temp %>%
      arrange(date_rand)
    
    fit <- tryCatch(
      coxph(Surv(time_event, event) ~ NEW_breast_molec_alt_f_drop +
              armAB_01_rand:NEW_breast_molec_alt_f_drop + chemo_line_f + disease_status_f,
            data = data_temp),
      error = function(e) NULL)
    
    S_rand <- NA_real_
    if (!is.null(fit)) {
      cf     <- coef(fit)
      int_nm <- grep("armAB_01_rand", names(cf), value = TRUE)
      int_nm <- int_nm[!is.na(cf[int_nm])]
      lt     <- sub(".*_drop([A-Za-z0-9]+)(:.*|$)", "\\1", int_nm)   # es. "A" "B" "D"
      
      if (length(lt) >= 2) {
        mf <- model.frame(fit)
        weights  <- table(mf$NEW_breast_molec_alt_f_drop)[lt]
        weights  <- as.numeric(weights / sum(weights))
        names(weights) <- lt
        
        ref    <- lt[length(lt)]
        others <- lt[-length(lt)]
        
        for (l in lt) {
          data_temp[[paste0("int", l)]] <-
            as.numeric(data_temp$NEW_breast_molec_alt_f_drop == l) * data_temp$armAB_01_rand
        }
        for (o in others) {
          data_temp[[paste0("Z", o)]] <-
            data_temp[[paste0("int", o)]] - (weights[o] / weights[ref]) * data_temp[[paste0("int", ref)]]
        }
        
        # int_cols <- paste0("int", lt)
        z_cols   <- paste0("Z", others)
        int_ref <- paste0("int", ref)
        
        f_red  <- as.formula(paste("Surv(time_event, event) ~ NEW_breast_molec_alt_f_drop +",
                                   paste(c(z_cols,   "chemo_line_f", "disease_status_f"), collapse = " + ")))
        f_full <- as.formula(paste("Surv(time_event, event) ~ NEW_breast_molec_alt_f_drop +",
                                   paste(c(z_cols, "chemo_line_f", "disease_status_f", int_ref), collapse = " + ")))
        
        S_rand <- tryCatch({
          m_red  <- coxph(f_red,  data = data_temp)
          m_full <- coxph(f_full, data = data_temp, init = c(coef(m_red), 0), control = coxph.control(iter.max = 0))
          if(any(is.null(coef(m_red)))) print("NA in reduced model")
          if(any(is.null(coef(m_full)))) print("NA in full model")
          # anova(m_red, m_full)$Chisq[2]
          unname(m_full$score)
        }, error = function(e) NA_real_)
      }
    }
  }
  
  
  tibble::tibble(S_rand) 
}

## Interim analysis:
## First randomization of the treatment labels and randomization test statistic calculation:
## the minimization procedure is re-run on full stage I data and the Z-score for
## the log-hazard ratio is compared to an interim futility threshold. If the trial
## continues to stage II, the minimization procedure is re-run on stage II data and
## the adjusted score test statistic is re-computed.
## Arguments:
## * 'data' is the full dataset (including patients treated with the discontinued
## treatments at stage I);
## * 'target' is the name of the molecular alteration on which to perform the
## interim futility analysis;
## * 'cov_list' is the list of variables names to consider as minimization factors;
## * 'crit_val' is the interim futility threshold;
## * 'L' is the planned number of valid randomizations;
## * 'max_attempts' is the maximum number of attempts to run to find the 'L'
## required valid randomizations;
## * all the other entries follow the 'Minirand' function from the package 'Minirand'.
simulate_one_rand_IA <- function(data, target, cov_list, crit_val=0, L=15000, 
                                 max_attempts=100000, ntrt, trtseq, prob_trt, ratio, covwt, p) {
  
  S_rand_L <- rep(NA_real_,L)
  l <- 1  # valid randomization counter
  attempts <- 0 # number of attempts counter
  
  while(l <= L){
    attempts <- attempts + 1
    
    if(attempts > max_attempts){
      warning(paste("Only ", l-1, " randomizations obtained out of ", L, " randomizations"))
      S_rand_L <- S_rand_L[1:(l-1)] # S_rand_L[seq_len(l-1)]
      break
    }
    
    # STAGE I
    
    data_stageI <- data[data$stage==1,]
    data_stageI <- data_stageI %>%
      arrange(date_rand)
    
    N1 <- nrow(data_stageI)
    
    covmat1 <- as.matrix(data_stageI[, cov_list])
    
    res1 <- numeric(N1)
    res1[1] <- sample(trtseq, prob=prob_trt, size=1)
    
    for(j in 2:N1) {
      res1[j] <- Minirand(
        covmat = covmat1,
        j = j,
        covwt = covwt,
        ratio = ratio,
        ntrt = ntrt,
        trtseq = trtseq,
        method = "Range",
        result = res1,
        p = p
      )
    }
    
    new_data_stageI <- data_stageI
    new_data_stageI <- new_data_stageI %>%
      arrange(date_rand)
    new_data_stageI$trt_rand <- res1
    new_data_stageI$armAB_rand <- ifelse(new_data_stageI$trt_rand == 1, "A", "B")
    new_data_stageI$armAB_rand <- as.factor(new_data_stageI$armAB_rand)
    new_data_stageI$armAB_rand <- relevel(new_data_stageI$armAB_rand, ref = "B")
    new_data_stageI$armAB_01_rand <- ifelse(new_data_stageI$armAB_rand == "A",1,0)
    
    data_interim <- new_data_stageI[new_data_stageI$TARGET==target,]
    if (nrow(data_interim) == 0)
      stop(paste("No rows found for TARGET =", target))
    
    data_interim <- data_interim %>%
      arrange(date_rand)
    
    cox_interim <- coxph(Surv(time_event,event) ~ armAB_01_rand + strata(histo_2cl_f) + 
                           strata(Disease_status_f) + strata(Smoking_status_f), data = data_interim)
    
    logHR_interim <- coef(cox_interim)        
    SE_logHR_interim <- sqrt(diag(vcov(cox_interim))) 
    
    Z_score_interim <- logHR_interim / SE_logHR_interim
    
    if(Z_score_interim<crit_val){
      
      data_stageII <- data[data$stage==2,]
      data_stageII <- data_stageII %>%
        arrange(date_rand)
      
      N2 <- nrow(data_stageII)
      
      covmat2 <- as.matrix(data_stageII[, cov_list])
      covmat_tot <- as.matrix(rbind(covmat1,covmat2))
      
      res_tot <- c(res1,numeric(N2))
      for(j in (N1+1):(N1+N2)) {
        res_tot[j] <- Minirand(
          covmat = covmat_tot,
          j = j,
          covwt = covwt,
          ratio = ratio,
          ntrt = ntrt,
          trtseq = trtseq,
          method = "Range",
          result = res_tot,
          p = p
        )
      }
      res2 <- res_tot[(N1 + 1):(N1 + N2)]
      
      new_data_stageII <- data_stageII
      new_data_stageII <- new_data_stageII %>%
        arrange(date_rand)
      new_data_stageII$trt_rand <- res2
      new_data_stageII$armAB_rand <- ifelse(new_data_stageII$trt_rand == 1, "A", "B")
      new_data_stageII$armAB_rand <- as.factor(new_data_stageII$armAB_rand)
      new_data_stageII$armAB_rand <- relevel(new_data_stageII$armAB_rand, ref = "B")
      new_data_stageII$armAB_01_rand <- ifelse(new_data_stageII$armAB_rand == "A",1,0)
      
      data_rand <- rbind(new_data_stageI,new_data_stageII)
      
      data_rand <- data_rand %>%
        arrange(stage,date_rand)
      
      data_rand$molecular_cat_f_drop <- droplevels(data_rand$molecular_cat_f)
      
      fit <- tryCatch(
        coxph(Surv(time_event, event) ~ molecular_cat_f_drop +
                armAB_01_rand:molecular_cat_f_drop + histo_2cl_f + Disease_status_f + Smoking_status_f,
              data = data_rand),
        error = function(e) NULL)
      
      if (!is.null(fit)) {
        cf     <- coef(fit)
        int_nm <- grep("armAB_01_rand", names(cf), value = TRUE)
        int_nm <- int_nm[!is.na(cf[int_nm])]
        lt     <- sub(".*molecular_cat_f_drop([^:]+).*", "\\1", int_nm)   # ex: Category A, , Category B, Category C, Category D
        
        if (length(lt) >= 2) {
          mf <- model.frame(fit)
          weights  <- table(mf$molecular_cat_f_drop)[lt]
          weights  <- as.numeric(weights / sum(weights))
          names(weights) <- lt
          
          ref    <- lt[length(lt)]
          others <- lt[-length(lt)]
          
          int_cols <- paste0("int", seq_len(length(lt)))
          z_cols   <- paste0("Z",   seq_len(length(lt) - 1))
          int_ref <- paste0("int", length(lt))
          
          for (i in seq_len(length(lt))) {
            data_rand[[int_cols[i]]] <-
              as.numeric(data_rand$molecular_cat_f_drop == lt[i]) * data_rand$armAB_01_rand
          }
          for (i in seq_len(length(lt) - 1)) {
            data_rand[[z_cols[i]]] <-
              data_rand[[int_cols[i]]] -
              (weights[lt[i]] / weights[ref]) * data_rand[[int_cols[length(lt)]]]
          }
          
          f_red  <- as.formula(paste("Surv(time_event, event) ~ molecular_cat_f_drop +",
                                     paste(c(z_cols,   "histo_2cl_f", "Disease_status_f", "Smoking_status_f"), collapse = " + ")))
          f_full <- as.formula(paste("Surv(time_event, event) ~ molecular_cat_f_drop +",
                                     paste(c(z_cols, "histo_2cl_f", "Disease_status_f", "Smoking_status_f", int_ref), collapse = " + ")))
          
          S_rand_L[l] <- tryCatch({
            m_red  <- coxph(f_red,  data = data_rand)
            m_full <- coxph(f_full, data = data_rand, init = c(coef(m_red), 0), control = coxph.control(iter.max = 0))
            if(any(is.null(coef(m_red)))) print("NA in reduced model")
            if(any(is.null(coef(m_full)))) print("NA in full model")
            # anova(m_red, m_full)$Chisq[2]
            unname(m_full$score)
          }, error = function(e) NA_real_)
        }
      }
      
      l <- l+1
    }
  }
  
  tibble::tibble(S_rand_L,total_attempts = attempts)
}


################################################################################


## Load the synthetic dataset generated from the SAFIR02-BREAST dataset (syn_BREAST.xlsx).
## It includes the following variables:
## * "ID": patient's ID
## * "TARGET": name of the biomarker;
## * "targeted_therapy": therapy tailored to the characteristics of the patient;
## * "armAB": treatment arm ('A' or 'B');
## * "breast_molec_alt": breast molecular alteration;
## * "chemo_line": line of chemotherapy at the time of randomization ;      
## * "disease_status": disease status at randomization;
## * "escat_class": escat classification;
## * "date_rand": date of randomization;
## * "PFS": follow-up;              
## * "PFS_event": event (1=event, 0=censored).      

syn_BREAST <- read_excel("syn_BREAST.xlsx") # change path
head(syn_BREAST)
nrow(syn_BREAST)

syn_BREAST$armAB <- as.factor(syn_BREAST$armAB)
syn_BREAST$armAB <- relevel(syn_BREAST$armAB, ref = "B")
syn_BREAST$armAB_01 <- ifelse(syn_BREAST$armAB=="A",1,0)

syn_BREAST$breast_molec_alt_f <- factor(syn_BREAST$breast_molec_alt)
syn_BREAST$chemo_line_f <- factor(syn_BREAST$chemo_line)
syn_BREAST$disease_status_f <- factor(syn_BREAST$disease_status)

syn_BREAST$event <- as.numeric(syn_BREAST$PFS_event)
syn_BREAST$time_event <- as.numeric(syn_BREAST$PFS)

# EXTERNAL MODIFICATIONS:
# drop of: AZD8931, AZD4547, AZD2014, and VANDETANIB
syn_BREAST %>%
  filter(armAB == "A") %>%
  with(table(targeted_therapy))

targets_of_interest <- syn_BREAST %>%
  filter(
    armAB == "A",
    targeted_therapy %in% c("AZD8931", "AZD4547", "AZD2014", "VANDETANIB")
  ) %>%
  pull(TARGET) %>%
  unique()  
print(targets_of_interest)

## ESCAT I/II population
db_escat <- syn_BREAST[syn_BREAST$escat_class %in% c("I","II"),]
nrow(db_escat)

# Unadjusted Cox models without AZD8931, AZD4547, AZD2014, and VANDETANIB
db_ext_escat <- db_escat[!(db_escat$TARGET %in% targets_of_interest),]
escatID <- db_ext_escat$ID
length(escatID)
table(db_ext_escat$armAB,db_ext_escat$event)
db_ext_escat <- db_ext_escat %>%
  arrange(date_rand)
summary(coxph(Surv(time_event, event) ~ armAB_01, data = db_ext_escat), conf.int = 0.90)

# Kaplan-Meier curve (without AZD8931, AZD4547, AZD2014, and VANDETANIB)
PFS_km_fit_ext_escat_arm <- survfit(Surv(time_event, event) ~ armAB_01, data = db_ext_escat, conf.int = 0.90)
PFS_plot_ext_escat <- ggsurvplot(PFS_km_fit_ext_escat_arm,
                                  risk.table = TRUE,
                                  censored = TRUE,
                                  risk.table.col = "strata",
                                  palette = c("red", "blue"),
                                  xlab = "PFS in months",
                                  xlim = c(0, 60),
                                  break.time.by = 6,
                                  ylab = "Survival probability",
                                  legend.title = "",
                                  legend.labs = c("Standard of care","Targeted therapy"))
print(PFS_plot_ext_escat)


# Randomization test

# Unstratified score (or log-rank) test
model_unadj_escat <- survdiff(Surv(time_event, event) ~ armAB_01, data = db_ext_escat)
S_obs_unadj_escat <- model_unadj_escat$chisq
pvalue_obs_unadj_escat <- model_unadj_escat$pvalue

# Parallelization
# future::availableCores() # check the number of available cores 
plan(multisession, workers = future::availableCores() - 2)

L <- 15000
cov_list <- c("breast_molec_alt","chemo_line","disease_status")
ntrt <- length(unique(db_ext_escat$armAB))
ratio <- c(2, 1) # treated:control
covwt <- rep(1, length(cov_list))
rslt1 <- furrr::future_map(1:L, ~{
  simulate_one_rand_EXT(data = db_escat, 
                        keepID = escatID,
                        cov_list = cov_list, 
                        ntrt = ntrt, 
                        trtseq = c(1,0),
                        prob_trt = c(0.5,0.5),
                        ratio = ratio, 
                        covwt = covwt,
                        test = "unadj",
                        p = 0.9
  ) %>%
    mutate(sim=.x)
}, .progress = TRUE, .options = furrr_options(seed = 1234))
rslt1_df <- list_rbind(rslt1) 

print(paste("Observed p-value: ",pvalue_obs_unadj_escat))
pvalue_rand_unadj_escat <- 1/L*sum(ifelse(abs(rslt1_df$S_rand) >= abs(S_obs_unadj_escat),1,0))
print(paste("Permutation p-value: ",pvalue_rand_unadj_escat))


## OVERALL population
db_ext_BREAST <- syn_BREAST[!(syn_BREAST$TARGET %in% targets_of_interest),]
extID <- db_ext_BREAST$ID
length(extID)
table(db_ext_BREAST$armAB,db_ext_BREAST$event)
db_ext_BREAST <- db_ext_BREAST %>%
  arrange(date_rand)

# Adjusted Cox models without AZD8931, AZD4547, AZD2014, and VANDETANIB
db_ext_BREAST$breast_molec_alt_f_drop <- droplevels(db_ext_BREAST$breast_molec_alt_f)
summary(coxph(Surv(time_event, event) ~ breast_molec_alt_f_drop + armAB_01:breast_molec_alt_f_drop + chemo_line_f + disease_status_f, data = db_ext_BREAST))

db_ext_BREAST$NEW_breast_molec_alt <- ifelse(db_ext_BREAST$breast_molec_alt=="C","A",db_ext_BREAST$breast_molec_alt)
db_ext_BREAST$NEW_breast_molec_alt_f <- factor(db_ext_BREAST$NEW_breast_molec_alt)
table(db_ext_BREAST$NEW_breast_molec_alt_f)
db_ext_BREAST$NEW_breast_molec_alt_f_drop <- droplevels(db_ext_BREAST$NEW_breast_molec_alt_f)
fit <- coxph(Surv(time_event, event) ~ NEW_breast_molec_alt_f_drop + armAB_01:NEW_breast_molec_alt_f_drop + chemo_line_f + disease_status_f,data = db_ext_BREAST)
summary(fit)
cf     <- coef(fit)
int_nm <- grep("armAB_01", names(cf), value = TRUE)     
ok     <- !is.na(cf[int_nm])
if (any(!ok)) message("Gruppi esclusi (coef NA): ", paste(int_nm[!ok], collapse = ", "))
int_nm <- int_nm[ok]
lt     <- sub(".*_drop([A-Za-z0-9]+)(:.*|$)", "\\1", int_nm)   # ex: "A" "B" "D"
names(int_nm) <- lt
mf <- model.frame(fit)
weights  <- table(mf$NEW_breast_molec_alt_f_drop)[lt]
weights  <- as.numeric(weights / sum(weights))
names(weights) <- lt
b   <- cf[int_nm]
V   <- vcov(fit)[int_nm, int_nm, drop = FALSE]
est <- sum(weights * b)
se  <- sqrt(drop(t(weights) %*% V %*% weights))
res <- c(HR    = exp(est),
         lower = exp(est - 1.96 * se),
         upper = exp(est + 1.96 * se),
         p     = 2 * pnorm(-abs(est / se)))
res

# Kaplan-Meier curve (without AZD8931, AZD4547, AZD2014, and VANDETANIB)
PFS_km_fit_ext_BREAST_arm <- survfit(Surv(time_event, event) ~ armAB_01, data = db_ext_BREAST, conf.int = 0.95)
PFS_plot_ext_BREAST <- ggsurvplot(PFS_km_fit_ext_BREAST_arm,
                              risk.table = TRUE,
                              censored = TRUE,
                              risk.table.col = "strata",
                              palette = c("red", "blue"),
                              xlab = "PFS in months",
                              xlim = c(0, 60),
                              break.time.by = 6,
                              ylab = "Survival probability",
                              legend.title = "",
                              legend.labs = c("Standard of care","Targeted therapy"))
print(PFS_plot_ext_BREAST)


# Randomization test

# Adjusted score test
ref    <- lt[length(lt)]       # reference: last estimable level (ex: D)
others <- lt[-length(lt)]

for (l in lt) {
  db_ext_BREAST[[paste0("int", l)]] <-
    as.numeric(db_ext_BREAST$NEW_breast_molec_alt_f_drop == l) * db_ext_BREAST$armAB_01
}
for (o in others) {
  db_ext_BREAST[[paste0("Z", o)]] <-
    db_ext_BREAST[[paste0("int", o)]] - (weights[o] / weights[ref]) * db_ext_BREAST[[paste0("int", ref)]]
}

# int_cols <- paste0("int", lt)
z_cols   <- paste0("Z", others)
int_ref <- paste0("int", ref)

f_red  <- as.formula(paste("Surv(time_event, event) ~ NEW_breast_molec_alt_f_drop +",
                           paste(c(z_cols,   "chemo_line_f", "disease_status_f"), collapse = " + ")))
f_full <- as.formula(paste("Surv(time_event, event) ~ NEW_breast_molec_alt_f_drop +",
                           paste(c(z_cols, "chemo_line_f", "disease_status_f", int_ref), collapse = " + ")))

model_red_BREAST  <- coxph(f_red,  data = db_ext_BREAST)
model_full_BREAST <- coxph(f_full, data = db_ext_BREAST, init = c(coef(model_red_BREAST), 0), control = coxph.control(iter.max = 0))

S_obs_adj_BREAST <- unname(model_full_BREAST$score)
pvalue_obs_adj_BREAST <- pchisq(S_obs_adj_BREAST, df = 1, lower.tail = FALSE)

# Parallelization
# future::availableCores() # check the number of available cores 
# plan(multisession, workers = future::availableCores() - 2)

# L <- 15000
cov_list <- c("breast_molec_alt","chemo_line","disease_status")
ntrt <- length(unique(db_ext_BREAST$armAB))
ratio <- c(2, 1) # treated:control
covwt <- rep(1, length(cov_list))
rslt2 <- furrr::future_map(1:L, ~{
  simulate_one_rand_EXT(data = db_ext_BREAST, 
                        keepID = extID,
                        cov_list = cov_list, 
                        ntrt = ntrt, 
                        trtseq = c(1,0),
                        prob_trt = c(0.5,0.5),
                        ratio = ratio, 
                        covwt = covwt,
                        test = "adj",
                        p = 0.9
                        ) %>%
    mutate(sim=.x)
}, .progress = TRUE, .options = furrr_options(seed = 1234))
rslt2_df <- list_rbind(rslt2) 

print(paste("Observed p-value: ",pvalue_obs_adj_BREAST))
pvalue_rand_adj_BREAST <- 1/L*sum(ifelse(abs(rslt2_df$S_rand) >= abs(S_obs_adj_BREAST),1,0))
print(paste("Permutation p-value: ",pvalue_rand_adj_BREAST))


################################################################################


## Load the synthetic dataset generated from the SAFIR02-LUNG dataset (syn_LUNG.xlsx).
## It includes the following variables:
## * "ID": patient's ID;
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

syn_LUNG <- read_excel("syn_LUNG.xlsx") # change path
head(syn_LUNG)
nrow(syn_LUNG)

syn_LUNG$armAB <- as.factor(ifelse(syn_LUNG$ARMCD_2cl == "Arm B: Standard maintenance therapy","B","A"))
syn_LUNG$armAB <- relevel(syn_LUNG$armAB, ref = "B")
syn_LUNG$armAB_01 <- ifelse(syn_LUNG$armAB=="A",1,0)

syn_LUNG$histo_2cl_f <- factor(syn_LUNG$histo_2cl)
syn_LUNG$Disease_status_f <- factor(syn_LUNG$Disease_status)
syn_LUNG$Smoking_status_f <- factor(syn_LUNG$Smoking_status)
syn_LUNG$molecular_cat_f <- factor(syn_LUNG$molecular_cat)

syn_LUNG$date_rand <- as.Date(syn_LUNG$date_rand, format = "%d-%m-%Y", tz = "UTC")
# syn_LUNG$event_date <- syn_LUNG$PFS + syn_LUNG$date_rand
syn_LUNG$event_date <- syn_LUNG$date_rand + round(syn_LUNG$PFS * 30.5)

syn_LUNG$PFS_event_num <- as.numeric(syn_LUNG$PFS_event)


# INTERIM FUTILITY ANALYSIS:
# an interim futility analysis is performed after having observed 60 events 
# in the most frequent subgroup (KRAS) (irrespective of the treatment arm)
table(syn_LUNG$TARGET)
syn_LUNG %>%
  filter(TARGET == "KRAS") %>%
  with(table(PFS_event_num))
syn_LUNG %>%
  filter(TARGET == "KRAS") %>%
  with(table(armAB))

# Unajusted Cox model without SELUMETINIB
syn_LUNG_noSELUMETINIB <- syn_LUNG[!(syn_LUNG$treatment %in% c("SELUMETINIB")),]
syn_LUNG_noSELUMETINIB <- syn_LUNG_noSELUMETINIB %>%
  arrange(date_rand)
summary(coxph(Surv(PFS, PFS_event_num) ~ armAB_01, data = syn_LUNG_noSELUMETINIB))

# KM curves for PFS without SELUMETINIB
PFS_km_fit_arm_noSELUMETINIB <- survfit(Surv(PFS,PFS_event_num) ~ armAB_01, data = syn_LUNG_noSELUMETINIB)
PFS_plot_noSELUMETINIB <- ggsurvplot(PFS_km_fit_arm_noSELUMETINIB,
                                     risk.table = TRUE,
                                     censored = TRUE,
                                     risk.table.col = "strata",
                                     palette = c("red", "blue"),
                                     xlab = "PFS in months",
                                     xlim = c(0, 42),
                                     break.time.by = 6,
                                     ylab = "Survival probability",
                                     legend.title = "",
                                     legend.labs = c("Standard of care","Targeted therapy"))
print(PFS_plot_noSELUMETINIB)

# definition of the stage on the date of the 60th observed events in KRAS subgroup
syn_LUNG_KRAS_T <- syn_LUNG[syn_LUNG$TARGET %in% c("KRAS"),]
events <- syn_LUNG_KRAS_T %>%
  filter(PFS_event_num == 1) %>%
  arrange(event_date)
date_60 <- events$event_date[60]
date_60

syn_LUNG$stage <- rep(NA,nrow(syn_LUNG))
syn_LUNG$event <- rep(NA,nrow(syn_LUNG))
syn_LUNG$time_event <- rep(NA,nrow(syn_LUNG))
for(i in 1:nrow(syn_LUNG)){
  if(syn_LUNG$date_rand[i]<=date_60 & syn_LUNG$event_date[i]>date_60){
    syn_LUNG$stage[i] <- 1
    syn_LUNG$event[i] <- 0
    syn_LUNG$time_event[i] <- as.numeric(date_60 - syn_LUNG$date_rand[i]) / 30.5
  }
  else if(syn_LUNG$date_rand[i]<=date_60 & syn_LUNG$event_date[i]<=date_60){
    syn_LUNG$stage[i] <- 1
    syn_LUNG$event[i] <- syn_LUNG$PFS_event_num[i]
    syn_LUNG$time_event[i] <- syn_LUNG$PFS[i]
  }
  else if(syn_LUNG$date_rand[i]>date_60){ # & syn_LUNG$event_date[i]>date_60
    syn_LUNG$stage[i] <- 2
    syn_LUNG$event[i] <- syn_LUNG$PFS_event_num[i]
    syn_LUNG$time_event[i] <- syn_LUNG$PFS[i]
  }
}
table(syn_LUNG$stage)
syn_LUNG %>%
  filter(TARGET == "KRAS",event == 1) %>%
  with(table(stage))

# futility decision on observed data
syn_LUNG_first60 <- syn_LUNG[syn_LUNG$stage==1 & syn_LUNG$TARGET == "KRAS",]
table(syn_LUNG_first60$event)
table(syn_LUNG_first60$event,syn_LUNG_first60$armAB_01)
table(syn_LUNG_first60$event,syn_LUNG_first60$treatment)

syn_LUNG_first60 <- syn_LUNG_first60 %>%
  arrange(stage,date_rand)
# summary(coxph(Surv(time_event,event) ~ armAB_01, data = syn_LUNG_first60))
cox_model <- coxph(Surv(time_event,event) ~ armAB_01 + strata(histo_2cl_f) + strata(Disease_status_f) + 
                     strata(Smoking_status_f), data = syn_LUNG_first60)
summary(cox_model)

logHR <- coef(cox_model)                # log hazard ratio
SE_logHR <- sqrt(diag(vcov(cox_model))) # standard error

Z_score <- logHR / SE_logHR
print(Z_score) # -0.66 < 0,0.34


# Randomization test

# Adjusted score test
syn_LUNG$molecular_cat_f_drop <- droplevels(syn_LUNG$molecular_cat_f)
fit <- coxph(Surv(PFS,PFS_event_num) ~ molecular_cat_f_drop + armAB_01:molecular_cat_f_drop + 
               histo_2cl_f + Disease_status_f + Smoking_status_f, data = syn_LUNG)
summary(fit)
cf     <- coef(fit)
int_nm <- grep("armAB_01", names(cf), value = TRUE)     
ok     <- !is.na(cf[int_nm])
if (any(!ok)) message("Gruppi esclusi (coef NA): ", paste(int_nm[!ok], collapse = ", "))
int_nm <- int_nm[ok]
lt     <- sub(".*molecular_cat_f_drop([^:]+).*", "\\1", int_nm) # ex: Category A, Category B, Category C, Category D
names(int_nm) <- lt
mf <- model.frame(fit)
weights  <- table(mf$molecular_cat_f_drop)[lt]
weights  <- as.numeric(weights / sum(weights))
names(weights) <- lt
b   <- cf[int_nm]
V   <- vcov(fit)[int_nm, int_nm, drop = FALSE]
est <- sum(weights * b)
se  <- sqrt(drop(t(weights) %*% V %*% weights))
res <- c(HR    = exp(est),
         lower = exp(est - 1.96 * se),
         upper = exp(est + 1.96 * se),
         p     = 2 * pnorm(-abs(est / se)))
res


ref    <- lt[length(lt)]       # reference: last estimable level (ex: D)
others <- lt[-length(lt)]

int_cols <- paste0("int", seq_len(length(lt)))
z_cols   <- paste0("Z",   seq_len(length(lt) - 1))
int_ref <- paste0("int", length(lt))

for (i in seq_len(length(lt))) {
  syn_LUNG[[int_cols[i]]] <-
    as.numeric(syn_LUNG$molecular_cat_f_drop == lt[i]) * syn_LUNG$armAB_01
}
for (i in seq_len(length(lt) - 1)) {
  syn_LUNG[[z_cols[i]]] <-
    syn_LUNG[[int_cols[i]]] -
    (weights[lt[i]] / weights[ref]) * syn_LUNG[[int_cols[length(lt)]]]
}

f_red  <- as.formula(paste("Surv(time_event, event) ~ molecular_cat_f_drop +",
                           paste(c(z_cols,   "histo_2cl_f", "Disease_status_f", "Smoking_status_f"), collapse = " + ")))
f_full <- as.formula(paste("Surv(time_event, event) ~ molecular_cat_f_drop +",
                           paste(c(z_cols, "histo_2cl_f", "Disease_status_f", "Smoking_status_f", int_ref), collapse = " + ")))

model_red_LUNG  <- coxph(f_red,  data = syn_LUNG)
model_full_LUNG <- coxph(f_full, data = syn_LUNG, init = c(coef(model_red_LUNG), 0), control = coxph.control(iter.max = 0))

S_obs_adj_LUNG <- unname(model_full_LUNG$score)
pvalue_obs_adj_LUNG <- pchisq(S_obs_adj_LUNG, df = 1, lower.tail = FALSE)

# Parallelization
# future::availableCores() # check the number of available cores 
# plan(multisession, workers = future::availableCores() - 2)

# L <- 15000
target <- "KRAS"
cov_list <- c("histo_2cl","Disease_status","Smoking_status","molecular_cat")
ntrt <- length(unique(syn_LUNG$armAB))
ratio <- c(2, 1) # treated:control
covwt <- rep(1, length(cov_list))
crit_val <- 0 # 0.34
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
                       trtseq = c(1,0),
                       prob_trt = c(0.5,0.5), 
                       covwt = covwt,
                       p = 0.9
                       ) %>%
    mutate(sim=.x)
  
}, .progress = TRUE, .options = furrr_options(seed = 1234))
rslt3_df <- list_rbind(rslt3)

print(paste("Observed p-value: ",pvalue_obs_adj_LUNG))
print(paste("Total number of attemps: ",rslt3_df$total_attempts[1]))
pvalue_rand_adj_LUNG <- (1/rslt3_df$total_attempts[1])*sum(ifelse(abs(rslt3_df$S_rand_L) >= abs(S_obs_adj_LUNG),1,0))
print(paste("Permutation p-value: ",pvalue_rand_adj_LUNG))

simulate_drops <- function(data, target, cov_list, ntrt, trtseq, prob_trt, ratio, covwt, p) {
  
  # STAGE I
  
  data_stageI <- data[data$stage==1,]
  data_stageI <- data_stageI %>%
    arrange(date_rand)
  
  N1 <- nrow(data_stageI)
  
  covmat1 <- as.matrix(data_stageI[, cov_list])
  
  res1 <- numeric(N1)
  res1[1] <- sample(trtseq, prob=prob_trt, size=1)
  
  for(j in 2:N1) {
    res1[j] <- Minirand(
      covmat = covmat1,
      j = j,
      covwt = covwt,
      ratio = ratio,
      ntrt = ntrt,
      trtseq = trtseq,
      method = "Range",
      result = res1,
      p = p
    )
  }
  
  new_data_stageI <- data_stageI
  new_data_stageI <- new_data_stageI %>%
    arrange(date_rand)
  new_data_stageI$trt_rand <- res1
  new_data_stageI$armAB_rand <- ifelse(new_data_stageI$trt_rand == 1, "A", "B")
  new_data_stageI$armAB_rand <- as.factor(new_data_stageI$armAB_rand)
  new_data_stageI$armAB_rand <- relevel(new_data_stageI$armAB_rand, ref = "B")
  new_data_stageI$armAB_01_rand <- ifelse(new_data_stageI$armAB_rand == "A",1,0)
  
  data_interim <- new_data_stageI[new_data_stageI$TARGET==target,]
  if (nrow(data_interim) == 0)
    stop(paste("No rows found for TARGET =", target))
  
  data_interim <- data_interim %>%
    arrange(date_rand)
  
  cox_interim <- coxph(Surv(time_event,event) ~ armAB_01_rand + strata(histo_2cl_f) + 
                         strata(Disease_status_f) + strata(Smoking_status_f), data = data_interim)
  
  logHR_interim <- coef(cox_interim)        
  SE_logHR_interim <- sqrt(diag(vcov(cox_interim))) 
  
  Z_score_interim <- logHR_interim / SE_logHR_interim
  
  
  tibble::tibble(Z_score_interim)
}

rslt4 <- furrr::future_map(1:15000, ~{
  simulate_drops(data=syn_LUNG,
                 target = target,
                 cov_list = cov_list,
                 ntrt = ntrt,
                 trtseq = c(1,0),
                 prob_trt = c(0.5,0.5), 
                 ratio = ratio,
                 covwt = covwt,
                 p = 0.9
  ) %>%
    mutate(sim=.x)
  
}, .progress = TRUE, .options = furrr_options(seed = 234))
rslt4_df <- list_rbind(rslt4)

print(paste("Probability of taking the same interim decision (KRAS subgroup not dropped): ", table(rslt4_df$Z_score_interim>=0))) # FALSE TRUE
print(paste("Probability of taking the same interim decision (KRAS subgroup not dropped): ", table(rslt4_df$Z_score_interim>=0.34))) # FALSE TRUE






