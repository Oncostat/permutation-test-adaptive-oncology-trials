# permutation-test-adaptive-oncology-trials

## Article: "Using permutation-based tests to evaluate an adaptive molecular treatment algorithm in randomized precision oncology trials" (under review)
Authors: Cristina Marelli, Xavier Paoletti, and Stefan Michiels

This repository GitHub includes:
* "simulation_studies.R", which allows to replicate the different simulation studies presented in the paper (examples reported at the end of the file or below) and perform additional simulation studies with different parameters;
* "permutation-test-adaptive-oncology-trials.Rshiny", an R Shiny app developed based on the functions in "simulation_studies.R". To run the app, please follow these steps:
1. Download the auxiliary file "aux_functions.R" and save it in the same folder as the main .Rshiny file (this is recommended);
2. Open the .Rshiny file in RStudio;
3. Click on the "Run App" button to launch the application.
* "clinical_trial_application.R", which allows to replicate the two clinical trial applications presented in the paper. To run these, you need to download "syn_BREAST.xlsx" and "syn_LUNG.xsx", which are two synthetic datasets generated from the SAFIR02-BREAST and the SAFIR02-LUNG datasets, respectively.

### Examples of how to use the functions in "simulation_studies.R":

#### 1. Continuous outcomes, EXT-DROP scenario, under H0

```r
future::availableCores() # check the number of available cores 
plan(multisession, workers = future::availableCores() - 2)

seed <- 1234
Ntrials <- 10000
Nperms <- 1000
Natts <- 10000

rslt_cont_EXT_DROP_H0 <- furrr::future_map(1:Ntrials, ~{
  simulate_one(type_outcome = "continuous",
               type_rand = "block",
               type_adapt = "EXT",
               stageII_approach = "DROP",
               N = 200, N1 = 100, M = 4, prevalence = c(1/4,1/4,1/4,1/4), trt_labels = c(1,2),
               probs = c(0.5,0.5), bs = 4, p_min = NULL, adapt_time = NULL, Tr = NULL, 
               Ts = NULL, delta = NULL, alpha = 1, lambda0 = NULL, sigma = c(0.5,0.5,0.5,0.5), 
               coefX = NULL, probX = NULL,
               tau = matrix(c(0,0,0,0,
                              0,0,0,0,
                              0,0,0,0,
                              0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE),
               drop_trt = 4, z_c = NULL,
               probs_orient = NULL,
               L = Nperms, max_atts = NULL) %>%
    mutate(sim=.x)
}, .progress = TRUE, .options = furrr_options(seed = seed))
rslt_cont_EXT_DROP_H0_df <- list_rbind(rslt_cont_EXT_DROP_H0)

print("Type I errors:")
print(paste("Naive t-test:", sum(rslt_cont_EXT_DROP_H0_df$rejectsNAIVET == 1, na.rm = TRUE)/
              sum(!is.na(rslt_cont_EXT_DROP_H0_df$rejectsNAIVET))*100, "%"))
print(paste("Asymptotic F-test:", sum(rslt_cont_EXT_DROP_H0_df$F1_rejectsWALD == 1, na.rm = TRUE)/
              sum(!is.na(rslt_cont_EXT_DROP_H0_df$F1_rejectsWALD))*100, "%"))
print(paste("Convergence issues in asymptotic:", sum(rslt_cont_EXT_DROP_H0_df$F1_nullWALD, na.rm = TRUE)/Ntrials*100, "%"))
print(paste("Permutation F-test:", sum(rslt_cont_EXT_DROP_H0_df$F1_rejectsPERM == 1, na.rm = TRUE)/
              sum(!is.na(rslt_cont_EXT_DROP_H0_df$F1_rejectsPERM))*100, "%"))
print(paste("Convergence issues in permutations:", sum(rslt_cont_EXT_DROP_H0_df$F1_nullsPERM, na.rm = TRUE)/
              ((Ntrials - sum(rslt_cont_EXT_DROP_H0_df$F1_nullWALD, na.rm = TRUE))*Nperms)*100, "%"))
```

#### 2. Continuous outcomes, EXT-ORIENT scenario, under H1

```r
rslt_cont_EXT_ORIENT_H1 <- furrr::future_map(1:Ntrials, ~{
  simulate_one(type_outcome = "continuous",
               type_rand = "block",
               type_adapt = "EXT",
               stageII_approach = "ORIENT",
               N = 200, N1 = 100, M = 4, prevalence = c(1/4,1/4,1/4,1/4), trt_labels = c(1,2),
               probs = c(0.5,0.5), bs = 4, p_min = NULL, adapt_time = NULL, Tr = NULL, 
               Ts = NULL, delta = NULL, alpha = 1, lambda0 = NULL, sigma = c(0.5,0.5,0.5,0.5), 
               coefX = NULL, probX = NULL,
               tau = matrix(c(0.4,0,0,0,
                              0,0.4,0,0,
                              0,0,0.4,0,
                              0.4,0.4,0.4,0.4), nrow = 4, ncol = 4, byrow = TRUE),
               drop_trt = 4, z_c = NULL, 
               probs_orient = matrix(c(0,0,0,0,
                                       0,0,0,0,
                                       0,0,0,0,
                                       1/3,1/3,1/3,0), nrow = 4, ncol = 4, byrow = TRUE),
               L = Nperms, max_atts = NULL) %>%
    mutate(sim=.x)
}, .progress = TRUE, .options = furrr_options(seed = seed))
rslt_cont_EXT_ORIENT_H1_df <- list_rbind(rslt_cont_EXT_ORIENT_H1)

print("Statistical powers:")
print(paste("Naive t-test:", sum(rslt_cont_EXT_ORIENT_H1_df$rejectsNAIVET == 1, na.rm = TRUE)/
              sum(!is.na(rslt_cont_EXT_ORIENT_H1_df$rejectsNAIVET))*100, "%"))
print(paste("Asymptotic F-test (F1):", sum(rslt_cont_EXT_ORIENT_H1_df$F1_rejectsWALD == 1, na.rm = TRUE)/
              sum(!is.na(rslt_cont_EXT_ORIENT_H1_df$F1_rejectsWALD))*100, "%"))
print(paste("Convergence issues in asymptotic (F1):", sum(rslt_cont_EXT_ORIENT_H1_df$F1_nullWALD, na.rm = TRUE)/Ntrials*100, "%"))
print(paste("Permutation F-test (F1):", sum(rslt_cont_EXT_ORIENT_H1_df$F1_rejectsPERM == 1, na.rm = TRUE)/
              sum(!is.na(rslt_cont_EXT_ORIENT_H1_df$F1_rejectsPERM))*100, "%"))
print(paste("Convergence issues in permutations (F1):", sum(rslt_cont_EXT_ORIENT_H1_df$F1_nullsPERM, na.rm = TRUE)/
              ((Ntrials - sum(rslt_cont_EXT_ORIENT_H1_df$F1_nullWALD, na.rm = TRUE))*Nperms)*100, "%"))
print(paste("Asymptotic F-test (F2):", sum(rslt_cont_EXT_ORIENT_H1_df$F2_rejectsWALD == 1, na.rm = TRUE)/
              sum(!is.na(rslt_cont_EXT_ORIENT_H1_df$F2_rejectsWALD))*100, "%"))
print(paste("Convergence issues in asymptotic (F2):", sum(rslt_cont_EXT_ORIENT_H1_df$F2_nullWALD, na.rm = TRUE)/Ntrials*100, "%"))
print(paste("Permutation F-test (F2):", sum(rslt_cont_EXT_ORIENT_H1_df$F2_rejectsPERM == 1, na.rm = TRUE)/
              sum(!is.na(rslt_cont_EXT_ORIENT_H1_df$F2_rejectsPERM))*100, "%"))
print(paste("Convergence issues in permutations (F2):", sum(rslt_cont_EXT_ORIENT_H1_df$F2_nullsPERM, na.rm = TRUE)/
              ((Ntrials - sum(rslt_cont_EXT_ORIENT_H1_df$F2_nullWALD, na.rm = TRUE))*Nperms)*100, "%"))
print(paste("Asymptotic F-test (F3):", sum(rslt_cont_EXT_ORIENT_H1_df$F3_rejectsWALD == 1, na.rm = TRUE)/
              sum(!is.na(rslt_cont_EXT_ORIENT_H1_df$F3_rejectsWALD))*100, "%"))
print(paste("Convergence issues in asymptotic (F3):", sum(rslt_cont_EXT_ORIENT_H1_df$F3_nullWALD, na.rm = TRUE)/Ntrials*100, "%"))
print(paste("Permutation F-test (F3):", sum(rslt_cont_EXT_ORIENT_H1_df$F3_rejectsPERM == 1, na.rm = TRUE)/
              sum(!is.na(rslt_cont_EXT_ORIENT_H1_df$F3_rejectsPERM))*100, "%"))
print(paste("Convergence issues in permutations (F3):", sum(rslt_cont_EXT_ORIENT_H1_df$F3_nullsPERM, na.rm = TRUE)/
              ((Ntrials - sum(rslt_cont_EXT_ORIENT_H1_df$F3_nullWALD, na.rm = TRUE))*Nperms)*100, "%"))
```

#### 3. Time-to-event outcomes, IA-DROP scenario H0

```r
rslt_surv_IA_DROP_H0 <- furrr::future_map(1:Ntrials, ~{
  simulate_one(type_outcome = "survival",
               type_rand = "minim",
               type_adapt = "IA",
               stageII_approach = "DROP",
               N = 200, N1 = 100, M = 4, prevalence = c(1/4,1/4,1/4,1/4), trt_labels = c(1,2),
               probs = c(0.5,0.5), bs = NULL, p_min = 0.9, adapt_time = 36, Tr = 72, 
               Ts = 84, delta = 0.01, alpha = NULL, lambda0 = log(2)/4, sigma = c(0.5,0.5,0.5,0.5), 
               coefX = 0.3, probX = 0.4,
               tau = matrix(c(0,0,0,0,
                              0,0,0,0,
                              0,0,0,0,
                              0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE),
               drop_trt = NULL, z_c = 0.9,
               probs_orient = NULL,
               L = Nperms, max_atts = Natts) %>%
    mutate(sim=.x)
}, .progress = TRUE, .options = furrr_options(seed = seed))
rslt_surv_IA_DROP_H0_df <- list_rbind(rslt_surv_IA_DROP_H0)

print("Dropping probabilities:")
print(table(rslt_surv_IA_DROP_H0_df$drop_trtIA)/Ntrials*100)

valid_adj <- !(rslt_surv_IA_DROP_H0_df$conv_issue_asympt == 1 | rslt_surv_IA_DROP_H0_df$na_issue_asympt == 1)
n_valid_adj <- sum(valid_adj, na.rm = TRUE)

print("Type I errors:")
print(paste("Asymptotic unadjusted score test:",
            sum(rslt_surv_IA_DROP_H0_df$rejects_unadj_asympt == 1, na.rm = TRUE)/Ntrials*100, "%"))
print(paste("Mean number of attempts in randomizations:", mean(rslt_surv_IA_DROP_H0_df$Nattempts, na.rm = T)))
print(paste("Randomization unadjusted score test:",
            sum(rslt_surv_IA_DROP_H0_df$rejects_unadj_rand == 1, na.rm = TRUE)/Ntrials*100, "%"))
print(paste("Asymptotic adjusted score test:",
            sum(rslt_surv_IA_DROP_H0_df$rejects_adj_asympt == 1, na.rm = TRUE)/n_valid_adj*100, "%"))
print(paste("Convergence issues in asymptotic:",
            (Ntrials - n_valid_adj)/Ntrials*100, "%"))
print(paste("Randomization adjusted score test:",
            sum(rslt_surv_IA_DROP_H0_df$rejects_adj_rand == 1, na.rm = TRUE)/
              sum(!is.na(rslt_surv_IA_DROP_H0_df$rejects_adj_rand))*100, "%"))
print(paste("Convergence issues in randomizations:",
            sum(rslt_surv_IA_DROP_H0_df$conv_issues_rand + rslt_surv_IA_DROP_H0_df$na_issues_rand, na.rm = TRUE)/
              (n_valid_adj*Nperms)*100, "%"))
```
