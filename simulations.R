## This R code allows to replicate the two simulation studies presented in the paper 
## "Using permutation-based tests to evaluate an adaptive molecular treatment
##  algorithm in randomized precision oncology trials" by Marelli et al.

## Install and upload packages and functions

# install.packages("car")
library(car)
# install.packages("tidyverse")
library(tidyverse)
# install.packages("furrr")
library(furrr)

## User-defined functions

## Stratified block randomization: stratification within strata and stages.
## Arguments:
## * the 'dataset' to be randomized must contain the 'stage' and 'strata' columns;
## * 'trt_seq' contains the labels of the possible treatments;
## * 'prob_each' is the vector of probabilities of assigning a certain treatment;
## * patients are assigned in block of fixed size 'blocksize'.
my_strat_block_rand <- function(dataset, trt_seq = c(1,2), prob_each = c(0.5,0.5), blocksize = 4) {
  
  # controls
  if (is.null(dataset)) 
    stop("Missing dataset")
  if (!"stage" %in% names(dataset)) 
    stop("The dataset must contains a column 'stage'")
  if (!"stratum" %in% names(dataset)) 
    stop("The dataset must contains a column 'stratum'")
  if (length(prob_each) != length(trt_seq)) 
    stop("'prob_each' and 'trt_seq' must have the same length")
  if (!is.numeric(prob_each) || any(prob_each < 0)) 
    stop("All entries of 'prob_each' must be non-negative numbers")
  if (abs(sum(prob_each) - 1) > 1e-10) 
    stop("The entries of 'prob_each' must sum to 1")
  if (length(blocksize)!=1 || !is.numeric(blocksize) || blocksize <= 0 || 
      blocksize != round(blocksize)) 
    stop("'blocksize' must be a positive integer")
  
  # FOR EACH STRATUM assign patients into BLOCKS of size 'blocksize'
  stratum_ID <- unique(dataset$stratum) 
  result_list <- list()
  
  for(i in 1:length(stratum_ID)){
    stratum_i <- dataset[dataset$stratum==stratum_ID[i],]
    n_stratum_i <- nrow(stratum_i)
    
    random_order <- sample(1:n_stratum_i)
    stratum_i_shuffled <- stratum_i[random_order, ]
    
    n_complete_blocks <- as.integer(n_stratum_i/blocksize) 
    n_remaining <- n_stratum_i %% blocksize
    
    block_assign <- c()
    if (n_complete_blocks > 0) {
      block_assign <- rep(1:n_complete_blocks, each = blocksize)
    }
    if (n_remaining > 0) {
      block_assign <- c(block_assign, rep(n_complete_blocks + 1, n_remaining))
    }
    
    # and FOR EACH BLOCK assign patients to the treatments based on 'prob_each' 
    block_ID <- unique(block_assign) 
    trt_assign <- numeric(n_stratum_i)
    
    for (j in 1:length(block_ID)){
      block_j <- which(block_assign == j)
      n_block_j <- length(block_j)
      
      n_each_trt <- floor(n_block_j * prob_each) 
      n_assigned <- sum(n_each_trt)
      n_remainder <- n_block_j - n_assigned
      
      if (n_remainder > 0) {
        prob_each_fixed <- ((prob_each * n_block_j) - n_each_trt) / n_remainder
        block_trt <- c(rep(trt_seq, n_each_trt),
                       sample(trt_seq, n_remainder, prob = prob_each_fixed, replace = F))
      } 
      else{
        block_trt <- rep(trt_seq, n_each_trt)
      }
      
      if(length(block_trt)==1){
        trt_assign[block_j] <- block_trt
      }
      else{
        trt_assign[block_j] <- sample(block_trt)
      }
    }
    
    stratum_rslt <- stratum_i_shuffled
    stratum_rslt$block <- paste0("M",stratum_ID[i], "_B", block_assign, 
                                 "_P", stratum_i_shuffled$stage[i])
    stratum_rslt$trt <- trt_assign
    
    result_list[[i]] <- stratum_rslt
  }
  
  dataset2 <- do.call(rbind, result_list)
  rownames(dataset2) <- NULL
  
  return(data.frame(dataset2))
}

## Orientation of patients belonging to the dropped strata through block randomization.
## Arguments:
## * the non-missing 'dataset' to orient;
## * 'strata_seq' contains the labels of the possible strata:
## * 'prob_each_orient' is the vector of probabilities of assigning a certain stratum.
my_block_rand <- function(dataset, strata_seq = c(1,2,3), prob_each_orient = c(1/3,1/3,1/3)) {
  
  # controls
  if (is.null(dataset)) 
    stop("Missing dataset")
  if (length(prob_each_orient) != length(strata_seq)) 
    stop("'prob_each_orient' and 'strata_seq' must have the same length")
  if (!is.numeric(prob_each_orient) || any(prob_each_orient < 0)) 
    stop("All entries of 'prob_each_orient' must be non-negative numbers")
  if (abs(sum(prob_each_orient) - 1) > 1e-10) 
    stop("The entries of 'prob_each_orient' must sum to 1")
  
  n_tot <- nrow(dataset)
  
  random_order <- sample(1:n_tot)
  data_shuffled <- dataset[random_order, ]
  
  blocksize <- length(strata_seq)
  
  n_complete_blocks <- as.integer(n_tot/blocksize) 
  n_remaining <- n_tot %% blocksize
  
  block_assign <- c()
  if (n_complete_blocks > 0) {
    block_assign <- rep(1:n_complete_blocks, each = blocksize)
  }
  if (n_remaining > 0) {
    block_assign <- c(block_assign, rep(n_complete_blocks + 1, n_remaining))
  }
  
  # for each block assign patients to strata based on 'prob_each_orient' 
  block_ID <- unique(block_assign) 
  stratum_assign <- numeric(n_tot)
  
  for (j in 1:length(block_ID)){
    block_j <- which(block_assign == j)
    n_block_j <- length(block_j)
    
    n_each_stratum <- floor(n_block_j * prob_each_orient) 
    n_assigned <- sum(n_each_stratum)
    n_remainder <- n_block_j - n_assigned
    
    if (n_remainder > 0) {
      prob_each_orient_fixed <- ((prob_each_orient * n_block_j) - n_each_stratum)/n_remainder
      block_stratum <- c(rep(strata_seq, n_each_stratum),
                         sample(strata_seq, n_remainder, prob = prob_each_orient_fixed, 
                                replace = F))
    } 
    else{
      block_stratum <- rep(strata_seq, n_each_stratum)
    }
    
    if(length(block_stratum)==1){
      stratum_assign[block_j] <- block_stratum
    }
    else{
      stratum_assign[block_j] <- sample(block_stratum)
    }
  }
  
  data_rslt <- data_shuffled
  data_rslt$stratum <- stratum_assign
  rownames(data_rslt) <- NULL
  
  return(data.frame(data_rslt))
}

## Permutation of treatment assignments within blocks.
## Arguments:
## * the 'dataset' to be permuted must contain the columns 'trt' and 'block';
## * a 'seed' can be set for reproducibility.
permute_trts_by_block <- function(dataset, seed = NULL) {
  
  # controls
  if (is.null(dataset)) 
    stop("Missing dataset")
  if (!"trt" %in% names(dataset)) 
    stop("The dataset must contains a column 'trt'")
  if (!"block" %in% names(dataset)) 
    stop("The dataset must contains a column 'block'")
  
  if (!is.null(seed)) set.seed(seed)
  
  permuted_dataset <- dataset
  
  unique_blocks <- unique(dataset$block)
  
  for (block_id in unique_blocks) {
    
    block_indices <- which(dataset$block == block_id)
    
    current_treatments <- dataset$trt[block_indices]
    
    if(length(block_indices)==1){
      permuted_treatments <- current_treatments
    }
    else{
      permuted_treatments <- sample(current_treatments)
    }
    
    permuted_dataset$trt[block_indices] <- permuted_treatments
  }
  
  return(permuted_dataset)
}

## F-statistic specification 1: no account for shifts in prognosis or treatment effect.
## Arguments: the 'dataset' to be analysed must contain the columns 'trt' and 'stratum'.
analysis1 <- function(dataset){
  
  ## controls
  if (is.null(dataset)) 
    stop("Missing dataset")
  if (!"trt" %in% names(dataset)) 
    stop("The dataset must contains a column 'trt'")
  if (!"stratum" %in% names(dataset)) 
    stop("The dataset must contains a column 'stratum'")
  
  # definition of dummy variables and corresponding weights
  stratum_levels <- unique(dataset$stratum)
  trt_levels <- unique(dataset$trt)
  for(s in stratum_levels) {
    for(t in trt_levels) {
      dummy_name <- paste0(ifelse(t == trt_levels[1], "C", "T"), s)
      dataset[[dummy_name]] <- as.numeric(dataset$trt == t & dataset$stratum == s)
    }
  }
  weights <- sapply(stratum_levels, function(s) {
    sum(dataset[[paste0("C", s)]]) + sum(dataset[[paste0("T", s)]])
  }) / nrow(dataset)
  names(weights) <- paste0("w", stratum_levels)
  dummy_vars <- paste0(rep(c("C", "T"), each = length(stratum_levels)),
                       rep(stratum_levels, 2))
  # linear regression analysis formulation
  formula_str <- paste0("~ 0 + ", paste(dummy_vars, collapse = " + "))
  # hypothesis
  hypothesis_terms <- sapply(stratum_levels, function(s) {
    w <- weights[paste0("w", s)]
    paste0(w, "*T", s, " - ", w, "*C", s)
  })
  hypothesis <- paste0(paste(hypothesis_terms, collapse = " + "), " = 0")
  
  return(list(db = dataset, formula = formula_str, hypoth = hypothesis))
  
}

## F-statistic specification 2: account for shift in prognosis only.
## Arguments: the 'dataset' to be analysed must contain the columns 'trt', 'stratum', 
## 'stratum_old', and 'stage'.
analysis2 <- function(dataset){
  
  ## controls
  if (is.null(dataset)) 
    stop("Missing dataset")
  if (!"trt" %in% names(dataset)) 
    stop("The dataset must contains a column 'trt'")
  if (!"stratum" %in% names(dataset)) 
    stop("The dataset must contains a column 'stratum'")
  if (!"stratum_old" %in% names(dataset)) 
    stop("The dataset must contains a column 'stratum_old'")
  if (!"stage" %in% names(dataset)) 
    stop("The dataset must contains a column 'stage'")
  
  # definition of dummy variables and corresponding weights
  stratum_levels <- unique(dataset$stratum)
  trt_levels <- unique(dataset$trt)
  stratum_old_levels <- sort(unique(dataset$stratum_old))
  for(s in stratum_levels) {
    for(t in trt_levels) {
      dummy_name <- paste0(ifelse(t == trt_levels[1], "C", "T"), s)
      dataset[[dummy_name]] <- as.numeric(dataset$trt == t & dataset$stratum == s)
    }
  }
  period_dummy_vars <- c()
  existing_combinations <- unique(dataset[dataset$stage == 2, c("stratum_old", "stratum")])
  for(i in 1:nrow(existing_combinations)) {
    s_old <- existing_combinations$stratum_old[i]
    s <- existing_combinations$stratum[i]
    dummy_name <- paste0("P", s_old, s, "II")
    dataset[[dummy_name]] <- as.numeric(dataset$stratum_old == s_old & 
                                          dataset$stratum == s & 
                                          dataset$stage == 2)
    period_dummy_vars <- c(period_dummy_vars, dummy_name)
  }
  weights <- sapply(stratum_levels, function(s) {
    sum(dataset[[paste0("C", s)]]) + sum(dataset[[paste0("T", s)]])
  }) / nrow(dataset)
  names(weights) <- paste0("w", stratum_levels)
  trt_dummy_vars <- paste0(rep(c("C", "T"), each = length(stratum_levels)),
                           rep(stratum_levels, 2))
  all_dummy_vars <- c(trt_dummy_vars, period_dummy_vars)
  # linear regression analysis formulation
  formula_str <- paste0("~ 0 + ", paste(all_dummy_vars, collapse = " + "))
  # hypothesis
  hypothesis_terms <- sapply(stratum_levels, function(s) {
    w <- weights[paste0("w", s)]
    paste0(w, "*T", s, " - ", w, "*C", s)
  })
  hypothesis <- paste0(paste(hypothesis_terms, collapse = " + "), " = 0")
  
  return(list(db = dataset, formula = formula_str, hypoth = hypothesis))
  
}


## F-statistic specification 3: account for shifts in prognosis and treatment effect.
## Arguments: the 'dataset' to be analysed must contain the columns 'trt', 'stratum', 
## 'stratum_old', and 'stage'.
analysis3 <- function(dataset){
  
  ## controls
  if (is.null(dataset)) 
    stop("Missing dataset")
  if (!"trt" %in% names(dataset)) 
    stop("The dataset must contains a column 'trt'")
  if (!"stratum" %in% names(dataset)) 
    stop("The dataset must contains a column 'stratum'")
  if (!"stratum_old" %in% names(dataset)) 
    stop("The dataset must contains a column 'stratum_old'")
  if (!"stage" %in% names(dataset)) 
    stop("The dataset must contains a column 'stage'")
  
  # definition of dummy variables and corresponding weights
  stratum_levels <- sort(unique(dataset$stratum))
  trt_levels <- sort(unique(dataset$trt))
  stratum_old_levels <- sort(unique(dataset$stratum_old))
  for(s in stratum_levels) {
    for(t in trt_levels) {
      dummy_name <- paste0(ifelse(t == trt_levels[1], "C", "T"), s, "I")
      dataset[[dummy_name]] <- as.numeric(dataset$trt == t & 
                                            dataset$stratum == s & 
                                            dataset$stage == 1)
    }
  }
  stage2_dummy_vars <- c()
  existing_combinations <- unique(dataset[dataset$stage == 2, c("stratum_old", "stratum")])
  for(i in 1:nrow(existing_combinations)) {
    s_old <- existing_combinations$stratum_old[i]
    s <- existing_combinations$stratum[i]
    for(t in trt_levels) {
      dummy_name <- paste0(ifelse(t == trt_levels[1], "C", "T"), s_old, s, "II")
      dataset[[dummy_name]] <- as.numeric(dataset$trt == t & 
                                            dataset$stratum_old == s_old & 
                                            dataset$stratum == s & 
                                            dataset$stage == 2)
      stage2_dummy_vars <- c(stage2_dummy_vars, dummy_name)
    }
  }
  n_strata <- list()
  for(s in stratum_levels) {
    n_strata[[paste0("n", s, "I")]] <- sum(dataset[[paste0("C", s, "I")]]) + 
      sum(dataset[[paste0("T", s, "I")]])
    
    n_II <- 0
    for(s_old in stratum_old_levels) {
      c_var <- paste0("C", s_old, s, "II")
      t_var <- paste0("T", s_old, s, "II")
      if(c_var %in% names(dataset)) n_II <- n_II + sum(dataset[[c_var]])
      if(t_var %in% names(dataset)) n_II <- n_II + sum(dataset[[t_var]])
    }
    n_strata[[paste0("n", s, "II")]] <- n_II
  }
  weights_I <- list()
  for(s in stratum_levels) {
    n_sI <- n_strata[[paste0("n", s, "I")]]
    if(n_sI > 0) {
      w <- ((sum(dataset[[paste0("C", s, "I")]]) + 
               sum(dataset[[paste0("T", s, "I")]])) / n_sI) / 2
      weights_I[[paste0("w", s, "I")]] <- w
    }
  }
  weights_II <- list()
  for(s in stratum_levels) {
    n_sII <- n_strata[[paste0("n", s, "II")]]
    for(s_old in stratum_old_levels) {
      c_var <- paste0("C", s_old, s, "II")
      t_var <- paste0("T", s_old, s, "II")
      if(c_var %in% names(dataset) && t_var %in% names(dataset)) {
        if(n_sII > 0) {
          w <- ((sum(dataset[[c_var]]) + sum(dataset[[t_var]])) / n_sII) / 2
          weights_II[[paste0("w", s_old, s, "II")]] <- w
        }
      }
    }
  }
  stage1_dummy_vars <- paste0(rep(c("C", "T"), each = length(stratum_levels)), 
                              rep(stratum_levels, 2), "I")
  all_dummy_vars <- c(stage1_dummy_vars, stage2_dummy_vars)
  # linear regression analysis formulation
  formula_str <- paste0("~ 0 + ", paste(all_dummy_vars, collapse = " + "))
  # hypothesis
  hypothesis_terms <- c()
  for(s in stratum_levels) {
    w_name <- paste0("w", s, "I")
    if(w_name %in% names(weights_I)) {
      w <- weights_I[[w_name]]
      term <- paste0(w, "*T", s, "I - ", w, "*C", s, "I")
      hypothesis_terms <- c(hypothesis_terms, term)
    }
  }
  for(s in stratum_levels) {
    for(s_old in stratum_old_levels) {
      w_name <- paste0("w", s_old, s, "II")
      if(w_name %in% names(weights_II)) {
        w <- weights_II[[w_name]]
        term <- paste0(w, "*T", s_old, s, "II - ", w, "*C", s_old, s, "II")
        hypothesis_terms <- c(hypothesis_terms, term)
      }
    }
  }
  hypothesis <- paste0(paste(hypothesis_terms, collapse = " + "), " = 0")
  
  return(list(db = dataset, formula = formula_str, hypoth = hypothesis))
  
}

## Permuted p-value estimation for unplanned changes.
## Arguments:
## * the 'dataset' represents the full dataset to be analysed and must contain 
## the columns 'trt', 'stratum', and 'y_H0' or 'y_H1';
## * 'stat_obs' is the observed test statistic;
## * 'hyp' equals 0 or 1 if under the null (H0) or alternative (H1) hypothesis, respectively;
## * 'L' is the number of permutations;
## * 'analysis' is the type of F-statistic specification to be considered.
perm_test_pval_EXT <- function(dataset, stat_obs, hyp = 0, L = 1000, analysis = 0){
  
  ## controls
  if (is.null(dataset)) 
    stop("Missing dataset")
  if (!"trt" %in% names(dataset)) 
    stop("The dataset must contains a column 'trt'")
  if (!"stratum" %in% names(dataset)) 
    stop("The dataset must contains a column 'stratum'")
  if (length(stat_obs)!=1 || !is.numeric(stat_obs)) 
    stop("'stat_obs' must be numeric")
  if (!(hyp %in% c(0,1))) 
    stop("'hyp' must be 0 or 1")
  if (length(L)!=1 || !is.numeric(L) || L <= 0 || L != round(L)) 
    stop("'L' must be a positive integer")
  if (!(analysis %in% c(0,1,2))) 
    stop("'analysis' must be 0 or 1 or 2")
  
  wald_L <- rep(0,L)
  I_L <- rep(0,L)
  dataset_copy <- dataset
  for (l in 1:L){
    # treatment labels permutation
    new_dataset <- permute_trts_by_block(dataset = dataset_copy) 
    new_dataset <- new_dataset[order(new_dataset$ID), ]
    # treatment effects estimation
    if (analysis==0){
      new_analysis <- analysis1(dataset = new_dataset)
    }
    else if (analysis==1){
      new_analysis <- analysis2(dataset = new_dataset)
    }
    else if (analysis==2){
      new_analysis <- analysis3(dataset = new_dataset)
    }
    new_dataset_dummy <- new_analysis$db
    new_formula_str <- new_analysis$formula
    if(hyp==0){
      if (!"y_H0" %in% names(new_dataset_dummy)) 
        stop("The dataset must contains a column 'y_H0'")
      model_l <- lm(as.formula(paste("y_H0", new_formula_str)), 
                    data = new_dataset_dummy)
    }
    else if(hyp==1){
      if (!"y_H1" %in% names(new_dataset_dummy)) 
        stop("The dataset must contains a column 'y_H1'")
      model_l <- lm(as.formula(paste("y_H1", new_formula_str)), 
                    data = new_dataset_dummy)
    }
    if(any(is.na(coef(model_l)))) {
      wald_L[l] <- NA
      I_L[l] <- NA
    } 
    else{
      new_hypothesis <- new_analysis$hypoth
      # computation of the test statistic at this iteration
      wald_L[l] <- linearHypothesis(model_l, new_hypothesis)$F[2] 
      # checking the condition |Sl|>=|Sobs| at this iteration
      I_L[l] <- ifelse(abs(wald_L[l])>=abs(stat_obs),1,0) 
    }
  }
  # computation two-sided p-value
  pvalue <- (1/sum(!is.na(I_L)))*sum(I_L, na.rm = T) 
  # numbers of permutations with model convergence issues
  nulls <- sum(is.na(I_L))
  
  return(list(pvalue = pvalue, nulls = nulls))
}

## Permuted p-value estimation for interim analysis.
## Arguments:
## * 'dataset1' and 'dataset2' represent the dataset before and after the interim 
## analysis and must contain the columns 'trt', 'stratum', and 'y_H0' or 'y_H1';
## * 'crit_val' is the pre-planned futility threshold (can be Inf if the interim 
## futility analysis is not planned);
## * 'drop_trt' is the treatment dropped at interim on obseved data;
## * 'stat_obs' is the observed test statistic;
## * 'hyp' equals 0 or 1 if under the null (H0) or alternative (H1) hypothesis, respectively;
## * 'L' is the number of permutations;
## * 'max_attempts' is the maximum number of attempts to compute the permutation 
## test (it can be Inf);
## * 'analysis' is the type of F-statistic specification to be considered.
perm_test_pval_IA <- function(dataset1, dataset2, crit_val = 0, drop_trt = 0, stat_obs, 
                              hyp = 0, L = 1000, max_attempts = 100000, analysis = 0){
  
  ## controls
  if (is.null(dataset1)) 
    stop("Missing dataset for stage I")
  if (!"trt" %in% names(dataset1)) 
    stop("The dataset1 must contains a column 'trt'")
  if (!"stratum" %in% names(dataset1)) 
    stop("The dataset1 must contains a column 'stratum'")
  if (is.null(dataset2)) 
    stop("Missing dataset for stage II")
  if (!"trt" %in% names(dataset2)) 
    stop("The dataset2 must contains a column 'trt'")
  if (!"stratum" %in% names(dataset2)) 
    stop("The dataset2 must contains a column 'stratum'")
  if (length(crit_val)!=1 || !is.numeric(crit_val) || is.nan(crit_val)) 
    stop("'crit_val' must be numeric or infinit")
  if (length(drop_trt)!=1 || !is.numeric(drop_trt)) 
    stop("'drop_trt' must be numeric")
  if (length(stat_obs)!=1 || !is.numeric(stat_obs)) 
    stop("'stat_obs' must be numeric")
  if (!(hyp %in% c(0,1))) 
    stop("'hyp' must be 0 or 1")
  if (length(L)!=1 || !is.numeric(L) || L <= 0 || L != round(L)) 
    stop("'L' must be a positive integer")
  if (length(max_attempts)!=1 || !is.numeric(max_attempts) || is.nan(max_attempts)) 
    stop("'max_attempts' must be numeric or infinit")
  if (!(analysis %in% c(0,1,2))) 
    stop("'analysis' must be 0 or 1 or 2")
  if (hyp==0){
    if (!"y_H0" %in% names(dataset1)) 
      stop("The dataset1 must contains a column 'y_H0'")
    if (!"y_H0" %in% names(dataset2)) 
      stop("The dataset2 must contains a column 'y_H0'")
  }
  else if (hyp==1){
    if (!"y_H1" %in% names(dataset1)) 
      stop("The dataset1 must contains a column 'y_H1'")
    if (!"y_H1" %in% names(dataset2)) 
      stop("The dataset2 must contains a column 'y_H1'")
  }
  
  wald_L <- rep(0,L)
  I_L <- rep(0,L)
  dataset_copy <- dataset1
  l <- 1  # valid permutations counter
  attempts <- 0 # number of attempts counter
  
  while(l <= L){
    attempts <- attempts + 1
    
    if(attempts > max_attempts){
      warning(paste("Only ", l-1, " permutations obtained out of ", L, " permutations."))
      wald_L <- wald_L[1:(l-1)]
      I_L <- I_L[1:(l-1)]
      break
    }
    
    # treatment labels permutation
    new_dataset1 <- permute_trts_by_block(dataset = dataset_copy) 
    new_dataset1 <- new_dataset1[order(new_dataset1$ID), ]
    
    new_trt_labels <- sort(unique(new_dataset1$trt))  
    new_ctrl_label <- new_trt_labels[1]
    new_trt_label <- new_trt_labels[2]
    new_M <- length(unique(new_dataset1$stratum))
    new_stddiff_H0 <- sapply(1:new_M, function(s) {
      new_stratum_ctrl <- new_dataset1[new_dataset1$stratum == s & 
                                         new_dataset1$trt == new_ctrl_label, ]
      new_stratum_trt <- new_dataset1[new_dataset1$stratum == s & 
                                        new_dataset1$trt == new_trt_label, ]
      new_mean_diff <- mean(new_stratum_trt$y_H0) - mean(new_stratum_ctrl$y_H0)
      new_sd_ctrl <- sd(new_stratum_ctrl$y_H0)
      return(new_mean_diff / new_sd_ctrl)
      # if the number of subjects and variability differ between the treatment and control groups, 
      # comment the two previous lines and uncomment the following: 
      # new_sd_pooled <- sqrt(((nrow(new_stratum_trt)-1)*var(new_stratum_trt$y_H0) +
      #                         (nrow(new_stratum_ctrl)-1)*var(new_stratum_ctrl$y_H0))/
      #                         (nrow(new_stratum_trt)+nrow(new_stratum_ctrl)-2))
      # new_denom <- new_sd_pooled*(1/nrow(new_stratum_trt)+1/nrow(new_stratum_ctrl))
      # return(new_mean_diff / new_denom)
    })
    new_stddiff_H1 <- sapply(1:new_M, function(s) {
      new_stratum_ctrl <- new_dataset1[new_dataset1$stratum == s & 
                                         new_dataset1$trt == new_ctrl_label, ]
      new_stratum_trt <- new_dataset1[new_dataset1$stratum == s & 
                                        new_dataset1$trt == new_trt_label, ]
      new_mean_diff <- mean(new_stratum_trt$y_H1) - mean(new_stratum_ctrl$y_H1)
      new_sd_ctrl <- sd(new_stratum_ctrl$y_H1)
      return(new_mean_diff / new_sd_ctrl)
      # if the number of subjects and variability differ between the treatment and control groups, 
      # comment the two previous lines and uncomment the following: 
      # new_sd_pooled <- sqrt(((nrow(new_stratum_trt)-1)*var(new_stratum_trt$y_H1) +
      #                         (nrow(new_stratum_ctrl)-1)*var(new_stratum_ctrl$y_H1))/
      #                         (nrow(new_stratum_trt)+nrow(new_stratum_ctrl)-2))
      # new_denom <- new_sd_pooled*(1/nrow(new_stratum_trt)+1/nrow(new_stratum_ctrl))
      # return(new_mean_diff / new_denom)
    })
    
    new_min_diff_H0 <- min(new_stddiff_H0)
    new_drops_H0 <- ifelse(new_stddiff_H0 == new_min_diff_H0 & 
                             new_stddiff_H0 <= crit_val, 1, 0)
    new_min_diff_H1 <- min(new_stddiff_H1)
    new_drops_H1 <- ifelse(new_stddiff_H1 == new_min_diff_H1 & 
                             new_stddiff_H1 <= crit_val, 1, 0)
    
    new_pooled_stratum_H0 <- which(new_drops_H0 == 0)
    new_pooled_stratum_H1 <- which(new_drops_H1 == 0)
    
    new_drop_trt_H0 <- ifelse(length(new_pooled_stratum_H0) == new_M, 0, 
                              setdiff(1:new_M, new_pooled_stratum_H0))
    new_drop_trt_H1 <- ifelse(length(new_pooled_stratum_H1) == new_M, 0, 
                              setdiff(1:new_M, new_pooled_stratum_H1))
    
    if(hyp==0){
      if(identical(drop_trt,new_drop_trt_H0)){
        new_dataset1_H0 <- new_dataset1[new_dataset1$stratum!=new_drop_trt_H0,]
        
        new_dataset2_H0 <- permute_trts_by_block(dataset = dataset2)
        new_dataset2_H0 <- new_dataset2_H0[order(new_dataset2_H0$ID), ]
        
        new_dataset_H0 <- rbind(new_dataset1_H0,new_dataset2_H0)
        
        # treatment effects estimation
        if (analysis==0){
          new_analysis_H0 <- analysis1(dataset = new_dataset_H0)
        }
        else if (analysis==1){
          new_analysis_H0 <- analysis2(dataset = new_dataset_H0)
        }
        else if (analysis==2){
          new_analysis_H0 <- analysis3(dataset = new_dataset_H0)
        }
        new_dataset_dummy_H0 <- new_analysis_H0$db
        new_formula_str_H0 <- new_analysis_H0$formula
        model_l <- lm(as.formula(paste("y_H0", new_formula_str_H0)), 
                      data = new_dataset_dummy_H0)
        if(any(is.na(coef(model_l)))) {
          wald_L[l] <- NA
          I_L[l] <- NA
        } 
        else{
          new_hypothesis_H0 <- new_analysis_H0$hypoth
          # computation of the test statistic at this iteration
          wald_L[l] <- linearHypothesis(model_l, new_hypothesis_H0)$F[2] 
          # checking the condition |Sl|>=|Sobs| at this iteration
          I_L[l] <- ifelse(abs(wald_L[l])>=abs(stat_obs),1,0) 
        }
        
        l <- l + 1
      }
    }
    
    else if(hyp==1){
      if(identical(drop_trt,new_drop_trt_H1)){
        new_dataset1_H1 <- new_dataset1[new_dataset1$stratum!=new_drop_trt_H1,]
        
        new_dataset2_H1 <- permute_trts_by_block(dataset = dataset2)
        new_dataset2_H1 <- new_dataset2_H1[order(new_dataset2_H1$ID), ]
        
        new_dataset_H1 <- rbind(new_dataset1_H1,new_dataset2_H1)
        
        # treatment effects estimation
        if (analysis==0){
          new_analysis_H1 <- analysis1(dataset = new_dataset_H1)
        }
        else if (analysis==1){
          new_analysis_H1 <- analysis2(dataset = new_dataset_H1)
        }
        else if (analysis==2){
          new_analysis_H1 <- analysis3(dataset = new_dataset_H1)
        }
        new_dataset_dummy_H1 <- new_analysis_H1$db
        new_formula_str_H1 <- new_analysis_H1$formula
        if (!"y_H1" %in% names(new_dataset_dummy_H1)) 
          stop("The dataset must contains a column 'y_H1'")
        model_l <- lm(as.formula(paste("y_H1", new_formula_str_H1)), 
                      data = new_dataset_dummy_H1)
        if(any(is.na(coef(model_l)))) {
          wald_L[l] <- NA
          I_L[l] <- NA
        } 
        else{
          new_hypothesis_H1 <- new_analysis_H1$hypoth
          # computation of the test statistic at this iteration
          wald_L[l] <- linearHypothesis(model_l, new_hypothesis_H1)$F[2] 
          # checking the condition |Sl|>=|Sobs| at this iteration
          I_L[l] <- ifelse(abs(wald_L[l])>=abs(stat_obs),1,0) 
        }
        
        l <- l + 1
      }
    }
  }
  
  # computation two-sided p-value
  pvalue <- (1/sum(!is.na(I_L)))*sum(I_L, na.rm = T) 
  # numbers of permutations with model convergence issues
  nulls <- sum(is.na(I_L))
  
  return(list(pvalue = pvalue, nulls = nulls, total_attempts = attempts))
}

## Scenario 1 (external event) and exclusion of patients - first simulation:
## one (fixed) experimental treatment is dropped due to an external event.
## Arguments:
## * 'N' is the total number of planned patients;
## * 'M' is the number of molecular alterations considered at the beginning of the trial;
## * 'prevalence' is the vector of molecular alterations' prevalences;
## * 'trt_labels' is the vector of treatment labels (first entry for control treatment);
## * 'probs' is the vector of treatment probabilities;
## * 'bs' is the block size to be used for block randomization;
## * 'alpha' is the control effect;
## * 'sigma' is the vector of prognoses;
## * 'tau_H0' and 'tau_H1' are the vectors of treatment effects under the null (H0) and
## alternative (H1) hypothesis, respectively;
## * 'drop_trt' is the fixed experimental treatment dropped due to the external event;
## * 'perms' is the number of permutations.
simulate_one1 <- function(N = 200, M = 4, prevalence = c(1/4,1/4,1/4,1/4), trt_labels = c(1,2), 
                          probs = c(0.5,0.5), bs = 4, alpha = 1, sigma = c(0.5,0.5,0.5,0.5),
                          tau_H0 = c(0,0,0,0), tau_H1 = c(0.40,0.40,0.40,0.40), 
                          drop_trt = 4, perms = 1000)
  {
  
  # controls
  if (length(N)!=1 || !is.numeric(N) || N <= 0 || N != round(N)) 
    stop("'N' must be a positive integer")
  if (length(M)!=1 || !is.numeric(M) || M <= 0 || M != round(M)) 
    stop("'M' must be a positive integer")
  if (length(alpha)!=1 || !is.numeric(alpha)) 
    stop("'alpha' must be numeric")
  if (length(sigma) != M || length(tau_H0) != M || length(tau_H1) != M) 
    stop("'sigma', 'tau_H0', and 'tau_H1' must have length equal to 'M'")
  if (!is.numeric(sigma)) 
    stop("All entries of 'sigma' must be numeric")
  if (!is.numeric(tau_H0)) 
    stop("All entries of 'tau_H0' must be numeric")
  if (!is.numeric(tau_H1)) 
    stop("All entries of 'tau_H1' must be numeric")
  if (!is.null(drop_trt)) {
    if (!is.numeric(drop_trt) || length(drop_trt)!=1 || drop_trt>M) 
      stop("'drop_trt' must be NULL or numeric lower or equal than M") 
  }
  
  # STAGE I
  
  #  data generation
  molec_alter1 <- sample(1:M, prob = prevalence, size = N/2, replace = T)
  stage1 <- rep(1,N/2)
  data1 <- data.frame(stratum = molec_alter1, stage = stage1)
  N1 <- nrow(data1)
  data1$ID <- 1:N1
  
  #  treatment assignment
  data_strPBR1 <- my_strat_block_rand(dataset = data1, trt_seq = trt_labels, 
                                      prob_each = probs, blocksize = bs)
  data_strPBR1 <- data_strPBR1[order(data_strPBR1$ID), ]
  
  # continuous outcome generation
  data_strPBR1$y_H0 <- rep(0,N1)
  data_strPBR1$y_H1 <- rep(0,N1)
  for(j in 1:N1){
    data_strPBR1$err[j] <- rnorm(1,0,1)
    data_strPBR1$y_H0[j] <- alpha + sigma[data_strPBR1$stratum[j]] + 
      tau_H0[data_strPBR1$stratum[j]]*(data_strPBR1$trt[j]-1) + data_strPBR1$err[j]
    data_strPBR1$y_H1[j] <- alpha + sigma[data_strPBR1$stratum[j]] + 
      tau_H1[data_strPBR1$stratum[j]]*(data_strPBR1$trt[j]-1) + data_strPBR1$err[j]
  }
  
  
  # STAGE II
  
  #  data generation
  data_strPBR1_kept <- data_strPBR1[data_strPBR1$stratum!=drop_trt,]
  N1_k <- nrow(data_strPBR1_kept)
  N2 <- N - N1_k  
  molec_alter2 <- sample(setdiff(1:M, drop_trt), prob = prevalence[-drop_trt], N2, 
                         replace = T)
  stage2 <- rep(2, N2)
  data2 <- data.frame(stratum = molec_alter2, stage = stage2)
  data2$ID <- 1:N2
  
  #  treatment assignment
  data_strPBR2 <- my_strat_block_rand(dataset = data2, trt_seq = trt_labels, 
                                      prob_each = probs, blocksize = bs)
  data_strPBR2 <- data_strPBR2[order(data_strPBR2$ID), ]
  
  # outcome generation
  data_strPBR2$y_H0 <- rep(0,N2)
  data_strPBR2$y_H1 <- rep(0,N2)
  for(j in 1:N2){
    data_strPBR2$err[j] <- rnorm(1,0,1)
    data_strPBR2$y_H0[j] <- alpha + sigma[data_strPBR2$stratum[j]] + 
      tau_H0[data_strPBR2$stratum[j]]*(data_strPBR2$trt[j]-1) + data_strPBR2$err[j]
    data_strPBR2$y_H1[j] <- alpha + sigma[data_strPBR2$stratum[j]] + 
      tau_H1[data_strPBR2$stratum[j]]*(data_strPBR2$trt[j]-1) + data_strPBR2$err[j]
  }
  
  data_strPBR <- rbind(data_strPBR1_kept,data_strPBR2)
  Ntot <- nrow(data_strPBR)
  
  
  # t-test (NAIVE)
  # H0
  outNAIVET_H0 <- t.test(y_H0 ~ factor(trt), data = data_strPBR)
  pvaluesNAIVET_H0 <- outNAIVET_H0$p.value
  rejectsNAIVET <- ifelse(pvaluesNAIVET_H0<=0.05,1,0)
  
  # H1
  outNAIVET_H1 <- t.test(y_H1 ~ factor(trt), data = data_strPBR)
  pvaluesNAIVET_H1 <- outNAIVET_H1$p.value
  notrejectsNAIVET <- ifelse(pvaluesNAIVET_H1>0.05,1,0)
  
  
  # F-test (F-specification 1)
  analysis <- analysis1(dataset = data_strPBR)
  data_strPBR_dummy <- analysis$db
  formula_str <- analysis$formula
  hypothesis <- analysis$hypoth
  
  # H0
  modelWALD_H0 <- lm(as.formula(paste("y_H0", formula_str)), data = data_strPBR_dummy)
  if(any(is.na(coef(modelWALD_H0)))) {
    pvaluesWALD_H0 <- NA
    rejectsWALD <- NA
    pvaluesPERM_H0 <- NA
    rejectsPERM <- NA
    nulls_H0 <- NA
  } 
  else {
    outWALD_H0 <- linearHypothesis(modelWALD_H0, hypothesis)
    S_obsWALD_H0 <- outWALD_H0$F[2] # observed statistics under H0
    pvaluesWALD_H0 <- outWALD_H0$`Pr(>F)`[2]
    rejectsWALD <- ifelse(pvaluesWALD_H0<=0.05,1,0)
    
    # Permutation test (H0)
    outPERM_H0 <- perm_test_pval_EXT(dataset = data_strPBR, stat_obs = S_obsWALD_H0, hyp = 0, L = perms)
    pvaluesPERM_H0 <- outPERM_H0$pvalue
    rejectsPERM <- ifelse(pvaluesPERM_H0<=0.05,1,0)
    nulls_H0 <- outPERM_H0$nulls
  }
  
  # H1
  modelWALD_H1 <- lm(as.formula(paste("y_H1", formula_str)), data = data_strPBR_dummy)
  if(any(is.na(coef(modelWALD_H1)))) {
    pvaluesWALD_H1 <- NA
    notrejectsWALD <- NA
    pvaluesPERM_H1 <- NA
    notrejectsPERM <- NA
    nulls_H1 <- NA
  } 
  else {
    outWALD_H1 <- linearHypothesis(modelWALD_H1, hypothesis)
    S_obsWALD_H1 <- outWALD_H1$F[2] # observed statistics under H1
    pvaluesWALD_H1 <- outWALD_H1$`Pr(>F)`[2]
    notrejectsWALD <- ifelse(pvaluesWALD_H1>0.05,1,0)
    
    # Permutation test (H1)
    outPERM_H1 <- perm_test_pval_EXT(dataset = data_strPBR, stat_obs = S_obsWALD_H1, hyp = 1, L = perms)
    pvaluesPERM_H1 <- outPERM_H1$pvalue
    notrejectsPERM <- ifelse(pvaluesPERM_H1>0.05,1,0)
    nulls_H1 <- outPERM_H1$nulls
  }
  
  
  tibble::tibble(pvaluesNAIVET_H0, rejectsNAIVET, pvaluesNAIVET_H1, notrejectsNAIVET,
                 pvaluesWALD_H0, rejectsWALD, pvaluesWALD_H1, notrejectsWALD,
                 pvaluesPERM_H0, rejectsPERM, nulls_H0, pvaluesPERM_H1, notrejectsPERM, nulls_H1
  )
  
}

## Scenario 1 (external event) and orientation of patients - first simulation:
## one (fixed) experimental treatment is dropped due to an external event.
## Arguments:
## * 'N' is the total number of planned patients;
## * 'M' is the number of molecular alterations considered at the beginning of the trial;
## * 'prevalence' is the vector of molecular alterations' prevalences;
## * 'trt_labels' is the vector of treatment labels (first entry for control treatment);
## * 'probs' is the vector of treatment probabilities;
## * 'bs' is the block size to be used for block randomization;
## * 'alpha' is the control effect;
## * 'sigma' is the vector of prognoses;
## * 'tau_H0' and 'tau_H1' are the matrix of treatment effects under the null (H0) and
## alternative (H1) hypothesis, respectively;
## * 'drop_trt' is the fixed experimental treatment dropped due to the external event;
## * 'probs_orient' is the vector of probabilities of assigning a certain stratum;
## * 'perms' is the number of permutations.
simulate_one2 <- function(N = 200, M = 4, prevalence = c(1/4,1/4,1/4,1/4), trt_labels = c(1,2),
                          probs = c(0.5,0.5),bs = 4, alpha = 1, sigma = c(0.5,0.5,0.5,0.5), 
                          tau_H0 = matrix(c(0,0,0,0,
                                            0,0,0,0,
                                            0,0,0,0,
                                            0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE), 
                          tau_H1 = matrix(c(0.40,0,0,0,
                                            0,0.40,0,0,
                                            0,0,0.40,0,
                                            0.40,0.40,0.40,0.40), nrow = 4, ncol = 4, byrow = TRUE), 
                          drop_trt = 4, probs_orient = c(1/3,1/3,1/3), perms = 1000)
  {
  
  # controls
  if (length(N)!=1 || !is.numeric(N) || N <= 0 || N != round(N)) 
    stop("'N' must be a positive integer")
  if (length(M)!=1 || !is.numeric(M) || M <= 0 || M != round(M)) 
    stop("'M' must be a positive integer")
  if (length(alpha)!=1 || !is.numeric(alpha)) 
    stop("'alpha' must be numeric")
  if (length(sigma) != M) 
    stop("'sigma' must have length equal to 'M'")
  if (nrow(tau_H0) != M || ncol(tau_H0) != M) 
    stop("'tau_H0' must have dimention equal to 'MxM'")
  if (nrow(tau_H1) != M || ncol(tau_H1) != M) 
    stop("'tau_H1' must have dimention equal to 'MxM'")
  if (!is.numeric(sigma)) 
    stop("All entries of 'sigma' must be numeric")
  if (!is.numeric(tau_H0)) 
    stop("All entries of 'tau_H0' must be numeric")
  if (!is.numeric(tau_H1)) 
    stop("All entries of 'tau_H1' must be numeric")
  if (!is.null(drop_trt)) {
    if (!is.numeric(drop_trt) || length(drop_trt)!=1 || drop_trt>M) 
      stop("'drop_trt' must be NULL or numeric lower or equal than M") # extension to dropping 2+ treatments?
  }
  if (is.null(drop_trt)) {
    if (length(probs_orient)!=M) 
      stop("'probs_orient' must have length equal to 'M' if no strata is dropped")
  }
  else{
    if (length(probs_orient)!=M-1) 
      stop("'probs_orient' must have length equal to 'M'-1 if one strata is dropped")
  }
  
  # STAGE I
  
  #  data generation
  molec_alter1 <- sample(1:M, prob = prevalence, size = N/2, replace = T)
  stage1 <- rep(1,N/2)
  data1 <- data.frame(stratum_old = molec_alter1, stratum = molec_alter1, stage = stage1)
  N1 <- nrow(data1)
  data1$ID <- 1:N1
  
  #  treatment assignment
  data_strPBR1 <- my_strat_block_rand(dataset = data1, trt_seq = trt_labels, 
                                      prob_each = probs, blocksize = bs)
  data_strPBR1 <- data_strPBR1[order(data_strPBR1$ID), ]
  
  # continuous outcome generation
  data_strPBR1$y_H0 <- rep(0,N1)
  data_strPBR1$y_H1 <- rep(0,N1)
  for(j in 1:N1){
    data_strPBR1$err[j] <- rnorm(1,0,1)
    data_strPBR1$y_H0[j] <- alpha + sigma[data_strPBR1$stratum_old[j]] + 
      tau_H0[data_strPBR1$stratum_old[j],data_strPBR1$stratum[j]]*(data_strPBR1$trt[j]-1) + 
      data_strPBR1$err[j]
    data_strPBR1$y_H1[j] <- alpha + sigma[data_strPBR1$stratum_old[j]] + 
      tau_H1[data_strPBR1$stratum_old[j],data_strPBR1$stratum[j]]*(data_strPBR1$trt[j]-1) + 
      data_strPBR1$err[j]
  }
  
  
  # STAGE II
  
  #  data generation
  data_strPBR1_kept <- data_strPBR1[data_strPBR1$stratum!=drop_trt,]
  N1_k <- nrow(data_strPBR1_kept)
  N2 <- N - N1_k  
  molec_alter2 <- sample(1:M, prob = prevalence, N2, replace = T)
  stage2 <- rep(2, N2)
  data2 <- data.frame(stratum_old = molec_alter2, stage = stage2)
  data2$ID <- 1:N2
  
  # orientation of patients belonging to the dropped stratum (block randomization)
  data2_maintain <- data2[data2$stratum_old!=drop_trt,]
  data2_maintain$stratum <- data2_maintain$stratum_old
  data2_toorient <- data2[data2$stratum_old==drop_trt,]
  data2_oriented <- my_block_rand(data2_toorient, strata_seq = setdiff(1:M, drop_trt), 
                                  prob_each_orient = probs_orient)
  data2_toassign <- rbind(data2_maintain,data2_oriented)
  data2_toassign <- data2_toassign[order(data2_toassign$ID), ]
  
  #  treatment assignment
  data_strPBR2 <- my_strat_block_rand(dataset = data2_toassign, trt_seq = trt_labels, 
                                      prob_each = probs, blocksize = bs)
  data_strPBR2 <- data_strPBR2[order(data_strPBR2$ID), ]
  
  # outcome generation
  data_strPBR2$y_H0 <- rep(0,N2)
  data_strPBR2$y_H1 <- rep(0,N2)
  for(j in 1:N2){
    data_strPBR2$err[j] <- rnorm(1,0,1)
    data_strPBR2$y_H0[j] <- alpha + sigma[data_strPBR2$stratum_old[j]] + 
      tau_H0[data_strPBR2$stratum_old[j],data_strPBR2$stratum[j]]*(data_strPBR2$trt[j]-1) + 
      data_strPBR2$err[j]
    data_strPBR2$y_H1[j] <- alpha + sigma[data_strPBR2$stratum_old[j]] + 
      tau_H1[data_strPBR2$stratum_old[j],data_strPBR2$stratum[j]]*(data_strPBR2$trt[j]-1) + 
      data_strPBR2$err[j]
  }
  
  data_strPBR <- rbind(data_strPBR1_kept,data_strPBR2)
  Ntot <- nrow(data_strPBR)
  
  
  # t-test (NAIVE)
  # H0
  outNAIVET_H0 <- t.test(y_H0 ~ factor(trt), data = data_strPBR)
  pvaluesNAIVET_H0 <- outNAIVET_H0$p.value
  rejectsNAIVET <- ifelse(pvaluesNAIVET_H0<=0.05,1,0)
  
  # H1
  outNAIVET_H1 <- t.test(y_H1 ~ factor(trt), data = data_strPBR)
  pvaluesNAIVET_H1 <- outNAIVET_H1$p.value
  notrejectsNAIVET <- ifelse(pvaluesNAIVET_H1>0.05,1,0)
  
  
  # F-specification 1
  
  # F-test
  M1_analysis <- analysis1(dataset = data_strPBR)
  M1_data_strPBR_dummy <- M1_analysis$db
  M1_formula_str <- M1_analysis$formula
  M1_hypothesis <- M1_analysis$hypoth
  
  # H0
  M1_modelWALD_H0 <- lm(as.formula(paste("y_H0", M1_formula_str)), data = M1_data_strPBR_dummy)
  if(any(is.na(coef(M1_modelWALD_H0)))) {
    M1_pvaluesWALD_H0 <- NA
    M1_rejectsWALD <- NA
    M1_pvaluesPERM_H0 <- NA
    M1_rejectsPERM <- NA
    M1_nulls_H0 <- NA
  } 
  else {
    M1_outWALD_H0 <- linearHypothesis(M1_modelWALD_H0, M1_hypothesis)
    M1_S_obsWALD_H0 <- M1_outWALD_H0$F[2] # observed statistics under H0
    M1_pvaluesWALD_H0 <- M1_outWALD_H0$`Pr(>F)`[2]
    M1_rejectsWALD <- ifelse(M1_pvaluesWALD_H0<=0.05,1,0)
    
    # Permutation test (H0)
    M1_outPERM_H0 <- perm_test_pval_EXT(dataset = data_strPBR, stat_obs = M1_S_obsWALD_H0, 
                                    hyp = 0, L = perms, analysis = 0)
    M1_pvaluesPERM_H0 <- M1_outPERM_H0$pvalue
    M1_rejectsPERM <- ifelse(M1_pvaluesPERM_H0<=0.05,1,0)
    M1_nulls_H0 <- M1_outPERM_H0$nulls
  }
  
  # H1
  M1_modelWALD_H1 <- lm(as.formula(paste("y_H1", M1_formula_str)), data = M1_data_strPBR_dummy)
  if(any(is.na(coef(M1_modelWALD_H1)))) {
    M1_pvaluesWALD_H1 <- NA
    M1_notrejectsWALD <- NA
    M1_pvaluesPERM_H1 <- NA
    M1_notrejectsPERM <- NA
    M1_nulls_H1 <- NA
  } 
  else {
    M1_outWALD_H1 <- linearHypothesis(M1_modelWALD_H1, M1_hypothesis)
    M1_S_obsWALD_H1 <- M1_outWALD_H1$F[2] # observed statistics under H1
    M1_pvaluesWALD_H1 <- M1_outWALD_H1$`Pr(>F)`[2]
    M1_notrejectsWALD <- ifelse(M1_pvaluesWALD_H1>0.05,1,0)
    
    # Permutation test (H1)
    M1_outPERM_H1 <- perm_test_pval_EXT(dataset = data_strPBR, stat_obs = M1_S_obsWALD_H1, 
                                    hyp = 1, L = perms, analysis = 0)
    M1_pvaluesPERM_H1 <- M1_outPERM_H1$pvalue
    M1_notrejectsPERM <- ifelse(M1_pvaluesPERM_H1>0.05,1,0)
    M1_nulls_H1 <- M1_outPERM_H1$nulls
  }
  
  
  # F-specification 2
  
  # F-test
  M2_analysis <- analysis2(dataset = data_strPBR)
  M2_data_strPBR_dummy <- M2_analysis$db
  M2_formula_str <- M2_analysis$formula
  M2_hypothesis <- M2_analysis$hypoth
  
  # H0
  M2_modelWALD_H0 <- lm(as.formula(paste("y_H0", M2_formula_str)), data = M2_data_strPBR_dummy)
  if(any(is.na(coef(M2_modelWALD_H0)))) {
    M2_pvaluesWALD_H0 <- NA
    M2_rejectsWALD <- NA
    M2_pvaluesPERM_H0 <- NA
    M2_rejectsPERM <- NA
    M2_nulls_H0 <- NA
  } 
  else {
    M2_outWALD_H0 <- linearHypothesis(M2_modelWALD_H0, M2_hypothesis)
    M2_S_obsWALD_H0 <- M2_outWALD_H0$F[2] # observed statistics under H0
    M2_pvaluesWALD_H0 <- M2_outWALD_H0$`Pr(>F)`[2]
    M2_rejectsWALD <- ifelse(M2_pvaluesWALD_H0<=0.05,1,0)
    
    # Permutation test (H0)
    M2_outPERM_H0 <- perm_test_pval_EXT(dataset = data_strPBR, stat_obs = M2_S_obsWALD_H0, 
                                    hyp = 0, L = perms, analysis = 1)
    M2_pvaluesPERM_H0 <- M2_outPERM_H0$pvalue
    M2_rejectsPERM <- ifelse(M2_pvaluesPERM_H0<=0.05,1,0)
    M2_nulls_H0 <- M2_outPERM_H0$nulls
  }
  
  # H1
  M2_modelWALD_H1 <- lm(as.formula(paste("y_H1", M2_formula_str)), data = M2_data_strPBR_dummy)
  if(any(is.na(coef(M2_modelWALD_H1)))) {
    M2_pvaluesWALD_H1 <- NA
    M2_notrejectsWALD <- NA
    M2_pvaluesPERM_H1 <- NA
    M2_notrejectsPERM <- NA
    M2_nulls_H1 <- NA
  } 
  else {
    M2_outWALD_H1 <- linearHypothesis(M2_modelWALD_H1, M2_hypothesis)
    M2_S_obsWALD_H1 <- M2_outWALD_H1$F[2] # observed statistics under H1
    M2_pvaluesWALD_H1 <- M2_outWALD_H1$`Pr(>F)`[2]
    M2_notrejectsWALD <- ifelse(M2_pvaluesWALD_H1>0.05,1,0)
    
    # Permutation test (H1)
    M2_outPERM_H1 <- perm_test_pval_EXT(dataset = data_strPBR, stat_obs = M2_S_obsWALD_H1, 
                                    hyp = 1, L = perms, analysis = 1)
    M2_pvaluesPERM_H1 <- M2_outPERM_H1$pvalue
    M2_notrejectsPERM <- ifelse(M2_pvaluesPERM_H1>0.05,1,0)
    M2_nulls_H1 <- M2_outPERM_H1$nulls
  }
  
  
  # F-specification 3
  
  # F-test
  M3_analysis <- analysis3(dataset = data_strPBR)
  M3_data_strPBR_dummy <- M3_analysis$db
  M3_formula_str <- M3_analysis$formula
  M3_hypothesis <- M3_analysis$hypoth
  
  # H0
  M3_modelWALD_H0 <- lm(as.formula(paste("y_H0", M3_formula_str)), data = M3_data_strPBR_dummy)
  if(any(is.na(coef(M3_modelWALD_H0)))) {
    M3_pvaluesWALD_H0 <- NA
    M3_rejectsWALD <- NA
    M3_pvaluesPERM_H0 <- NA
    M3_rejectsPERM <- NA
    M3_nulls_H0 <- NA
  } 
  else {
    M3_outWALD_H0 <- linearHypothesis(M3_modelWALD_H0, M3_hypothesis)
    M3_S_obsWALD_H0 <- M3_outWALD_H0$F[2] # observed statistics under H0
    M3_pvaluesWALD_H0 <- M3_outWALD_H0$`Pr(>F)`[2]
    M3_rejectsWALD <- ifelse(M3_pvaluesWALD_H0<=0.05,1,0)
    
    # Permutation test (H0)
    M3_outPERM_H0 <- perm_test_pval_EXT(dataset = data_strPBR, stat_obs = M3_S_obsWALD_H0, 
                                        hyp = 0, L = perms, analysis = 2)
    M3_pvaluesPERM_H0 <- M3_outPERM_H0$pvalue
    M3_rejectsPERM <- ifelse(M3_pvaluesPERM_H0<=0.05,1,0)
    M3_nulls_H0 <- M3_outPERM_H0$nulls
  }
  
  # H1
  M3_modelWALD_H1 <- lm(as.formula(paste("y_H1", M3_formula_str)), data = M3_data_strPBR_dummy)
  if(any(is.na(coef(M3_modelWALD_H1)))) {
    M3_pvaluesWALD_H1 <- NA
    M3_notrejectsWALD <- NA
    M3_pvaluesPERM_H1 <- NA
    M3_notrejectsPERM <- NA
    M3_nulls_H1 <- NA
  } 
  else {
    M3_outWALD_H1 <- linearHypothesis(M3_modelWALD_H1, M3_hypothesis)
    M3_S_obsWALD_H1 <- M3_outWALD_H1$F[2] # observed statistics under H1
    M3_pvaluesWALD_H1 <- M3_outWALD_H1$`Pr(>F)`[2]
    M3_notrejectsWALD <- ifelse(M3_pvaluesWALD_H1>0.05,1,0)
    
    # Permutation test (H1)
    M3_outPERM_H1 <- perm_test_pval_EXT(dataset = data_strPBR, stat_obs = M3_S_obsWALD_H1, 
                                        hyp = 1, L = perms, analysis = 2)
    M3_pvaluesPERM_H1 <- M3_outPERM_H1$pvalue
    M3_notrejectsPERM <- ifelse(M3_pvaluesPERM_H1>0.05,1,0)
    M3_nulls_H1 <- M3_outPERM_H1$nulls
  }
  
  
  tibble::tibble(pvaluesNAIVET_H0, rejectsNAIVET, pvaluesNAIVET_H1, notrejectsNAIVET,
                 M1_pvaluesWALD_H0, M1_rejectsWALD, M1_pvaluesWALD_H1, M1_notrejectsWALD,
                 M1_pvaluesPERM_H0, M1_rejectsPERM, M1_nulls_H0, 
                 M1_pvaluesPERM_H1, M1_notrejectsPERM, M1_nulls_H1,
                 M2_pvaluesWALD_H0, M2_rejectsWALD, M2_pvaluesWALD_H1, M2_notrejectsWALD,
                 M2_pvaluesPERM_H0, M2_rejectsPERM, M2_nulls_H0, 
                 M2_pvaluesPERM_H1, M2_notrejectsPERM, M2_nulls_H1,
                 M3_pvaluesWALD_H0, M3_rejectsWALD, M3_pvaluesWALD_H1, M3_notrejectsWALD,
                 M3_pvaluesPERM_H0, M3_rejectsPERM, M3_nulls_H0, 
                 M3_pvaluesPERM_H1, M3_notrejectsPERM, M3_nulls_H1
  )
  
}

## Scenario 2 (interim analysis) and exclusion of patients - first simulation:
## the least effective experimental treatment is droppedafter an interim futility analysis.
## Arguments:
## * 'N' is the total number of planned patients;
## * 'M' is the number of molecular alterations considered at the beginning of the trial;
## * 'prevalence' is the vector of molecular alterations' prevalences;
## * 'trt_labels' is the vector of treatment labels (first entry for control treatment);
## * 'probs' is the vector of treatment probabilities;
## * 'bs' is the block size to be used for block randomization;
## * 'alpha' is the control effect;
## * 'sigma' is the vector of prognoses;
## * 'tau_H0' and 'tau_H1' are the vectors of treatment effects under the null (H0) and
## alternative (H1) hypothesis, respectively;
## * 'z_c' is the pre-planned futility threshold (can be Inf if the interim futility
## analysis is not planned);
## * 'perms' is the number of permutations;
## * 'max_atts' is the maximum number of attempts to compute the permutation test
## (it can be Inf).
simulate_one3 <- function(N = 200, M = 4, prevalence = c(1/4,1/4,1/4,1/4), trt_labels = c(1,2),
                          probs = c(0.5,0.5), bs = 4, alpha = 1, sigma = c(0.5,0.5,0.5,0.5), 
                          tau_H0 = c(0,0,0,0), tau_H1 = c(0.40,0.40,0.40,-0.40),
                          z_c = 0, perms = 1000, max_atts = 100000)
  {
  
  # controls
  if (length(N)!=1 || !is.numeric(N) || N <= 0 || N != round(N)) 
    stop("'N' must be a positive integer")
  if (length(M)!=1 || !is.numeric(M) || M <= 0 || M != round(M)) 
    stop("'M' must be a positive integer")
  if (length(alpha)!=1 || !is.numeric(alpha)) stop("'alpha' must be numeric")
  if (length(sigma) != M || length(tau_H0) != M || length(tau_H1) != M) 
    stop("'sigma', 'tau_H0', and 'tau_H1' must have length equal to 'M'")
  if (!is.numeric(sigma)) 
    stop("All entries of 'sigma' must be numeric")
  if (!is.numeric(tau_H0)) 
    stop("All entries of 'tau_H0' must be numeric")
  if (!is.numeric(tau_H1)) 
    stop("All entries of 'tau_H1' must be numeric")
  if (length(z_c)!=1 || !is.numeric(z_c) || is.nan(z_c)) 
    stop("'z_c' must be numeric or infinit")
  
  # STAGE I
  
  #  data generation
  molec_alter1 <- sample(1:M, prob = prevalence, size = N/2, replace = T)
  stage1 <- rep(1,N/2)
  data1 <- data.frame(stratum = molec_alter1, stage = stage1)
  N1 <- nrow(data1)
  data1$ID <- 1:N1
  
  #  treatment assignment
  data_strPBR1 <- my_strat_block_rand(dataset = data1, trt_seq = trt_labels, 
                                      prob_each = probs, blocksize = bs)
  data_strPBR1 <- data_strPBR1[order(data_strPBR1$ID), ]
  
  # continuous outcome generation
  data_strPBR1$y_H0 <- rep(0,N1)
  data_strPBR1$y_H1 <- rep(0,N1)
  for(j in 1:N1){
    data_strPBR1$err[j] <- rnorm(1,0,1)
    data_strPBR1$y_H0[j] <- alpha + sigma[data_strPBR1$stratum[j]] + 
      tau_H0[data_strPBR1$stratum[j]]*(data_strPBR1$trt[j]-1) + data_strPBR1$err[j]
    data_strPBR1$y_H1[j] <- alpha + sigma[data_strPBR1$stratum[j]] + 
      tau_H1[data_strPBR1$stratum[j]]*(data_strPBR1$trt[j]-1) + data_strPBR1$err[j]
  }
  
  ctrl_label <- trt_labels[1]
  trt_label <- trt_labels[2]
  stddiff_H0 <- sapply(1:M, function(s) {
    stratum_ctrl <- data_strPBR1[data_strPBR1$stratum == s & data_strPBR1$trt == ctrl_label, ]
    stratum_trt <- data_strPBR1[data_strPBR1$stratum == s & data_strPBR1$trt == trt_label, ]
    mean_diff <- mean(stratum_trt$y_H0) - mean(stratum_ctrl$y_H0)
    sd_ctrl <- sd(stratum_ctrl$y_H0)
    return(mean_diff / sd_ctrl)
    # if the number of subjects and variability differ between the treatment and control groups, 
    # comment the two previous lines and uncomment the following: 
    # sd_pooled <- sqrt(((nrow(stratum_trt)-1)*var(stratum_trt$y_H0) +
    #                    (nrow(stratum_ctrl)-1)*var(stratum_ctrl$y_H0))/
    #                    (nrow(stratum_trt)+nrow(stratum_ctrl)-2))
    # denom <- sd_pooled*(1/nrow(stratum_trt)+1/nrow(stratum_ctrl))
    # return(mean_diff / denom)
  })
  stddiff_H1 <- sapply(1:M, function(s) {
    stratum_ctrl <- data_strPBR1[data_strPBR1$stratum == s & data_strPBR1$trt == ctrl_label, ]
    stratum_trt <- data_strPBR1[data_strPBR1$stratum == s & data_strPBR1$trt == trt_label, ]
    mean_diff <- mean(stratum_trt$y_H1) - mean(stratum_ctrl$y_H1)
    sd_ctrl <- sd(stratum_ctrl$y_H1)
    return(mean_diff / sd_ctrl)
    # if the number of subjects and variability differ between the treatment and control groups, 
    # comment the two previous lines and uncomment the following: 
    # sd_pooled <- sqrt(((nrow(stratum_trt)-1)*var(stratum_trt$y_H1) +
    #                    (nrow(stratum_ctrl)-1)*var(stratum_ctrl$y_H1))/
    #                    (nrow(stratum_trt)+nrow(stratum_ctrl)-2))
    # denom <- sd_pooled*(1/nrow(stratum_trt)+1/nrow(stratum_ctrl))
    # return(mean_diff / denom)
  })
  
  min_diff_H0 <- min(stddiff_H0)
  drops_H0 <- ifelse(stddiff_H0 == min_diff_H0 & stddiff_H0 <= z_c, 1, 0)
  min_diff_H1 <- min(stddiff_H1)
  drops_H1 <- ifelse(stddiff_H1 == min_diff_H1 & stddiff_H1 <= z_c, 1, 0)
  
  pooled_stratum_H0 <- which(drops_H0 == 0)
  pooled_stratum_H1 <- which(drops_H1 == 0)
  
  drop_trt_H0 <- ifelse(length(pooled_stratum_H0) == M, 0, setdiff(1:M, pooled_stratum_H0))
  drop_trt_H1 <- ifelse(length(pooled_stratum_H1) == M, 0, setdiff(1:M, pooled_stratum_H1))
  
  if (drop_trt_H0 == 0) {
    new_prevalence_H0 <- prevalence
  } else {
    new_prevalence_H0 <- prevalence[-drop_trt_H0]
  }
  
  if (drop_trt_H1 == 0) {
    new_prevalence_H1 <- prevalence
  } else {
    new_prevalence_H1 <- prevalence[-drop_trt_H1]
  }
  
  
  # STAGE II
  
  # H0
  
  #  data generation
  data_strPBR1_kept_H0 <- data_strPBR1[data_strPBR1$stratum!=drop_trt_H0,]
  N1_k_H0 <- nrow(data_strPBR1_kept_H0)
  N2_H0 <- N - N1_k_H0  
  molec_alter2_H0 <- sample(setdiff(1:M, drop_trt_H0), prob = new_prevalence_H0, 
                            N2_H0, replace = T)
  stage2_H0 <- rep(2, N2_H0)
  data2_H0 <- data.frame(stratum = molec_alter2_H0, stage = stage2_H0)
  data2_H0$ID <- 1:N2_H0
  
  #  treatment assignment
  data_strPBR2_H0 <- my_strat_block_rand(dataset = data2_H0, trt_seq = trt_labels, 
                                         prob_each = probs, blocksize = bs)
  data_strPBR2_H0 <- data_strPBR2_H0[order(data_strPBR2_H0$ID), ]
  
  # outcome generation
  data_strPBR2_H0$y_H0 <- rep(0,N2_H0)
  data_strPBR2_H0$y_H1 <- rep(0,N2_H0)
  for(j in 1:N2_H0){
    data_strPBR2_H0$err[j] <- rnorm(1,0,1)
    data_strPBR2_H0$y_H0[j] <- alpha + sigma[data_strPBR2_H0$stratum[j]] + 
      tau_H0[data_strPBR2_H0$stratum[j]]*(data_strPBR2_H0$trt[j]-1) + data_strPBR2_H0$err[j]
    data_strPBR2_H0$y_H1[j] <- alpha + sigma[data_strPBR2_H0$stratum[j]] + 
      tau_H1[data_strPBR2_H0$stratum[j]]*(data_strPBR2_H0$trt[j]-1) + data_strPBR2_H0$err[j]
  }
  
  data_strPBR_H0 <- rbind(data_strPBR1_kept_H0,data_strPBR2_H0)
  Ntot_H0 <- nrow(data_strPBR_H0)
  
  
  # H1
  
  #  data generation
  data_strPBR1_kept_H1 <- data_strPBR1[data_strPBR1$stratum!=drop_trt_H1,]
  N1_k_H1 <- nrow(data_strPBR1_kept_H1)
  N2_H1 <- N - N1_k_H1  
  molec_alter2_H1 <- sample(setdiff(1:M, drop_trt_H1), prob = new_prevalence_H1, 
                            N2_H1, replace = T)
  stage2_H1 <- rep(2, N2_H1)
  data2_H1 <- data.frame(stratum = molec_alter2_H1, stage = stage2_H1)
  data2_H1$ID <- 1:N2_H1
  
  #  treatment assignment
  data_strPBR2_H1 <- my_strat_block_rand(dataset = data2_H1, trt_seq = trt_labels, 
                                         prob_each = probs, blocksize = bs)
  data_strPBR2_H1 <- data_strPBR2_H1[order(data_strPBR2_H1$ID), ]
  
  # outcome generation
  data_strPBR2_H1$y_H0 <- rep(0,N2_H1)
  data_strPBR2_H1$y_H1 <- rep(0,N2_H1)
  for(j in 1:N2_H1){
    data_strPBR2_H1$err[j] <- rnorm(1,0,1)
    data_strPBR2_H1$y_H0[j] <- alpha + sigma[data_strPBR2_H1$stratum[j]] + 
      tau_H0[data_strPBR2_H1$stratum[j]]*(data_strPBR2_H1$trt[j]-1) + data_strPBR2_H1$err[j]
    data_strPBR2_H1$y_H1[j] <- alpha + sigma[data_strPBR2_H1$stratum[j]] + 
      tau_H1[data_strPBR2_H1$stratum[j]]*(data_strPBR2_H1$trt[j]-1) + data_strPBR2_H1$err[j]
  }
  
  data_strPBR_H1 <- rbind(data_strPBR1_kept_H1,data_strPBR2_H1)
  Ntot_H1 <- nrow(data_strPBR_H1)
  
  
  # t-test (NAIVE)
  # H0
  outNAIVET_H0 <- t.test(y_H0 ~ factor(trt), data = data_strPBR_H0)
  pvaluesNAIVET_H0 <- outNAIVET_H0$p.value
  rejectsNAIVET <- ifelse(pvaluesNAIVET_H0<=0.05,1,0)
  
  # H1
  outNAIVET_H1 <- t.test(y_H1 ~ factor(trt), data = data_strPBR_H1)
  pvaluesNAIVET_H1 <- outNAIVET_H1$p.value
  notrejectsNAIVET <- ifelse(pvaluesNAIVET_H1>0.05,1,0)
  
  
  # F-test (F-specification 1)
  
  # H0
  analysis_H0 <- analysis1(dataset = data_strPBR_H0)
  data_strPBR_dummy_H0 <- analysis_H0$db
  formula_str_H0 <- analysis_H0$formula
  hypothesis_H0 <- analysis_H0$hypoth
  
  modelWALD_H0 <- lm(as.formula(paste("y_H0", formula_str_H0)), data = data_strPBR_dummy_H0)
  if(any(is.na(coef(modelWALD_H0)))) {
    pvaluesWALD_H0 <- NA
    rejectsWALD <- NA
    pvaluesPERM_H0 <- NA
    rejectsPERM <- NA
    nulls_H0 <- NA
    attempts_H0 <- NA
  } 
  else {
    outWALD_H0 <- linearHypothesis(modelWALD_H0, hypothesis_H0)
    S_obsWALD_H0 <- outWALD_H0$F[2] # observed statistics under H0
    pvaluesWALD_H0 <- outWALD_H0$`Pr(>F)`[2]
    rejectsWALD <- ifelse(pvaluesWALD_H0<=0.05,1,0)
    
    # Permutation test (H0)
    outPERM_H0 <- perm_test_pval_IA(dataset1 = data_strPBR1, dataset2 = data_strPBR2_H0, crit_val = z_c, drop_trt = drop_trt_H0, stat_obs = S_obsWALD_H0, hyp = 0,
                                    L = perms, max_attempts = max_atts)
    pvaluesPERM_H0 <- outPERM_H0$pvalue
    rejectsPERM <- ifelse(pvaluesPERM_H0<=0.05,1,0)
    nulls_H0 <- outPERM_H0$nulls
    attempts_H0 <- outPERM_H0$total_attempts
  }
  
  # H1
  analysis_H1 <- analysis1(dataset = data_strPBR_H1)
  data_strPBR_dummy_H1 <- analysis_H1$db
  formula_str_H1 <- analysis_H1$formula
  hypothesis_H1 <- analysis_H1$hypoth
  
  modelWALD_H1 <- lm(as.formula(paste("y_H1", formula_str_H1)), data = data_strPBR_dummy_H1)
  if(any(is.na(coef(modelWALD_H1)))) {
    pvaluesWALD_H1 <- NA
    notrejectsWALD <- NA
    pvaluesPERM_H1 <- NA
    notrejectsPERM <- NA
    nulls_H1 <- NA
    attempts_H1 <- NA
  } 
  else {
    outWALD_H1 <- linearHypothesis(modelWALD_H1, hypothesis_H1)
    S_obsWALD_H1 <- outWALD_H1$F[2] # observed statistics under H1
    pvaluesWALD_H1 <- outWALD_H1$`Pr(>F)`[2]
    notrejectsWALD <- ifelse(pvaluesWALD_H1>0.05,1,0)
    
    # Permutation test (H1)
    outPERM_H1 <- perm_test_pval_IA(dataset1 = data_strPBR1, dataset2 = data_strPBR2_H1, crit_val = z_c, drop_trt = drop_trt_H1, stat_obs = S_obsWALD_H1, hyp = 1,
                                    L = perms, max_attempts = max_atts)
    pvaluesPERM_H1 <- outPERM_H1$pvalue
    notrejectsPERM <- ifelse(pvaluesPERM_H1>0.05,1,0)
    nulls_H1 <- outPERM_H1$nulls
    attempts_H1 <- outPERM_H1$total_attempts
  }
  
  
  tibble::tibble(drop1_H0=drops_H0[1],drop2_H0=drops_H0[2],drop3_H0=drops_H0[3],drop4_H0=drops_H0[4],
                 drop1_H1=drops_H1[1],drop2_H1=drops_H1[2],drop3_H1=drops_H1[3],drop4_H1=drops_H1[4],
                 pvaluesNAIVET_H0, rejectsNAIVET, pvaluesNAIVET_H1, notrejectsNAIVET,
                 pvaluesWALD_H0, rejectsWALD, pvaluesWALD_H1, notrejectsWALD,
                 pvaluesPERM_H0, rejectsPERM, nulls_H0, attempts_H0, 
                 pvaluesPERM_H1, notrejectsPERM, nulls_H1, attempts_H1
  )
  
}

## Scenario 2 (interim analysis) and orientation of patients - first simulation:
## the least effective experimental treatment is droppedafter an interim futility analysis.
## Arguments:
## * 'N' is the total number of planned patients;
## * 'M' is the number of molecular alterations considered at the beginning of the trial;
## * 'prevalence' is the vector of molecular alterations' prevalences;
## * 'trt_labels' is the vector of treatment labels (first entry for control treatment);
## * 'probs' is the vector of treatment probabilities;
## * 'bs' is the block size to be used for block randomization;
## * 'alpha' is the control effect;
## * 'sigma' is the vector of prognoses;
## * 'tau_H0' and 'tau_H1' are the matrix of treatment effects under the null (H0) and
## alternative (H1) hypothesis, respectively;
## * 'z_c' is the pre-planned futility threshold (can be Inf if the interim futility
## analysis is not planned);
## * 'probs_orient' is the vector of probabilities of assigning a certain stratum;
## * 'perms' is the number of permutations;
## * 'max_atts' is the maximum number of attempts to compute the permutation test
## (it can be Inf).
simulate_one4 <- function(N = 200, M = 4, prevalence = c(1/4,1/4,1/4,1/4), trt_labels = c(1,2), 
                          probs = c(0.5,0.5), bs = 4, alpha = 1,  sigma = c(0.5,0.5,0.5,0.5), 
                          tau_H0 = matrix(c(0,0,0,0,
                                            0,0,0,0,
                                            0,0,0,0,
                                            0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE), 
                          tau_H1 = matrix(c(0.40,0.40,0.40,-0.40,
                                            0.40,0.40,0.40,-0.40,
                                            0.40,0.40,0.40,-0.40,
                                            0.40,0.40,0.40,-0.40), nrow = 4, ncol = 4, byrow = TRUE), 
                          z_c = 0, probs_orient = c(1/3,1/3,1/3), perms = 1000, max_atts = 100000)
  {
  
  # controls
  if (length(N)!=1 || !is.numeric(N) || N <= 0 || N != round(N)) 
    stop("'N' must be a positive integer")
  if (length(M)!=1 || !is.numeric(M) || M <= 0 || M != round(M)) 
    stop("'M' must be a positive integer")
  if (length(alpha)!=1 || !is.numeric(alpha)) 
    stop("'alpha' must be numeric")
  if (length(sigma) != M) stop("'sigma' must have length equal to 'M'")
  if (nrow(tau_H0) != M || ncol(tau_H0) != M) 
    stop("'tau_H0' must have dimention equal to 'MxM'")
  if (nrow(tau_H1) != M || ncol(tau_H1) != M) 
    stop("'tau_H1' must have dimention equal to 'MxM'")
  if (!is.numeric(sigma)) 
    stop("All entries of 'sigma' must be numeric")
  if (!is.numeric(tau_H0)) 
    stop("All entries of 'tau_H0' must be numeric")
  if (!is.numeric(tau_H1)) 
    stop("All entries of 'tau_H1' must be numeric")
  if (length(z_c)!=1 || !is.numeric(z_c) || is.nan(z_c)) 
    stop("'z_c' must be numeric or infinit")
  if (length(probs_orient)!=M-1) 
    stop("'probs_orient' must have length equal to 'M'-1")
  
  # STAGE I
  
  #  data generation
  molec_alter1 <- sample(1:4, prob = prevalence, size = N/2, replace = T)
  stage1 <- rep(1,N/2)
  data1 <- data.frame(stratum_old = molec_alter1, stratum = molec_alter1, stage = stage1)
  N1 <- nrow(data1)
  data1$ID <- 1:N1
  
  #  treatment assignment
  data_strPBR1 <- my_strat_block_rand(dataset = data1, trt_seq = trt_labels, 
                                      prob_each = probs, blocksize = bs)
  data_strPBR1 <- data_strPBR1[order(data_strPBR1$ID), ]
  
  # continuous outcome generation
  data_strPBR1$y_H0 <- rep(0,N1)
  data_strPBR1$y_H1 <- rep(0,N1)
  for(j in 1:N1){
    data_strPBR1$err[j] <- rnorm(1,0,1)
    data_strPBR1$y_H0[j] <- alpha + sigma[data_strPBR1$stratum_old[j]] + 
      tau_H0[data_strPBR1$stratum_old[j],data_strPBR1$stratum[j]]*(data_strPBR1$trt[j]-1) + 
      data_strPBR1$err[j]
    data_strPBR1$y_H1[j] <- alpha + sigma[data_strPBR1$stratum_old[j]] + 
      tau_H1[data_strPBR1$stratum_old[j],data_strPBR1$stratum[j]]*(data_strPBR1$trt[j]-1) + 
      data_strPBR1$err[j]
  }
  
  ctrl_label <- trt_labels[1]
  trt_label <- trt_labels[2]
  stddiff_H0 <- sapply(1:M, function(s) {
    stratum_ctrl <- data_strPBR1[data_strPBR1$stratum == s & 
                                   data_strPBR1$trt == ctrl_label, ]
    stratum_trt <- data_strPBR1[data_strPBR1$stratum == s & 
                                  data_strPBR1$trt == trt_label, ]
    mean_diff <- mean(stratum_trt$y_H0) - mean(stratum_ctrl$y_H0)
    sd_ctrl <- sd(stratum_ctrl$y_H0)
    return(mean_diff / sd_ctrl)
    # if the number of subjects and variability differ between the treatment and control groups, 
    # comment the two previous lines and uncomment the following: 
    # sd_pooled <- sqrt(((nrow(stratum_trt)-1)*var(stratum_trt$y_H0) +
    #                    (nrow(stratum_ctrl)-1)*var(stratum_ctrl$y_H0))/
    #                    (nrow(stratum_trt)+nrow(stratum_ctrl)-2))
    # denom <- sd_pooled*(1/nrow(stratum_trt)+1/nrow(stratum_ctrl))
    # return(mean_diff / denom)
  })
  stddiff_H1 <- sapply(1:M, function(s) {
    stratum_ctrl <- data_strPBR1[data_strPBR1$stratum == s & 
                                   data_strPBR1$trt == ctrl_label, ]
    stratum_trt <- data_strPBR1[data_strPBR1$stratum == s & 
                                  data_strPBR1$trt == trt_label, ]
    mean_diff <- mean(stratum_trt$y_H1) - mean(stratum_ctrl$y_H1)
    sd_ctrl <- sd(stratum_ctrl$y_H1)
    return(mean_diff / sd_ctrl)
    # if the number of subjects and variability differ between the treatment and control groups, 
    # comment the two previous lines and uncomment the following: 
    # sd_pooled <- sqrt(((nrow(stratum_trt)-1)*var(stratum_trt$y_H1) +
    #                    (nrow(stratum_ctrl)-1)*var(stratum_ctrl$y_H1))/
    #                    (nrow(stratum_trt)+nrow(stratum_ctrl)-2))
    # denom <- sd_pooled*(1/nrow(stratum_trt)+1/nrow(stratum_ctrl))
    # return(mean_diff / denom)
  })
  
  min_diff_H0 <- min(stddiff_H0)
  drops_H0 <- ifelse(stddiff_H0 == min_diff_H0 & stddiff_H0 <= z_c, 1, 0)
  min_diff_H1 <- min(stddiff_H1)
  drops_H1 <- ifelse(stddiff_H1 == min_diff_H1 & stddiff_H1 <= z_c, 1, 0)
  
  pooled_stratum_H0 <- which(drops_H0 == 0)
  pooled_stratum_H1 <- which(drops_H1 == 0)
  
  drop_trt_H0 <- ifelse(length(pooled_stratum_H0) == M, 0, setdiff(1:M, pooled_stratum_H0))
  drop_trt_H1 <- ifelse(length(pooled_stratum_H1) == M, 0, setdiff(1:M, pooled_stratum_H1))
  
  if (drop_trt_H0 == 0) {
    new_prevalence_H0 <- prevalence
  } else {
    new_prevalence_H0 <- prevalence[-drop_trt_H0]
  }
  
  if (drop_trt_H1 == 0) {
    new_prevalence_H1 <- prevalence
  } else {
    new_prevalence_H1 <- prevalence[-drop_trt_H1]
  }
  
  
  # STAGE II
  
  # H0
  
  #  data generation
  data_strPBR1_kept_H0 <- data_strPBR1[data_strPBR1$stratum!=drop_trt_H0,]
  N1_k_H0 <- nrow(data_strPBR1_kept_H0)
  N2_H0 <- N - N1_k_H0  
  molec_alter2_H0 <- sample(1:M, prob = prevalence, N2_H0, replace = T)
  stage2_H0 <- rep(2, N2_H0)
  data2_H0 <- data.frame(stratum_old = molec_alter2_H0, stage = stage2_H0)
  data2_H0$ID <- 1:N2_H0
  
  if (drop_trt_H0 != 0) {
    # orientation of patients belonging to the dropped stratum (block randomization)
    data2_maintain_H0 <- data2_H0[data2_H0$stratum_old!=drop_trt_H0,]
    data2_maintain_H0$stratum <- data2_maintain_H0$stratum_old
    data2_toorient_H0 <- data2_H0[data2_H0$stratum_old==drop_trt_H0,]
    data2_oriented_H0 <- my_block_rand(data2_toorient_H0, strata_seq = setdiff(1:M, drop_trt_H0), 
                                       prob_each_orient = probs_orient)
    data2_toassign_H0 <- rbind(data2_maintain_H0,data2_oriented_H0)
    data2_toassign_H0 <- data2_toassign_H0[order(data2_toassign_H0$ID), ]
  }
  else {
    data2_toassign_H0 <- data2_H0
    data2_toassign_H0$stratum <- data2_toassign_H0$stratum_old
  }
  
  #  treatment assignment
  data_strPBR2_H0 <- my_strat_block_rand(dataset = data2_toassign_H0, trt_seq = trt_labels, 
                                         prob_each = probs, blocksize = bs)
  data_strPBR2_H0 <- data_strPBR2_H0[order(data_strPBR2_H0$ID), ]
  
  # outcome generation
  data_strPBR2_H0$y_H0 <- rep(0,N2_H0)
  data_strPBR2_H0$y_H1 <- rep(0,N2_H0)
  for(j in 1:N2_H0){
    data_strPBR2_H0$err[j] <- rnorm(1,0,1)
    data_strPBR2_H0$y_H0[j] <- alpha + sigma[data_strPBR2_H0$stratum_old[j]] + 
      tau_H0[data_strPBR2_H0$stratum_old[j],data_strPBR2_H0$stratum[j]]*(data_strPBR2_H0$trt[j]-1) + 
      data_strPBR2_H0$err[j]
    data_strPBR2_H0$y_H1[j] <- alpha + sigma[data_strPBR2_H0$stratum_old[j]] + 
      tau_H1[data_strPBR2_H0$stratum_old[j],data_strPBR2_H0$stratum[j]]*(data_strPBR2_H0$trt[j]-1) + 
      data_strPBR2_H0$err[j]
  }
  
  data_strPBR_H0 <- rbind(data_strPBR1_kept_H0,data_strPBR2_H0)
  Ntot_H0 <- nrow(data_strPBR_H0)
  
  
  # H1
  
  #  data generation
  data_strPBR1_kept_H1 <- data_strPBR1[data_strPBR1$stratum!=drop_trt_H1,]
  N1_k_H1 <- nrow(data_strPBR1_kept_H1)
  N2_H1 <- N - N1_k_H1  
  molec_alter2_H1 <- sample(1:M, prob = prevalence, N2_H1, replace = T)
  stage2_H1 <- rep(2, N2_H1)
  data2_H1 <- data.frame(stratum_old = molec_alter2_H1, stage = stage2_H1)
  data2_H1$ID <- 1:N2_H1
  
  if (drop_trt_H1 != 0) {
    # orientation of patients belonging to the dropped stratum (block randomization)
    data2_maintain_H1 <- data2_H1[data2_H1$stratum_old!=drop_trt_H1,]
    data2_maintain_H1$stratum <- data2_maintain_H1$stratum_old
    data2_toorient_H1 <- data2_H1[data2_H1$stratum_old==drop_trt_H1,]
    data2_oriented_H1 <- my_block_rand(data2_toorient_H1, strata_seq = setdiff(1:M, drop_trt_H1), 
                                       prob_each_orient = probs_orient)
    data2_toassign_H1 <- rbind(data2_maintain_H1,data2_oriented_H1)
    data2_toassign_H1 <- data2_toassign_H1[order(data2_toassign_H1$ID), ]
  }
  else {
    data2_toassign_H1 <- data2_H1
    data2_toassign_H1$stratum <- data2_toassign_H1$stratum_old
  }
  
  #  treatment assignment
  data_strPBR2_H1 <- my_strat_block_rand(dataset = data2_toassign_H1, trt_seq = trt_labels, 
                                         prob_each = probs, blocksize = bs)
  data_strPBR2_H1 <- data_strPBR2_H1[order(data_strPBR2_H1$ID), ]
  
  # outcome generation
  data_strPBR2_H1$y_H0 <- rep(0,N2_H1)
  data_strPBR2_H1$y_H1 <- rep(0,N2_H1)
  for(j in 1:N2_H1){
    data_strPBR2_H1$err[j] <- rnorm(1,0,1)
    data_strPBR2_H1$y_H0[j] <- alpha + sigma[data_strPBR2_H1$stratum_old[j]] + 
      tau_H0[data_strPBR2_H1$stratum_old[j],data_strPBR2_H1$stratum[j]]*(data_strPBR2_H1$trt[j]-1) + 
      data_strPBR2_H1$err[j]
    data_strPBR2_H1$y_H1[j] <- alpha + sigma[data_strPBR2_H1$stratum_old[j]] + 
      tau_H1[data_strPBR2_H1$stratum_old[j],data_strPBR2_H1$stratum[j]]*(data_strPBR2_H1$trt[j]-1) + 
      data_strPBR2_H1$err[j]
  }
  
  data_strPBR_H1 <- rbind(data_strPBR1_kept_H1,data_strPBR2_H1)
  Ntot_H1 <- nrow(data_strPBR_H1)
  
  
  # t-test (NAIVE)
  # H0
  outNAIVET_H0 <- t.test(y_H0 ~ factor(trt), data = data_strPBR_H0)
  pvaluesNAIVET_H0 <- outNAIVET_H0$p.value
  rejectsNAIVET <- ifelse(pvaluesNAIVET_H0<=0.05,1,0)
  
  # H1
  outNAIVET_H1 <- t.test(y_H1 ~ factor(trt), data = data_strPBR_H1)
  pvaluesNAIVET_H1 <- outNAIVET_H1$p.value
  notrejectsNAIVET <- ifelse(pvaluesNAIVET_H1>0.05,1,0)
  
  
  # F-specification 1
  
  # F-test
  
  # H0
  M1_analysis_H0 <- analysis1(dataset = data_strPBR_H0)
  M1_data_strPBR_dummy_H0 <- M1_analysis_H0$db
  M1_formula_str_H0 <- M1_analysis_H0$formula
  M1_hypothesis_H0 <- M1_analysis_H0$hypoth
  
  M1_modelWALD_H0 <- lm(as.formula(paste("y_H0", M1_formula_str_H0)), data = M1_data_strPBR_dummy_H0)
  if(any(is.na(coef(M1_modelWALD_H0)))) {
    M1_pvaluesWALD_H0 <- NA
    M1_rejectsWALD <- NA
    M1_pvaluesPERM_H0 <- NA
    M1_rejectsPERM <- NA
    M1_nulls_H0 <- NA
    M1_attempts_H0 <- NA
  } 
  else {
    M1_outWALD_H0 <- linearHypothesis(M1_modelWALD_H0, M1_hypothesis_H0)
    M1_S_obsWALD_H0 <- M1_outWALD_H0$F[2] # observed statistics under H0
    M1_pvaluesWALD_H0 <- M1_outWALD_H0$`Pr(>F)`[2]
    M1_rejectsWALD <- ifelse(M1_pvaluesWALD_H0<=0.05,1,0)
    
    # Permutation test (H0)
    M1_outPERM_H0 <- perm_test_pval_IA(dataset1 = data_strPBR1, dataset2 = data_strPBR2_H0, 
                                       crit_val = z_c, drop_trt = drop_trt_H0, stat_obs = M1_S_obsWALD_H0, 
                                       hyp = 0, L = perms, max_attempts = max_atts, analysis = 0)
    M1_pvaluesPERM_H0 <- M1_outPERM_H0$pvalue
    M1_rejectsPERM <- ifelse(M1_pvaluesPERM_H0<=0.05,1,0)
    M1_nulls_H0 <- M1_outPERM_H0$nulls
    M1_attempts_H0 <- M1_outPERM_H0$total_attempts
  }
  
  # H1
  M1_analysis_H1 <- analysis1(dataset = data_strPBR_H1)
  M1_data_strPBR_dummy_H1 <- M1_analysis_H1$db
  M1_formula_str_H1 <- M1_analysis_H1$formula
  M1_hypothesis_H1 <- M1_analysis_H1$hypoth
  
  M1_modelWALD_H1 <- lm(as.formula(paste("y_H1", M1_formula_str_H1)), data = M1_data_strPBR_dummy_H1)
  if(any(is.na(coef(M1_modelWALD_H1)))) {
    M1_pvaluesWALD_H1 <- NA
    M1_notrejectsWALD <- NA
    M1_pvaluesPERM_H1 <- NA
    M1_notrejectsPERM <- NA
    M1_nulls_H1 <- NA
    M1_attempts_H1 <- NA
  } 
  else {
    M1_outWALD_H1 <- linearHypothesis(M1_modelWALD_H1, M1_hypothesis_H1)
    M1_S_obsWALD_H1 <- M1_outWALD_H1$F[2] # observed statistics under H1
    M1_pvaluesWALD_H1 <- M1_outWALD_H1$`Pr(>F)`[2]
    M1_notrejectsWALD <- ifelse(M1_pvaluesWALD_H1>0.05,1,0)
    
    # Permutation test (H1)
    M1_outPERM_H1 <- perm_test_pval_IA(dataset1 = data_strPBR1, dataset2 = data_strPBR2_H1, 
                                       crit_val = z_c, drop_trt = drop_trt_H1, stat_obs = M1_S_obsWALD_H1, 
                                       hyp = 1, L = perms, max_attempts = max_atts, analysis = 0)
    M1_pvaluesPERM_H1 <- M1_outPERM_H1$pvalue
    M1_notrejectsPERM <- ifelse(M1_pvaluesPERM_H1>0.05,1,0)
    M1_nulls_H1 <- M1_outPERM_H1$nulls
    M1_attempts_H1 <- M1_outPERM_H1$total_attempts
  }
  
  
  # F-specification 2
  
  # F-test
  
  # H0
  M2_analysis_H0 <- analysis2(dataset = data_strPBR_H0)
  M2_data_strPBR_dummy_H0 <- M2_analysis_H0$db
  M2_formula_str_H0 <- M2_analysis_H0$formula
  M2_hypothesis_H0 <- M2_analysis_H0$hypoth
  
  M2_modelWALD_H0 <- lm(as.formula(paste("y_H0", M2_formula_str_H0)), data = M2_data_strPBR_dummy_H0)
  if(any(is.na(coef(M2_modelWALD_H0)))) {
    M2_pvaluesWALD_H0 <- NA
    M2_rejectsWALD <- NA
    M2_pvaluesPERM_H0 <- NA
    M2_rejectsPERM <- NA
    M2_nulls_H0 <- NA
    M2_attempts_H0 <- NA
  } 
  else {
    M2_outWALD_H0 <- linearHypothesis(M2_modelWALD_H0, M2_hypothesis_H0)
    M2_S_obsWALD_H0 <- M2_outWALD_H0$F[2] # observed statistics under H0
    M2_pvaluesWALD_H0 <- M2_outWALD_H0$`Pr(>F)`[2]
    M2_rejectsWALD <- ifelse(M2_pvaluesWALD_H0<=0.05,1,0)
    
    # Permutation test (H0)
    M2_outPERM_H0 <- perm_test_pval_IA(dataset1 = data_strPBR1, dataset2 = data_strPBR2_H0, 
                                       crit_val = z_c, drop_trt = drop_trt_H0, stat_obs = M2_S_obsWALD_H0, 
                                       hyp = 0, L = perms, max_attempts = max_atts, analysis = 1)
    M2_pvaluesPERM_H0 <- M2_outPERM_H0$pvalue
    M2_rejectsPERM <- ifelse(M2_pvaluesPERM_H0<=0.05,1,0)
    M2_nulls_H0 <- M2_outPERM_H0$nulls
    M2_attempts_H0 <- M2_outPERM_H0$total_attempts
  }
  
  # H1
  M2_analysis_H1 <- analysis2(dataset = data_strPBR_H1)
  M2_data_strPBR_dummy_H1 <- M2_analysis_H1$db
  M2_formula_str_H1 <- M2_analysis_H1$formula
  M2_hypothesis_H1 <- M2_analysis_H1$hypoth
  
  M2_modelWALD_H1 <- lm(as.formula(paste("y_H1", M2_formula_str_H1)), data = M2_data_strPBR_dummy_H1)
  if(any(is.na(coef(M2_modelWALD_H1)))) {
    M2_pvaluesWALD_H1 <- NA
    M2_notrejectsWALD <- NA
    M2_pvaluesPERM_H1 <- NA
    M2_notrejectsPERM <- NA
    M2_nulls_H1 <- NA
    M2_attempts_H1 <- NA
  } 
  else {
    M2_outWALD_H1 <- linearHypothesis(M2_modelWALD_H1, M2_hypothesis_H1)
    M2_S_obsWALD_H1 <- M2_outWALD_H1$F[2] # observed statistics under H1
    M2_pvaluesWALD_H1 <- M2_outWALD_H1$`Pr(>F)`[2]
    M2_notrejectsWALD <- ifelse(M2_pvaluesWALD_H1>0.05,1,0)
    
    # Permutation test (H1)
    M2_outPERM_H1 <- perm_test_pval_IA(dataset1 = data_strPBR1, dataset2 = data_strPBR2_H1, 
                                       crit_val = z_c, drop_trt = drop_trt_H1, stat_obs = M2_S_obsWALD_H1, 
                                       hyp = 1, L = perms, max_attempts = max_atts, analysis = 1)
    M2_pvaluesPERM_H1 <- M2_outPERM_H1$pvalue
    M2_notrejectsPERM <- ifelse(M2_pvaluesPERM_H1>0.05,1,0)
    M2_nulls_H1 <- M2_outPERM_H1$nulls
    M2_attempts_H1 <- M2_outPERM_H1$total_attempts
  }
  
  
  # F-specification 3
  
  # F-test
  
  # H0
  M3_analysis_H0 <- analysis3(dataset = data_strPBR_H0)
  M3_data_strPBR_dummy_H0 <- M3_analysis_H0$db
  M3_formula_str_H0 <- M3_analysis_H0$formula
  M3_hypothesis_H0 <- M3_analysis_H0$hypoth
  
  M3_modelWALD_H0 <- lm(as.formula(paste("y_H0", M3_formula_str_H0)), data = M3_data_strPBR_dummy_H0)
  if(any(is.na(coef(M3_modelWALD_H0)))) {
    M3_pvaluesWALD_H0 <- NA
    M3_rejectsWALD <- NA
    M3_pvaluesPERM_H0 <- NA
    M3_rejectsPERM <- NA
    M3_nulls_H0 <- NA
    M3_attempts_H0 <- NA
  } 
  else {
    M3_outWALD_H0 <- linearHypothesis(M3_modelWALD_H0, M3_hypothesis_H0)
    M3_S_obsWALD_H0 <- M3_outWALD_H0$F[2] # observed statistics under H0
    M3_pvaluesWALD_H0 <- M3_outWALD_H0$`Pr(>F)`[2]
    M3_rejectsWALD <- ifelse(M3_pvaluesWALD_H0<=0.05,1,0)
    
    # Permutation test (H0)
    M3_outPERM_H0 <- perm_test_pval_IA(dataset1 = data_strPBR1, dataset2 = data_strPBR2_H0, 
                                       crit_val = z_c, drop_trt = drop_trt_H0, stat_obs = M3_S_obsWALD_H0, 
                                       hyp = 0, L = perms, max_attempts = max_atts, analysis = 2)
    M3_pvaluesPERM_H0 <- M3_outPERM_H0$pvalue
    M3_rejectsPERM <- ifelse(M3_pvaluesPERM_H0<=0.05,1,0)
    M3_nulls_H0 <- M3_outPERM_H0$nulls
    M3_attempts_H0 <- M3_outPERM_H0$total_attempts
  }
  
  # H1
  M3_analysis_H1 <- analysis3(dataset = data_strPBR_H1)
  M3_data_strPBR_dummy_H1 <- M3_analysis_H1$db
  M3_formula_str_H1 <- M3_analysis_H1$formula
  M3_hypothesis_H1 <- M3_analysis_H1$hypoth
  
  M3_modelWALD_H1 <- lm(as.formula(paste("y_H1", M3_formula_str_H1)), data = M3_data_strPBR_dummy_H1)
  if(any(is.na(coef(M3_modelWALD_H1)))) {
    M3_pvaluesWALD_H1 <- NA
    M3_notrejectsWALD <- NA
    M3_pvaluesPERM_H1 <- NA
    M3_notrejectsPERM <- NA
    M3_nulls_H1 <- NA
    M3_attempts_H1 <- NA
  } 
  else {
    M3_outWALD_H1 <- linearHypothesis(M3_modelWALD_H1, M3_hypothesis_H1)
    M3_S_obsWALD_H1 <- M3_outWALD_H1$F[2] # observed statistics under H1
    M3_pvaluesWALD_H1 <- M3_outWALD_H1$`Pr(>F)`[2]
    M3_notrejectsWALD <- ifelse(M3_pvaluesWALD_H1>0.05,1,0)
    
    # Permutation test (H1)
    M3_outPERM_H1 <- perm_test_pval_IA(dataset1 = data_strPBR1, dataset2 = data_strPBR2_H1, 
                                       crit_val = z_c, drop_trt = drop_trt_H1, stat_obs = M3_S_obsWALD_H1, 
                                       hyp = 1, L = perms, max_attempts = max_atts, analysis = 2)
    M3_pvaluesPERM_H1 <- M3_outPERM_H1$pvalue
    M3_notrejectsPERM <- ifelse(M3_pvaluesPERM_H1>0.05,1,0)
    M3_nulls_H1 <- M3_outPERM_H1$nulls
    M3_attempts_H1 <- M3_outPERM_H1$total_attempts
  }
  
  
  tibble::tibble(drop1_H0=drops_H0[1],drop2_H0=drops_H0[2],drop3_H0=drops_H0[3],drop4_H0=drops_H0[4],
                 drop1_H1=drops_H1[1],drop2_H1=drops_H1[2],drop3_H1=drops_H1[3],drop4_H1=drops_H1[4],
                 pvaluesNAIVET_H0, rejectsNAIVET, pvaluesNAIVET_H1, notrejectsNAIVET,
                 M1_pvaluesWALD_H0, M1_rejectsWALD, M1_pvaluesWALD_H1, M1_notrejectsWALD,
                 M1_pvaluesPERM_H0, M1_rejectsPERM, M1_nulls_H0, M1_attempts_H0, 
                 M1_pvaluesPERM_H1, M1_notrejectsPERM, M1_nulls_H1, M1_attempts_H1,
                 M2_pvaluesWALD_H0, M2_rejectsWALD, M2_pvaluesWALD_H1, M2_notrejectsWALD,
                 M2_pvaluesPERM_H0, M2_rejectsPERM, M2_nulls_H0, M2_attempts_H0, 
                 M2_pvaluesPERM_H1, M2_notrejectsPERM, M2_nulls_H1, M2_attempts_H1,
                 M3_pvaluesWALD_H0, M3_rejectsWALD, M3_pvaluesWALD_H1, M3_notrejectsWALD,
                 M3_pvaluesPERM_H0, M3_rejectsPERM, M3_nulls_H0, M3_attempts_H0, 
                 M3_pvaluesPERM_H1, M3_notrejectsPERM, M3_nulls_H1, M3_attempts_H1
  )
  
}


## Parallelization
future::availableCores() # check the number of available cores 
plan(multisession, workers = future::availableCores() - 2)

## Running simulation studies
# Scenario 1 and exclusion of patients
Ntrials <- 10000
rslt1 <- furrr::future_map(1:Ntrials, ~{
  simulate_one1(N = 200, M = 4, prevalence = c(1/4,1/4,1/4,1/4), trt_labels = c(1,2), 
                probs = c(0.5,0.5), bs = 4, alpha = 1, sigma = c(0.5,0.5,0.5,0.5), tau_H0 = c(0,0,0,0),
                tau_H1 = c(0.40,0.40,0.40,0.40), drop_trt = 4, perms = 1000) %>%
    mutate(sim=.x)
  }, .progress = TRUE, .options = furrr_options(seed = 1234))
rslt1_df <- list_rbind(rslt1)

# Scenario 1 and orientation of patients
rslt2 <- furrr::future_map(1:Ntrials, ~{
  simulate_one2(N = 200, M = 4, prevalence = c(1/4,1/4,1/4,1/4), trt_labels = c(1,2), 
                probs = c(0.5,0.5), bs = 4, alpha = 1, sigma = c(0.5,0.5,0.5,0.5),
                tau_H0 = matrix(c(0,0,0,0,
                                 0,0,0,0,
                                 0,0,0,0,
                                 0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE), 
                tau_H1 = matrix(c(0.40,0,0,0,
                                 0,0.40,0,0,
                                 0,0,0.40,0,
                                 0.40,0.40,0.40,0.40), nrow = 4, ncol = 4, byrow = TRUE), 
                drop_trt = 4, probs_orient = c(1/3,1/3,1/3),  perms = 1000) %>%
    mutate(sim=.x)
  }, .progress = TRUE, .options = furrr_options(seed = 1234))
rslt2_df <- list_rbind(rslt2)

# Scenario 2 and exclusion of patients
rslt3 <- furrr::future_map(1:Ntrials, ~{
  simulate_one3(N = 200, M = 4, prevalence = c(1/4,1/4,1/4,1/4), trt_labels = c(1,2), 
                probs = c(0.5,0.5), bs = 4, alpha = 1,  sigma = c(0.5,0.5,0.5,0.5), 
                tau_H0 = c(0,0,0,0), tau_H1 = c(0.40,0.40,0.40,-0.40), z_c = 0,
                perms = 1000, max_atts = 100000) %>%
    mutate(sim=.x)
  }, .progress = TRUE, .options = furrr_options(seed = 1234))
rslt3_df <- list_rbind(rslt3)

# Scenario 2 and orientation of patients
rslt4 <- furrr::future_map(1:Ntrials, ~{
  simulate_one4(N = 200, M = 4, prevalence = c(1/4,1/4,1/4,1/4), trt_labels = c(1,2), 
               probs = c(0.5,0.5), bs = 4, alpha = 1, sigma = c(0.5,0.5,0.5,0.5),
               tau_H0 = matrix(c(0,0,0,0,
                                 0,0,0,0,
                                 0,0,0,0,
                                 0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE), 
               tau_H1 = matrix(c(0.40,0.40,0.40,-0.40,
                                 0.40,0.40,0.40,-0.40,
                                 0.40,0.40,0.40,-0.40,
                                 0.40,0.40,0.40,-0.40), nrow = 4, ncol = 4, byrow = TRUE), 
               z_c = 0, probs_orient = c(1/3,1/3,1/3), perms = 1000, max_atts = 100000) %>%
    mutate(sim=.x)
  }, .progress = TRUE, .options = furrr_options(seed = 1234))
rslt4_df <- list_rbind(rslt4)




