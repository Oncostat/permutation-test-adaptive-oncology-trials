## User-defined functions

## Stratified block randomization procedure: stratification within strata and stages.
## It returns the dataset randomized (and indication of blocks).
## Arguments:
## * 'dataset_torand': the dataset to be randomized, which must contain the columns 
## 'stratum_assign' and 'stage';
## * 'trt_seq': vector of labels of possible treatments;
## * 'prob_each': vector of probabilities of assigning each treatment;
## * 'blocksize': size of blocks patients are assigned into.
my_strat_block_rand <- function(dataset_torand, trt_seq = c(1,2), prob_each = c(0.5,0.5), 
                                blocksize = 4) {
  
  # controls
  if (!is.data.frame(dataset_torand))
    stop("dataset_torand must be a data.frame", call. = FALSE)
  if (nrow(dataset_torand) == 0)
    stop("dataset_torand must contain at least one subject", call. = FALSE)
  required_cols <- c("stratum_assign","stage")
  missing_cols <- setdiff(required_cols, names(dataset_torand))
  if(length(missing_cols) > 0)
    stop(paste0("dataset_torand must contain the following columns: ",
                paste(missing_cols, collapse = ", ")), call. = FALSE)
  if(length(trt_seq) < 2 || anyNA(trt_seq) || anyDuplicated(trt_seq) > 0)
    stop("trt_seq must contain at least 2 distinct, non-missing labels", call. = FALSE)
  if (length(prob_each) != length(trt_seq) || !is.numeric(prob_each) || 
      any(!is.finite(prob_each)) || any(prob_each < 0) || any(prob_each > 1))
    stop(paste0("prob_each must be a numeric vector of ", length(trt_seq), " probabilities"), 
         call. = FALSE)
  if(abs(sum(prob_each) - 1) > 1e-8)
    stop(paste0("prob_each must sum to 1 (current sum: ", sum(prob_each), ")"), call. = FALSE)
  if(length(blocksize) != 1 || !is.numeric(blocksize) || !is.finite(blocksize) ||
     blocksize != round(blocksize) || blocksize < 1)
    stop("blocksize must be a positive integer", call. = FALSE)
  
  stratum_counts <- table(dataset_torand$stratum_assign)
  if(any(blocksize > stratum_counts)){
    warning("blocksize is greater than the number of subjects in one or more strata. In those strata, block randomization will effectively act as simple randomization.", call. = FALSE)
  }
  
  
  # FOR EACH STRATUM assign patients into BLOCKS of size 'blocksize'
  stratum_assign_ID <- unique(dataset_torand$stratum_assign) 
  result_list <- list()
  
  for(i in 1:length(stratum_assign_ID)){
    stratum_assign_i <- dataset_torand[dataset_torand$stratum_assign==stratum_assign_ID[i],]
    n_stratum_assign_i <- nrow(stratum_assign_i)
    
    random_order <- sample(1:n_stratum_assign_i)
    stratum_assign_i_shuffled <- stratum_assign_i[random_order, ]
    
    n_complete_blocks <- as.integer(n_stratum_assign_i/blocksize) 
    n_remaining <- n_stratum_assign_i %% blocksize
    
    block_assign <- c()
    if (n_complete_blocks > 0) {
      block_assign <- rep(1:n_complete_blocks, each = blocksize)
    }
    if (n_remaining > 0) {
      block_assign <- c(block_assign, rep(n_complete_blocks + 1, n_remaining))
    }
    
    # and FOR EACH BLOCK assign patients to the treatments based on 'prob_each' 
    block_ID <- unique(block_assign) 
    trt_assign <- numeric(n_stratum_assign_i)
    
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
    
    stratum_assign_rslt <- stratum_assign_i_shuffled
    stratum_assign_rslt$block <- paste0("M",stratum_assign_ID[i], "_B", block_assign, 
                                        "_P", stratum_assign_i_shuffled$stage[i])
    stratum_assign_rslt$trt <- trt_assign
    
    result_list[[i]] <- stratum_assign_rslt
  }
  
  datasetstrPBR <- do.call(rbind, result_list)
  rownames(datasetstrPBR) <- NULL
  
  return(data.frame(datasetstrPBR))
}

## Orientation of patients belonging to the dropped stratum through block randomization.
## It returns the dataset oriented.
## Arguments:
## * 'dataset_toorient': the dataset to orient;
## * 'strata_seq': vector of labels of remaining strata;
## * 'probs_orient': vector of probabilities of assigning a subject to each stratum.
my_block_rand <- function(dataset_toorient, strata_seq = c(1,2,3), probs_orient) {
  
  # controls
  if (!is.data.frame(dataset_toorient))
    stop("dataset_toorient must be a data.frame", call. = FALSE)
  if (nrow(dataset_toorient) == 0)
    stop("dataset_toorient must contain at least one subject", call. = FALSE)
  if (length(strata_seq) < 1 || anyNA(strata_seq) || anyDuplicated(strata_seq) > 0)
    stop("strata_seq must contain at least 1 distinct, non-missing label", call. = FALSE)
  if (length(probs_orient) != length(strata_seq) || !is.numeric(probs_orient) || 
      any(!is.finite(probs_orient)) || any(probs_orient < 0) || any(probs_orient > 1))
    stop(paste0("probs_orient must be a numeric vector of ", length(strata_seq),
                " probabilities"), call. = FALSE)
  if(abs(sum(probs_orient) - 1) > 1e-8)
    stop(paste0("probs_orient must sum to 1 (current sum: ", sum(probs_orient), ")"), 
         call. = FALSE)
  
  
  n_tot <- nrow(dataset_toorient)
  
  random_order <- sample(1:n_tot)
  data_shuffled <- dataset_toorient[random_order, ]
  
  blocksize <- length(strata_seq) # blocks of size equal to the number of remaining strata
  
  n_complete_blocks <- as.integer(n_tot/blocksize) 
  n_remaining <- n_tot %% blocksize
  
  block_assign <- c()
  if (n_complete_blocks > 0) {
    block_assign <- rep(1:n_complete_blocks, each = blocksize)
  }
  if (n_remaining > 0) {
    block_assign <- c(block_assign, rep(n_complete_blocks + 1, n_remaining))
  }
  
  # for each block assign patients to strata based on 'probs_orient' 
  block_ID <- unique(block_assign) 
  stratum_assign <- numeric(n_tot)
  
  for (j in 1:length(block_ID)){
    block_j <- which(block_assign == j)
    n_block_j <- length(block_j)
    
    n_each_stratum <- floor(n_block_j * probs_orient) 
    n_assigned <- sum(n_each_stratum)
    n_remainder <- n_block_j - n_assigned
    
    if (n_remainder > 0) {
      probs_orient_fixed <- ((probs_orient * n_block_j) - n_each_stratum)/n_remainder
      block_stratum <- c(rep(strata_seq, n_each_stratum),
                         sample(strata_seq, n_remainder, prob = probs_orient_fixed, 
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
  data_rslt$stratum_assign <- stratum_assign
  rownames(data_rslt) <- NULL
  
  return(data.frame(data_rslt))
}

## Data generation for continuous setting (stratified block randomization procedure).
## It returns the generated dataset.
## Arguments:
## * 'type_adapt': type of adaptation ('EXT' or 'IA');
## * 'stageII_approach': approach chosen for patients belonging to the dropped stratum
## in stage II ('DROP' or 'ORIENT');
## * 'N': total patients planned to be enrolled;
## * 'N1': total patients enrolled at stage I (before adaptation);
## * 'M': number of molecular alteration considered;
## * 'prevalence': vector of prevalence of each molecular alteration;
## * 'trt_labels': vector of treatment labels. The first entry must be the label 
## of the control arm;
## * 'probs': vector of probabilities of assigning each treatment;
## * 'bs': size of blocks patients are assigned into by the stratified block randomization;
## * 'alpha': control effect;
## * 'sigma': vector of alteration-specific prognostic effects;
## * 'tau': matrix of alteration*stratum-specific treatment effects;
## * 'drop_trt': molecular alteration to be dropped in the 'EXT' scenario;
## * 'z_c': futility threshold to be used in the 'IA' scenario';
## * the other parameters are described in the previous functions.
data_generation_cont <- function(type_adapt, stageII_approach, 
                                 N, N1, M, prevalence, trt_labels, probs, bs, 
                                 alpha, sigma, tau, drop_trt, z_c, probs_orient){
  
  # STAGE I
  
  # data generation
  molec_alter1 <- sample(1:M, prob = prevalence, size = N1, replace = T)
  stage1 <- rep(1,N1)
  data1 <- data.frame(stratum_orig = molec_alter1, stratum_assign = molec_alter1, 
                      stage = stage1)
  data1$ID <- 1:N1
  
  # treatment assignment (stratified block randomization)
  data_strPBR1 <- my_strat_block_rand(dataset_torand = data1, trt_seq = trt_labels, 
                                      prob_each = probs, blocksize = bs)
  data_strPBR1 <- data_strPBR1[order(data_strPBR1$ID), ]
  
  # continuous outcome generation
  data_strPBR1$y <- rep(NA_real_,N1)
  for(j in 1:N1){
    data_strPBR1$err[j] <- rnorm(1,0,1)
    data_strPBR1$y[j] <- alpha + sigma[data_strPBR1$stratum_orig[j]] + 
      tau[data_strPBR1$stratum_orig[j],data_strPBR1$stratum_assign[j]]*(data_strPBR1$trt[j]-1) + 
      data_strPBR1$err[j]
  }
  
  
  # INTERIM ANALYSIS
  if(type_adapt=="IA"){
    ctrl_label <- trt_labels[1]
    trt_label <- trt_labels[2]
    stddiff <- sapply(1:M, function(s_assign) {
      stratum_ctrl <- data_strPBR1[data_strPBR1$stratum_assign == s_assign & data_strPBR1$trt == ctrl_label, ]
      stratum_trt <- data_strPBR1[data_strPBR1$stratum_assign == s_assign & data_strPBR1$trt == trt_label, ]
      mean_diff <- mean(stratum_trt$y) - mean(stratum_ctrl$y)
      sd_ctrl <- sd(stratum_ctrl$y) # extension required for different variances (heteroschedasticity)
      return(mean_diff / sd_ctrl)
    })
    
    min_diff <- min(stddiff)
    drops <- ifelse(stddiff == min_diff & stddiff <= z_c, 1, 0)
    
    pooled_stratum <- which(drops == 0)
    
    drop_trt <- ifelse(length(pooled_stratum) == M, 0, setdiff(1:M, pooled_stratum))
  }
  
  
  # STAGE II
  
  # data generation
  data_strPBR1_keep <- data_strPBR1[data_strPBR1$stratum_orig!=drop_trt,]
  N1_k <- nrow(data_strPBR1_keep)
  N2 <- N - N1_k  
  if(stageII_approach=="DROP"){
    molec_alter2 <- sample(setdiff(1:M, drop_trt), prob = prevalence[!(seq_along(prevalence) %in% drop_trt)], 
                           N2, replace = T)
  }
  else if(stageII_approach=="ORIENT"){
    molec_alter2 <- sample(1:M, prob = prevalence, N2, replace = T)
  }
  stage2 <- rep(2, N2)
  data2 <- data.frame(stratum_orig = molec_alter2, stage = stage2)
  data2$ID <- 1:N2
  
  if(stageII_approach=="DROP" | (type_adapt=="IA" & stageII_approach=="ORIENT" & drop_trt==0)){
    data2_toassign <- data2
    data2_toassign$stratum_assign <- data2$stratum_orig
  }
  else if((type_adapt=="EXT" & stageII_approach=="ORIENT") | 
          (type_adapt=="IA" & stageII_approach=="ORIENT" & drop_trt!=0)){ 
    data2_maintain <- data2[data2$stratum_orig!=drop_trt,]
    data2_maintain$stratum_assign <- data2_maintain$stratum_orig
    data2_toorient <- data2[data2$stratum_orig==drop_trt,]
    # orientation of patients belonging to the dropped stratum (block randomization)
    data2_oriented <- my_block_rand(dataset_toorient = data2_toorient, strata_seq = setdiff(1:M, drop_trt), 
                                    probs_orient)
    data2_toassign <- rbind(data2_maintain,data2_oriented)
  }
  data2_toassign <- data2_toassign[order(data2_toassign$ID), ]
  
  # treatment assignment (stratified block randomization)
  data_strPBR2 <- my_strat_block_rand(dataset_torand = data2_toassign, trt_seq = trt_labels, 
                                      prob_each = probs, blocksize = bs)
  data_strPBR2 <- data_strPBR2[order(data_strPBR2$ID), ]
  
  # continuous outcome generation
  data_strPBR2$y <- rep(NA_real_,N2)
  for(j in 1:N2){
    data_strPBR2$err[j] <- rnorm(1,0,1)
    data_strPBR2$y[j] <- alpha + sigma[data_strPBR2$stratum_orig[j]] + 
      tau[data_strPBR2$stratum_orig[j],data_strPBR2$stratum_assign[j]]*(data_strPBR2$trt[j]-1) + 
      data_strPBR2$err[j]
  }
  
  data_strPBR <- rbind(data_strPBR1_keep,data_strPBR2)
  
  return(list(data_strPBR1 = data_strPBR1, data_strPBR2 = data_strPBR2, data_strPBR = data_strPBR, 
              drop_trtIA = drop_trt))
}

## Data generation for survival setting (minimization procedure).
## It returns the generated dataset.
## Arguments:
## * 'p_min': treatment assignment probability in minimization procedure;
## * 'adapt_time': time of adaptation (should be consistent with N1);
## * 'Tr': last recruitment time;
## * 'Ts': total study duration;
## * 'delta': censoring rate;
## * 'lambda0': baseline hazard;
## * 'coefX': efecct of covariate X (additional minimization factor);
## * 'probX': success probability of covariate X (additional minimization factor);
## * the other parameters are described in the previous functions.
data_generation_surv <- function(type_adapt, stageII_approach, 
                                 N, N1, M, prevalence, trt_labels, probs, p_min, 
                                 adapt_time, Tr, Ts, delta, lambda0, sigma, coefX, 
                                 probX, tau, drop_trt, z_c, probs_orient){
  
  # STAGE I
  
  # data generation
  molec_alter1 <- sample(1:M, prob = prevalence, size = N1, replace = T)
  X1 <- rbinom(N1, 1, probX)
  covmat1 <- data.frame(S = molec_alter1, X = X1)
  stage1 <- rep(1,N1)
  
  # entry time = time of randomization
  entry1 <- sort(runif(N1, 0, adapt_time))
  
  data1 <- data.frame(ID = 1:N1, stratum_orig = molec_alter1, stratum_assign = molec_alter1, 
                      covX = X1, stage = stage1, t0 = entry1)
  
  # treatment assignment (minimization)
  trt1 <- numeric(N1)
  trt1[1] <- sample(trt_labels, prob=probs, size=1)
  for (j in 2:N1) {
    trt1[j] <- Minirand(covmat = covmat1, j = j, covwt = rep(1, ncol(covmat1)), 
                        ratio = probs, ntrt = length(trt_labels), 
                        trtseq = trt_labels, method = "Range", result = trt1, p = p_min)
  }
  data1$trt <- trt1
  
  # time-to-event outcome generation
  t_censamm1 <- adapt_time - entry1
  xb1 <- numeric(N1)
  for(j in 1:N1){
    idx <- data1$stratum_orig[j]
    xb1[j] <- sigma[idx] + tau[idx,idx] * (data1$trt[j]-1) + coefX * data1$covX[j]
  }
  t_event1 <- -log(runif(N1)) / (lambda0 * exp(xb1))
  t_drop1 <- rexp(N1,rate=delta)
  time1 <- pmin(t_event1, t_censamm1, t_drop1)
  event1 <- as.numeric(t_event1 <= pmin(t_drop1,t_censamm1))
  
  data1$t_eventIA <- time1
  data1$eventIA <- event1
  data1$t_real <- t_event1
  data1$t_drop <- t_drop1
  
  data1 <- data1 %>%
    arrange(t0)
  
  
  # INTERIM ANALYSIS
  if(type_adapt=="IA"){
    Zscore_IA <- numeric(M)
    for (s_assign in 1:M) {
      sub_df <- subset(data1, stratum_assign == s_assign)
      sub_df <- sub_df %>%
        arrange(t0)
      sub_df$trt_bin <- sub_df$trt - 1
      sub_df$covX_f <- factor(sub_df$covX)
      cox_IA <- coxph(Surv(t_eventIA, eventIA) ~ trt_bin + strata(covX_f), data=sub_df)
      logHR_IA <- as.numeric(coef(cox_IA)["trt_bin"])        
      SE_logHR_IA <- sqrt(as.numeric(vcov(cox_IA)["trt_bin","trt_bin"]))
      Zscore_IA[s_assign] <- logHR_IA / SE_logHR_IA
    }
    worst_stratum <- which.max(Zscore_IA)
    if (Zscore_IA[worst_stratum] > z_c) {
      drop_trt <- worst_stratum
    }
    else{
      drop_trt <- 0
    } 
  }
  
  
  # STAGE II
  
  # data generation
  data1$t_event <- pmin(data1$t_real, Ts - data1$t0, data1$t_drop) # will not be used if ID not in data1_keep
  data1$event <- as.numeric(data1$t_real <= pmin(Ts - data1$t0, data1$t_drop))  # will not be used if ID not in data1_keep
  data1_keep <- data1[data1$stratum_orig!=drop_trt,]
  data1_keep$t_eventIA <- NULL
  data1_keep$eventIA <- NULL
  data1_keep$t_real <- NULL
  data1_keep$t_drop <- NULL
  N1_k <- nrow(data1_keep)
  N2 <- N - N1_k  
  if(stageII_approach=="DROP"){
    molec_alter2 <- sample(setdiff(1:M, drop_trt), prob = prevalence[!(seq_along(prevalence) %in% drop_trt)], 
                           N2, replace = T)
  }
  else if(stageII_approach=="ORIENT"){
    molec_alter2 <- sample(1:M, prob = prevalence, N2, replace = T)
  }
  X2 <- rbinom(N2, 1, probX)
  covmat2 <- data.frame(S = molec_alter2, X = X2)
  stage2 <- rep(2, N2)
  
  # entry time = time of randomization
  entry2 <- sort(runif(N2, adapt_time, Tr))
  
  data2_pre <- data.frame(ID = (N1 + 1):(N1 + N2), stratum_orig = molec_alter2, 
                          covX = X2, stage = stage2, t0 = entry2)
  
  if(stageII_approach=="DROP" | (type_adapt=="IA" & stageII_approach=="ORIENT" & drop_trt==0)){
    data2 <- data2_pre
    data2$stratum_assign <- data2$stratum_orig
  }
  if((type_adapt=="EXT" & stageII_approach=="ORIENT") | 
     (type_adapt=="IA" & stageII_approach=="ORIENT" & drop_trt!=0)){
    data2_maintain <- data2_pre[data2_pre$stratum_orig!=drop_trt,]
    data2_maintain$stratum_assign <- data2_maintain$stratum_orig
    data2_toorient <- data2_pre[data2_pre$stratum_orig==drop_trt,]
    # orientation of patients belonging to the dropped stratum (block randomization)
    data2_oriented <- my_block_rand(dataset_toorient = data2_toorient, strata_seq = setdiff(1:M, drop_trt),
                                    probs_orient)
    data2 <- rbind(data2_maintain,data2_oriented)
  }
  data2 <- data2[order(data2$ID), ]
  
  # treatment assignment (minimization)
  covmat_tot <- as.matrix(rbind(covmat1,covmat2))
  trt_tot <- c(trt1,numeric(N2))
  for (j in (N1+1):(N1+N2)) {
    trt_tot[j] <- Minirand(covmat = covmat_tot, j = j, covwt = rep(1, ncol(covmat_tot)), 
                           ratio = probs, ntrt = length(trt_labels), 
                           trtseq = trt_labels, method = "Range", result = trt_tot, p = p_min)
  }
  trt2 <- trt_tot[(N1 + 1):(N1 + N2)]
  data2$trt <- trt2
  
  # time-to-event outcome generation
  t_censamm2 <- Ts - entry2
  xb2 <- numeric(N2)
  for(j in 1:N2){
    idx1 <- data2$stratum_orig[j]
    idx2 <- data2$stratum_assign[j]
    xb2[j] <- sigma[idx1] + tau[idx1,idx2] * (data2$trt[j]-1) + coefX * data2$covX[j]
  }
  t_event2 <- -log(runif(N2)) / (lambda0 * exp(xb2))
  t_drop2 <- rexp(N2,rate=delta)
  time2 <- pmin(t_event2, t_censamm2,t_drop2)
  event2 <- as.numeric(t_event2 <= pmin(t_drop2,t_censamm2))
  
  data2$t_event <- time2
  data2$event <- event2
  
  data2 <- data2 %>%
    arrange(t0)
  
  data_tot <- rbind(data1_keep,data2)
  data_tot <- data_tot %>%
    arrange(stage, t0)
  
  return(list(data1 = data1, data2 = data2, data_tot = data_tot, drop_trtIA = drop_trt))
}

## F-statistic specifications: 
## F1: no account for shifts in prognostic or treatment effects;
## F2: accounts for shifts in prognostic effect but not in treatment effect;
## F3: accounts for both shifts in prognostic and treatment effects.
## It returns the dataset and two strings: one for model construction and one for 
## hypothesis testing.
## Arguments: 
## * 'dataset': dataset to be analysed, which must contain the columns 'stratum_assign', 
## 'trt', stage', and 'stratum_orig';
## * 'analysis': F-specification to perform (1, 2, or 3).
Fanalysis <- function(dataset, analysis = 1){
  
  # controls
  if (!is.data.frame(dataset))
    stop("dataset must be a data.frame", call. = FALSE)
  if (nrow(dataset) == 0)
    stop("dataset must contain at least one subject", call. = FALSE)
  required_cols <- c("stratum_assign","trt","stage","stratum_orig")
  missing_cols <- setdiff(required_cols, names(dataset))
  if(length(missing_cols) > 0)
    stop(paste0("dataset must contain the following columns: ",
                paste(missing_cols, collapse = ", ")), call. = FALSE)
  if (length(analysis) != 1 || !is.numeric(analysis) || !analysis %in% c(1, 2, 3))
    stop("analysis must be 1, 2 or 3", call. = FALSE)
  
  
  if(analysis==1){
    # definition of dummy variables and corresponding weights
    stratum_assign_levels <- unique(dataset$stratum_assign)
    trt_levels <- unique(dataset$trt)
    for(s_assign in stratum_assign_levels) {
      for(t in trt_levels) {
        dummy_name <- paste0(ifelse(t == trt_levels[1], "C", "T"), s_assign)
        dataset[[dummy_name]] <- as.numeric(dataset$trt == t & dataset$stratum_assign == s_assign)
      }
    }
    weights <- sapply(stratum_assign_levels, function(s_assign) {
      sum(dataset[[paste0("C", s_assign)]]) + sum(dataset[[paste0("T", s_assign)]])
    }) / nrow(dataset)
    names(weights) <- paste0("w", stratum_assign_levels)
    dummy_vars <- paste0(rep(c("C", "T"), each = length(stratum_assign_levels)),
                         rep(stratum_assign_levels, 2))
    # linear regression hypothesis formulation
    formula_str <- paste0("~ 0 + ", paste(dummy_vars, collapse = " + "))
    hypothesis_terms <- sapply(stratum_assign_levels, function(s_assign) {
      w <- weights[paste0("w", s_assign)]
      paste0(w, "*T", s_assign, " - ", w, "*C", s_assign)
    })
  }
  else if(analysis==2){
    # definition of dummy variables and corresponding weights
    stratum_assign_levels <- unique(dataset$stratum_assign)
    trt_levels <- unique(dataset$trt)
    for(s_assign in stratum_assign_levels) {
      for(t in trt_levels) {
        dummy_name <- paste0(ifelse(t == trt_levels[1], "C", "T"), s_assign)
        dataset[[dummy_name]] <- as.numeric(dataset$trt == t & dataset$stratum_assign == s_assign)
      }
    }
    period_dummy_vars <- c()
    existing_combinations <- unique(dataset[dataset$stage == 2, c("stratum_orig", "stratum_assign")])
    for(i in 1:nrow(existing_combinations)) {
      s_orig <- existing_combinations$stratum_orig[i]
      s_assign <- existing_combinations$stratum_assign[i]
      dummy_name <- paste0("P", s_orig, s_assign, "II")
      dataset[[dummy_name]] <- as.numeric(dataset$stratum_orig == s_orig & 
                                            dataset$stratum_assign == s_assign & 
                                            dataset$stage == 2)
      period_dummy_vars <- c(period_dummy_vars, dummy_name)
    }
    weights <- sapply(stratum_assign_levels, function(s_assign) {
      sum(dataset[[paste0("C", s_assign)]]) + sum(dataset[[paste0("T", s_assign)]])
    }) / nrow(dataset)
    names(weights) <- paste0("w", stratum_assign_levels)
    trt_dummy_vars <- paste0(rep(c("C", "T"), each = length(stratum_assign_levels)),
                             rep(stratum_assign_levels, 2))
    all_dummy_vars <- c(trt_dummy_vars, period_dummy_vars)
    # linear regression hypothesis formulation
    formula_str <- paste0("~ 0 + ", paste(all_dummy_vars, collapse = " + "))
    hypothesis_terms <- sapply(stratum_assign_levels, function(s_assign) {
      w <- weights[paste0("w", s_assign)]
      paste0(w, "*T", s_assign, " - ", w, "*C", s_assign)
    })
  }
  else if(analysis==3){
    # definition of dummy variables and corresponding weights
    stratum_assign_levels <- sort(unique(dataset$stratum_assign))
    trt_levels <- sort(unique(dataset$trt))
    stratum_orig_levels <- sort(unique(dataset$stratum_orig))
    for(s_assign in stratum_assign_levels) {
      for(t in trt_levels) {
        dummy_name <- paste0(ifelse(t == trt_levels[1], "C", "T"), s_assign, "I")
        dataset[[dummy_name]] <- as.numeric(dataset$trt == t & 
                                              dataset$stratum_assign == s_assign & 
                                              dataset$stage == 1)
      }
    }
    stage2_dummy_vars <- c()
    existing_combinations <- unique(dataset[dataset$stage == 2, c("stratum_orig", "stratum_assign")])
    for(i in 1:nrow(existing_combinations)) {
      s_orig <- existing_combinations$stratum_orig[i]
      s_assign <- existing_combinations$stratum_assign[i]
      for(t in trt_levels) {
        dummy_name <- paste0(ifelse(t == trt_levels[1], "C", "T"), s_orig, s_assign, "II")
        dataset[[dummy_name]] <- as.numeric(dataset$trt == t & 
                                              dataset$stratum_orig == s_orig & 
                                              dataset$stratum_assign == s_assign & 
                                              dataset$stage == 2)
        stage2_dummy_vars <- c(stage2_dummy_vars, dummy_name)
      }
    }
    n_strata <- list()
    for(s_assign in stratum_assign_levels) {
      n_strata[[paste0("n", s_assign, "I")]] <- sum(dataset[[paste0("C", s_assign, "I")]]) + 
        sum(dataset[[paste0("T", s_assign, "I")]])
      
      n_II <- 0
      for(s_orig in stratum_orig_levels) {
        c_var <- paste0("C", s_orig, s_assign, "II")
        t_var <- paste0("T", s_orig, s_assign, "II")
        if(c_var %in% names(dataset)) n_II <- n_II + sum(dataset[[c_var]])
        if(t_var %in% names(dataset)) n_II <- n_II + sum(dataset[[t_var]])
      }
      n_strata[[paste0("n", s_assign, "II")]] <- n_II
    }
    weights_I <- list()
    for(s_assign in stratum_assign_levels) {
      n_sI <- n_strata[[paste0("n", s_assign, "I")]]
      if(n_sI > 0) {
        w <- ((sum(dataset[[paste0("C", s_assign, "I")]]) + 
                 sum(dataset[[paste0("T", s_assign, "I")]])) / n_sI) / 2
        weights_I[[paste0("w", s_assign, "I")]] <- w
      }
    }
    weights_II <- list()
    for(s_assign in stratum_assign_levels) {
      n_sII <- n_strata[[paste0("n", s_assign, "II")]]
      for(s_orig in stratum_orig_levels) {
        c_var <- paste0("C", s_orig, s_assign, "II")
        t_var <- paste0("T", s_orig, s_assign, "II")
        if(c_var %in% names(dataset) && t_var %in% names(dataset)) {
          if(n_sII > 0) {
            w <- ((sum(dataset[[c_var]]) + sum(dataset[[t_var]])) / n_sII) / 2
            weights_II[[paste0("w", s_orig, s_assign, "II")]] <- w
          }
        }
      }
    }
    stage1_dummy_vars <- paste0(rep(c("C", "T"), each = length(stratum_assign_levels)), 
                                rep(stratum_assign_levels, 2), "I")
    all_dummy_vars <- c(stage1_dummy_vars, stage2_dummy_vars)
    # linear regression hypothesis formulation
    formula_str <- paste0("~ 0 + ", paste(all_dummy_vars, collapse = " + "))
    hypothesis_terms <- c()
    for(s_assign in stratum_assign_levels) {
      w_name <- paste0("w", s_assign, "I")
      if(w_name %in% names(weights_I)) {
        w <- weights_I[[w_name]]
        term <- paste0(w, "*T", s_assign, "I - ", w, "*C", s_assign, "I")
        hypothesis_terms <- c(hypothesis_terms, term)
      }
    }
    for(s_assign in stratum_assign_levels) {
      for(s_orig in stratum_orig_levels) {
        w_name <- paste0("w", s_orig, s_assign, "II")
        if(w_name %in% names(weights_II)) {
          w <- weights_II[[w_name]]
          term <- paste0(w, "*T", s_orig, s_assign, "II - ", w, "*C", s_orig, s_assign, "II")
          hypothesis_terms <- c(hypothesis_terms, term)
        }
      }
    }
  }
  
  hypothesis <- paste0(paste(hypothesis_terms, collapse = " + "), " = 0")
  
  return(list(db = dataset, formula = formula_str, hypoth = hypothesis))
  
}

## Permutation of treatment labels within blocks.
## It returns the permuted dataset.
## Arguments:
## * 'dataset': dataset to be permuted, which must contain the columns 'trt' and 'block';
## * 'seed': can be set for reproducibility (not recommended in the simulation study).
permute_trts_by_block <- function(dataset, seed = NULL) {
  
  # controls
  if (!is.data.frame(dataset))
    stop("dataset must be a data.frame", call. = FALSE)
  if (nrow(dataset) == 0)
    stop("dataset must contain at least one subject", call. = FALSE)
  required_cols <- c("trt","block")
  missing_cols <- setdiff(required_cols, names(dataset))
  if(length(missing_cols) > 0)
    stop(paste0("dataset must contain the following columns: ",
                paste(missing_cols, collapse = ", ")), call. = FALSE)
  if (!is.null(seed) &&
      (length(seed) != 1 || !is.numeric(seed) || !is.finite(seed) || seed != round(seed)))
    stop("seed must be either NULL or a finite integer", call. = FALSE)
  
  
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

## Permutation p-value estimation for external event ('EXT') scenario.
## It returns the p-values and the number of convergence issues for F1, F2, and F3.
## Arguments:
## * 'dataset': dataset to be analysed (excluding the dropped stratum), which must 
## contain the columns 'ID' and 'y';
## * 'stat_F1_obs': observed F1 test statistic;
## * 'stat_F2_obs': observed F2 test statistic;
## * 'stat_F3_obs': observed F3 test statistic;
## * 'L': number of permutations (or randomizations) to perform.
perm_test_pval_EXT <- function(dataset, stat_F1_obs = NA, stat_F2_obs = NA, stat_F3_obs = NA, L){
  
  # controls
  if (!is.data.frame(dataset))
    stop("dataset must be a data.frame", call. = FALSE)
  if (nrow(dataset) == 0)
    stop("dataset must contain at least one subject", call. = FALSE)
  required_cols <- c("ID","y")
  missing_cols <- setdiff(required_cols, names(dataset))
  if(length(missing_cols) > 0)
    stop(paste0("dataset must contain the following columns: ",
                paste(missing_cols, collapse = ", ")), call. = FALSE)
  is_na_or_finite <- function(x)
    length(x) == 1 && (is.na(x) || (is.numeric(x) && is.finite(x)))
  for (nm in c("stat_F1_obs", "stat_F2_obs", "stat_F3_obs")) {
    if (!is_na_or_finite(get(nm)))
      stop(sprintf("'%s' must be either NA or a finite number", nm), call. = FALSE)
  }
  
  
  wald_F1_L <- rep(NA_real_,L)
  I_F1_L <- rep(NA_real_,L)
  wald_F2_L <- rep(NA_real_,L)
  I_F2_L <- rep(NA_real_,L)
  wald_F3_L <- rep(NA_real_,L)
  I_F3_L <- rep(NA_real_,L)
  dataset_copy <- dataset
  for (l in 1:L){
    # treatment labels permutation
    new_dataset <- permute_trts_by_block(dataset = dataset_copy) 
    new_dataset <- new_dataset[order(new_dataset$ID), ]
    
    # treatment effects estimation
    if(!is.na(stat_F1_obs)){
      new_analysis_F1 <- Fanalysis(dataset = new_dataset, analysis = 1)
      new_dataset_dummy_F1 <- new_analysis_F1$db
      new_formula_str_F1 <- new_analysis_F1$formula
      model_F1_l <- lm(as.formula(paste("y", new_formula_str_F1)), 
                       data = new_dataset_dummy_F1)
      if(!any(is.na(coef(model_F1_l)))) {
        new_hypothesis_F1 <- new_analysis_F1$hypoth
        # computation of the test statistic at this iteration
        wald_F1_L[l] <- linearHypothesis(model_F1_l, new_hypothesis_F1)$F[2] 
        # checking the condition |Sl|>=|Sobs| at this iteration
        I_F1_L[l] <- ifelse(abs(wald_F1_L[l])>=abs(stat_F1_obs),1,0) 
      }
    }
    if(!is.na(stat_F2_obs)){
      new_analysis_F2 <- Fanalysis(dataset = new_dataset, analysis = 2)
      new_dataset_dummy_F2 <- new_analysis_F2$db
      new_formula_str_F2 <- new_analysis_F2$formula
      model_F2_l <- lm(as.formula(paste("y", new_formula_str_F2)), 
                       data = new_dataset_dummy_F2)
      if(!any(is.na(coef(model_F2_l)))) {
        new_hypothesis_F2 <- new_analysis_F2$hypoth
        # computation of the test statistic at this iteration
        wald_F2_L[l] <- linearHypothesis(model_F2_l, new_hypothesis_F2)$F[2] 
        # checking the condition |Sl|>=|Sobs| at this iteration
        I_F2_L[l] <- ifelse(abs(wald_F2_L[l])>=abs(stat_F2_obs),1,0) 
      }
    }
    if(!is.na(stat_F3_obs)){
      new_analysis_F3 <- Fanalysis(dataset = new_dataset, analysis = 3)
      new_dataset_dummy_F3 <- new_analysis_F3$db
      new_formula_str_F3 <- new_analysis_F3$formula
      model_F3_l <- lm(as.formula(paste("y", new_formula_str_F3)), 
                       data = new_dataset_dummy_F3)
      if(!any(is.na(coef(model_F3_l)))) {
        new_hypothesis_F3 <- new_analysis_F3$hypoth
        # computation of the test statistic at this iteration
        wald_F3_L[l] <- linearHypothesis(model_F3_l, new_hypothesis_F3)$F[2] 
        # checking the condition |Sl|>=|Sobs| at this iteration
        I_F3_L[l] <- ifelse(abs(wald_F3_L[l])>=abs(stat_F3_obs),1,0) 
      }
    }
  }
  if(!is.na(stat_F1_obs)){
    # computation two-sided p-value
    pvalue_F1 <- (1/sum(!is.na(I_F1_L)))*sum(I_F1_L, na.rm = T) 
    # numbers of permutations with model convergence issues
    nulls_F1 <- sum(is.na(I_F1_L))
  }
  else{
    pvalue_F1 <- NA
    nulls_F1 <- NA
  }
  if(!is.na(stat_F2_obs)){
    # computation two-sided p-value
    pvalue_F2 <- (1/sum(!is.na(I_F2_L)))*sum(I_F2_L, na.rm = T)
    # numbers of permutations with model convergence issues
    nulls_F2 <- sum(is.na(I_F2_L))
  }
  else{
    pvalue_F2 <- NA
    nulls_F2 <- NA
  }
  if(!is.na(stat_F3_obs)){
    # computation two-sided p-value
    pvalue_F3 <- (1/sum(!is.na(I_F3_L)))*sum(I_F3_L, na.rm = T) 
    # numbers of permutations with model convergence issues
    nulls_F3 <- sum(is.na(I_F3_L))
  }
  else{
    pvalue_F3 <- NA
    nulls_F3 <- NA
  }
  
  return(list(F1_pvalue = pvalue_F1, F1_nulls = nulls_F1, 
              F2_pvalue = pvalue_F2, F2_nulls = nulls_F2, 
              F3_pvalue = pvalue_F3, F3_nulls = nulls_F3))
}

## Permutation p-value estimation for interim analysis ('IA') scenario.
## It returns the number of attemps, the p-values and the number of convergence 
## issues for F1, F2, and F3.
## Arguments:
## * 'dataset1' and 'dataset2': stage-specific datasets (excluding the dropped stratum), 
## which must contain the columns 'ID', 'trt', 'stratum_assign', 'y', and 'stratum_orig';
## * 'drop_trtIA': treatment dropped by the interim analysis on the observed data;
## * 'max_attemps': maximum number of permutations (or randomizations) to reach  
## 'L' valid permutations (or randomizations);
## * the other parameters are described in the previous functions.
perm_test_pval_IA <- function(dataset1, dataset2, z_c, drop_trtIA, stat_F1_obs = NA, 
                              stat_F2_obs = NA, stat_F3_obs = NA, L, max_atts){
  
  # controls
  if (!is.data.frame(dataset1))
    stop("dataset1 must be a data.frame", call. = FALSE)
  if (nrow(dataset1) == 0)
    stop("dataset1 must contain at least one subject", call. = FALSE)
  if (!is.data.frame(dataset2))
    stop("dataset2 must be a data.frame", call. = FALSE)
  if (nrow(dataset2) == 0)
    stop("dataset2 must contain at least one subject", call. = FALSE)
  required_cols <- c("ID","trt","stratum_assign","y","stratum_orig")
  missing_cols1 <- setdiff(required_cols, names(dataset1))
  if(length(missing_cols1) > 0)
    stop(paste0("dataset1 must contain the following columns: ",
                paste(missing_cols1, collapse = ", ")), call. = FALSE)
  missing_cols2 <- setdiff(required_cols, names(dataset2))
  if(length(missing_cols2) > 0)
    stop(paste0("dataset2 must contain the following columns: ",
                paste(missing_cols2, collapse = ", ")), call. = FALSE)
  if (length(drop_trtIA) != 1 || !is.numeric(drop_trtIA) || !is.finite(drop_trtIA) ||
      drop_trtIA != round(drop_trtIA) || drop_trtIA < 0)
    stop("drop_trtIA must be a non-negative integer", call. = FALSE)
  is_na_or_finite <- function(x)
    length(x) == 1 && (is.na(x) || (is.numeric(x) && is.finite(x)))
  for (nm in c("stat_F1_obs", "stat_F2_obs", "stat_F3_obs")) {
    if (!is_na_or_finite(get(nm)))
      stop(sprintf("'%s' must be either NA or a finite number", nm), call. = FALSE)
  }
  
  
  wald_F1_L <- rep(NA_real_,L)
  I_F1_L <- rep(NA_real_,L)
  wald_F2_L <- rep(NA_real_,L)
  I_F2_L <- rep(NA_real_,L)
  wald_F3_L <- rep(NA_real_,L)
  I_F3_L <- rep(NA_real_,L)
  dataset_copy <- dataset1
  
  l <- 1  # valid permutations counter
  attempts <- 0 # number of attempts counter
  
  while(l <= L){
    attempts <- attempts + 1
    
    if(attempts > max_atts){
      warning(paste("Only ", l-1, " permutations obtained out of ", L, " permutations."))
      wald_F1_L <- wald_F1_L[seq_len(l-1)]
      I_F1_L <- I_F1_L[seq_len(l-1)]
      wald_F2_L <- wald_F2_L[seq_len(l-1)]
      I_F2_L <- I_F2_L[seq_len(l-1)]
      wald_F3_L <- wald_F3_L[seq_len(l-1)]
      I_F3_L <- I_F3_L[seq_len(l-1)]
      break
    }
    
    # treatment labels permutation
    new_dataset1 <- permute_trts_by_block(dataset = dataset_copy) 
    new_dataset1 <- new_dataset1[order(new_dataset1$ID), ]
    
    new_trt_labels <- sort(unique(new_dataset1$trt))  
    new_ctrl_label <- new_trt_labels[1]
    new_trt_label <- new_trt_labels[2]
    new_M <- length(unique(new_dataset1$stratum_assign))
    new_stddiff <- sapply(1:new_M, function(s_assign) {
      new_stratum_ctrl <- new_dataset1[new_dataset1$stratum_assign == s_assign & 
                                         new_dataset1$trt == new_ctrl_label, ]
      new_stratum_trt <- new_dataset1[new_dataset1$stratum_assign == s_assign & 
                                        new_dataset1$trt == new_trt_label, ]
      new_mean_diff <- mean(new_stratum_trt$y) - mean(new_stratum_ctrl$y)
      new_sd_pooled <- sqrt(((nrow(new_stratum_trt)-1)*var(new_stratum_trt$y) +
                               (nrow(new_stratum_ctrl)-1)*var(new_stratum_ctrl$y))/
                              (nrow(new_stratum_trt)+nrow(new_stratum_ctrl)-2))
      new_denom <- new_sd_pooled*(1/nrow(new_stratum_trt)+1/nrow(new_stratum_ctrl))
      return(new_mean_diff / new_denom)
    })
    
    new_min_diff <- min(new_stddiff)
    new_drops <- ifelse(new_stddiff == new_min_diff & new_stddiff <= z_c, 1, 0)
    
    new_pooled_stratum <- which(new_drops == 0)
    
    new_drop_trt <- ifelse(length(new_pooled_stratum) == new_M, 0, setdiff(1:new_M, new_pooled_stratum))
    
    if(identical(drop_trtIA,new_drop_trt)){ # conditioning on the interim decision
      new_dataset1 <- new_dataset1[new_dataset1$stratum_orig!=new_drop_trt,]
      
      new_dataset2 <- permute_trts_by_block(dataset = dataset2)
      new_dataset2 <- new_dataset2[order(new_dataset2$ID), ]
      
      new_dataset <- rbind(new_dataset1,new_dataset2)
      
      # treatment effects estimation
      if(!is.na(stat_F1_obs)){
        new_analysis_F1 <- Fanalysis(dataset = new_dataset, analysis = 1)
        new_dataset_dummy_F1 <- new_analysis_F1$db
        new_formula_str_F1 <- new_analysis_F1$formula
        model_F1_l <- lm(as.formula(paste("y", new_formula_str_F1)), 
                         data = new_dataset_dummy_F1)
        if(!any(is.na(coef(model_F1_l)))) {
          new_hypothesis_F1 <- new_analysis_F1$hypoth
          # computation of the test statistic at this iteration
          wald_F1_L[l] <- linearHypothesis(model_F1_l, new_hypothesis_F1)$F[2] 
          # checking the condition |Sl|>=|Sobs| at this iteration
          I_F1_L[l] <- ifelse(abs(wald_F1_L[l])>=abs(stat_F1_obs),1,0) 
        }
      }
      if(!is.na(stat_F2_obs)){
        new_analysis_F2 <- Fanalysis(dataset = new_dataset, analysis = 2)
        new_dataset_dummy_F2 <- new_analysis_F2$db
        new_formula_str_F2 <- new_analysis_F2$formula
        model_F2_l <- lm(as.formula(paste("y", new_formula_str_F2)), 
                         data = new_dataset_dummy_F2)
        if(!any(is.na(coef(model_F2_l)))) {
          new_hypothesis_F2 <- new_analysis_F2$hypoth
          # computation of the test statistic at this iteration
          wald_F2_L[l] <- linearHypothesis(model_F2_l, new_hypothesis_F2)$F[2] 
          # checking the condition |Sl|>=|Sobs| at this iteration
          I_F2_L[l] <- ifelse(abs(wald_F2_L[l])>=abs(stat_F2_obs),1,0) 
        }
      }
      if(!is.na(stat_F3_obs)){
        new_analysis_F3 <- Fanalysis(dataset = new_dataset, analysis = 3)
        new_dataset_dummy_F3 <- new_analysis_F3$db
        new_formula_str_F3 <- new_analysis_F3$formula
        model_F3_l <- lm(as.formula(paste("y", new_formula_str_F3)), 
                         data = new_dataset_dummy_F3)
        if(!any(is.na(coef(model_F3_l)))) {
          new_hypothesis_F3 <- new_analysis_F3$hypoth
          # computation of the test statistic at this iteration
          wald_F3_L[l] <- linearHypothesis(model_F3_l, new_hypothesis_F3)$F[2] 
          # checking the condition |Sl|>=|Sobs| at this iteration
          I_F3_L[l] <- ifelse(abs(wald_F3_L[l])>=abs(stat_F3_obs),1,0) 
        }
      }
      l <- l + 1
    }
  }
  if(!is.na(stat_F1_obs)){
    # computation two-sided p-value
    pvalue_F1 <- (1/sum(!is.na(I_F1_L)))*sum(I_F1_L, na.rm = T) 
    # numbers of permutations with model convergence issues
    nulls_F1 <- sum(is.na(I_F1_L))
  }
  else{
    pvalue_F1 <- NA
    nulls_F1 <- NA
  }
  if(!is.null(stat_F2_obs)){
    # computation two-sided p-value
    pvalue_F2 <- (1/sum(!is.na(I_F2_L)))*sum(I_F2_L, na.rm = T)
    # numbers of permutations with model convergence issues
    nulls_F2 <- sum(is.na(I_F2_L))
  }
  else{
    pvalue_F2 <- NA
    nulls_F2 <- NA
  }
  if(!is.null(stat_F3_obs)){
    # computation two-sided p-value
    pvalue_F3 <- (1/sum(!is.na(I_F3_L)))*sum(I_F3_L, na.rm = T) 
    # numbers of permutations with model convergence issues
    nulls_F3 <- sum(is.na(I_F3_L))
  }
  else{
    pvalue_F3 <- NA
    nulls_F3 <- NA
  }
  
  return(list(F1_pvalue = pvalue_F1, F1_nulls = nulls_F1, 
              F2_pvalue = pvalue_F2, F2_nulls = nulls_F2, 
              F3_pvalue = pvalue_F3, F3_nulls = nulls_F3,
              tot_attempts = attempts))
}

## Randomization p-value estimation for external event ('EXT') scenario.
## It returns the p-values and the number of convergence issues for unadjusted and 
## adjusted score tests.
## Arguments:
## * 'dataset': full dataset to be analysed (including the dropped stratum), which 
## must contain the colums 'ID', 'stratum_assign', 'covX', t_event', and 'event';
## * 'keepID': ID of subjects kept in the observed data;
## * 'stat_unadj_obs': observed unadjusted score test;
## * 'stat_adj_obs': observed adjusted score test;
## * 'cov_list': list of variables names to consider as covariates in the minimization procedure;
## * the other parameters are described in the previous functions.
rand_test_pval_EXT <- function(dataset, keepID, stat_unadj_obs = NA, stat_adj_obs = NA, 
                               L, cov_list, trt_labels, probs, p_min) {
  
  # controls
  if (!is.data.frame(dataset))
    stop("dataset must be a data.frame", call. = FALSE)
  if (nrow(dataset) == 0)
    stop("dataset must contain at least one subject", call. = FALSE)
  required_cols <- c("ID","stratum_assign","covX","t_event","event")
  missing_cols <- setdiff(required_cols, names(dataset))
  if(length(missing_cols) > 0)
    stop(paste0("dataset must contain the following columns: ",
                paste(missing_cols, collapse = ", ")), call. = FALSE)
  if (!is.numeric(keepID) || length(keepID) == 0 || any(!is.finite(keepID)))
    stop("keepID must be a non-empty numeric vector with finite values", call. = FALSE)
  if(anyDuplicated(keepID) > 0)
    stop("keepID must not contain duplicated IDs", call. = FALSE)
  missing_ID <- setdiff(keepID, dataset$ID)
  if(length(missing_ID) > 0)
    stop(paste0("These keepID values are not in dataset$ID: ",
                paste(missing_ID, collapse = ", ")), call. = FALSE)
  is_na_or_finite <- function(x)
    length(x) == 1 && (is.na(x) || (is.numeric(x) && is.finite(x)))
  for (nm in c("stat_unadj_obs", "stat_adj_obs")) {
    if (!is_na_or_finite(get(nm)))
      stop(sprintf("'%s' must be either NA or a finite number", nm), call. = FALSE)
  }
  if (!is.character(cov_list) || length(cov_list) == 0 || anyNA(cov_list) || 
      anyDuplicated(cov_list) > 0)
    stop("cov_list must be a non-empty character vector without duplicates", call. = FALSE)
  missing_cov <- setdiff(cov_list, colnames(dataset))
  if(length(missing_cov) > 0)
    stop(paste0("These cov_list entries are not columns of dataset: ",
                paste(missing_cov, collapse = ", ")), call. = FALSE)
  
  
  N_full <- nrow(dataset)
  
  covmat_full <- as.matrix(dataset[, cov_list])
  
  S_unadj_L <- rep(NA_real_,L)
  S_adj_L <- rep(NA_real_,L)
  I_unadj_L <- rep(NA_real_,L)
  I_adj_L <- rep(NA_real_,L)
  
  conv_issues_rand <- 0
  na_issues_rand <- 0
  
  for (l in 1:L){
    
    res <- numeric(N_full)
    res[1] <- sample(trt_labels, prob=probs, size=1)
    for(j in 2:N_full) {
      res[j] <- Minirand(
        covmat = covmat_full,
        j = j,
        covwt = rep(1,ncol(covmat_full)),
        ratio = probs,
        ntrt = length(trt_labels),
        trtseq = trt_labels,
        method = "Range",
        result = res,
        p = p_min
      )
    }
    
    data_temp <- dataset
    data_temp$trt_rand <- res
    data_analysis <- data_temp %>% filter(ID %in% keepID)
    data_analysis$trt_rand_bin <- data_analysis$trt_rand - 1
    data_analysis$stratum_assign_f <- factor(data_analysis$stratum_assign)
    data_analysis$stratum_assign_f <- droplevels(data_analysis$stratum_assign_f)
    data_analysis$covX_f <- factor(data_analysis$covX)
    data_analysis <- data_analysis %>%
      arrange(stage, t0)
    
    # unadjusted score test (or log-rank test)
    # computation of the test statistic at this iteration
    S_unadj_L[l] <- survdiff(Surv(t_event,event) ~ trt_rand_bin, data = data_analysis)$chisq
    # checking the condition |Sl|>=|Sobs| at this iteration
    I_unadj_L[l] <- ifelse(abs(S_unadj_L[l])>=abs(stat_unadj_obs),1,0) 
    
    # adjusted score test
    if(!is.na(stat_adj_obs)){
      new_strata_levels <- levels(droplevels(data_analysis$stratum_assign_f))
      new_K <- length(new_strata_levels) 
      new_ref_stratum_assign <- new_strata_levels[new_K]
      new_weights <- prop.table(table(data_analysis$stratum_assign_f)) # or prop.table(table(data_analysis$stratum_assign_f[data_analysis$event == 1]))
      new_w_ref <- as.numeric(new_weights[new_ref_stratum_assign])
      # new_int_names <- c()
      new_z_names <- c()
      new_int_ref <- as.numeric(data_analysis$stratum_assign_f == new_ref_stratum_assign) * data_analysis$trt_rand_bin
      data_analysis[[paste0("int", new_K)]] <- new_int_ref
      # new_int_names <- c(new_int_names, paste0("int", new_K))
      new_int_ref_name <- paste0("int", new_K)
      for (new_k in 1:(new_K - 1)) {
        new_s_k <- new_strata_levels[new_k]
        new_w_k <- as.numeric(new_weights[new_s_k])
        new_int_k_name <- paste0("int", new_k)
        new_int_k <- as.numeric(data_analysis$stratum_assign_f == new_s_k) * data_analysis$trt_rand_bin
        data_analysis[[new_int_k_name]] <- new_int_k
        # new_int_names <- c(new_int_names, new_int_k_name)
        new_z_k_name <- paste0("Z", new_k)
        data_analysis[[new_z_k_name]] <- new_int_k - (new_w_k / new_w_ref) * new_int_ref
        new_z_names <- c(new_z_names, new_z_k_name)
      }
      # new_int_names <- sort(new_int_names)
      new_formule_red <- as.formula(
        paste("Surv(t_event, event) ~ stratum_assign_f +", 
              paste(new_z_names, collapse = " + "), "+ covX_f")
      )
      new_formule_full <- as.formula(
        paste("Surv(t_event, event) ~ stratum_assign_f +", 
              paste(new_z_names, collapse = " + "), "+ covX_f +", new_int_ref_name)
      )
      
      conv_perm <- FALSE
      na_perm   <- FALSE
      
      new_model_red <- withCallingHandlers(
        tryCatch({
          coxph(new_formule_red, data = data_analysis)
        }, error = function(e) NULL),
        warning = function(w) {
          if (grepl("Loglik converged before variable", w$message) || 
              grepl("coefficient may be infinite", w$message)) {
            has_warning_rand <<- TRUE
            invokeRestart("muffleWarning")
          }
        }
      )
      new_model_full <- NULL
      if (is.null(new_model_red)) {
        conv_perm <- TRUE                          
      } else if (anyNA(coef(new_model_red))) {
        na_perm <- TRUE                          
      } else if (!conv_perm) {
        new_model_full <- tryCatch(
          suppressWarnings(
            coxph(new_formule_full,
                  data = data_analysis,
                  init = c(coef(new_model_red), 0),
                  control = coxph.control(iter.max = 0))
          ),
          error = function(e) NULL
        )
        if (is.null(new_model_full)) {
          conv_perm <- TRUE                        
        } else if (anyNA(coef(new_model_full)) || is.na(new_model_full$score)) {
          na_perm <- TRUE                          
        }
      }
      
      if (conv_perm) conv_issues_rand <- conv_issues_rand + 1
      if (na_perm)   na_issues_rand   <- na_issues_rand + 1
      
      if (!conv_perm && !na_perm) {
        # computation of the test statistic at this iteration
        S_adj_L[l] <- unname(new_model_full$score)
        # checking the condition |Sl|>=|Sobs| at this iteration
        I_adj_L[l] <- ifelse(abs(S_adj_L[l]) >= abs(stat_adj_obs), 1, 0)
      }
    }
  }
  
  # computation two-sided p-value
  pvalue_unadj <- (1/sum(!is.na(I_unadj_L)))*sum(I_unadj_L, na.rm = T)
  
  valid_permutations <- sum(!is.na(I_adj_L))
  if (valid_permutations > 0) {
    pvalue_adj <- sum(I_adj_L, na.rm = TRUE) / valid_permutations
  } else {
    pvalue_adj <- NA_real_
  }
  
  tibble::tibble(pvalue_unadj,pvalue_adj,conv_issues_rand,na_issues_rand)
}

## Randomized p-value estimation for interim analysis.
## Arguments:
## * 'dataset1' and 'dataset2': full stage-specific datasets (excluding the dropped stratum). 
## dataset 1 must contain the columns 'ID', 'covX', 't_eventIA', 'eventIA', 'stratum_assign', 
## 't_event', and 'event'.
## dataset 2 must contain the columns 'ID', 'covX', 'stratum_assign', 't_event', 
## and 'event';
## * the other parameters are described in the previous functions.
rand_test_pval_IA <- function(M, dataset1, dataset2, keepID, z_c, drop_trtIA, stat_unadj_obs = NA, 
                              stat_adj_obs = NA, L, max_atts, cov_list, trt_labels, probs, p_min){
  
  # controls
  if (!is.data.frame(dataset1))
    stop("dataset1 must be a data.frame", call. = FALSE)
  if (nrow(dataset1) == 0)
    stop("dataset1 must contain at least one subject", call. = FALSE)
  if (!is.data.frame(dataset2))
    stop("dataset2 must be a data.frame", call. = FALSE)
  if (nrow(dataset2) == 0)
    stop("dataset2 must contain at least one subject", call. = FALSE)
  required_cols1 <- c("ID","covX","t_eventIA","eventIA","stratum_assign","t_event","event")
  missing_cols1 <- setdiff(required_cols1, names(dataset1))
  if(length(missing_cols1) > 0)
    stop(paste0("dataset1 must contain the following columns: ",
                paste(missing_cols1, collapse = ", ")), call. = FALSE)
  required_cols2 <- c("ID","covX","stratum_assign","t_event","event")
  missing_cols2 <- setdiff(required_cols2, names(dataset2))
  if(length(missing_cols2) > 0)
    stop(paste0("dataset2 must contain the following columns: ",
                paste(missing_cols2, collapse = ", ")), call. = FALSE)
  if (!is.numeric(keepID) || length(keepID) == 0 || any(!is.finite(keepID)))
    stop("keepID must be a non-empty numeric vector with finite values", call. = FALSE)
  if(anyDuplicated(keepID) > 0)
    stop("keepID must not contain duplicated IDs", call. = FALSE)
  missing_ID <- setdiff(keepID, c(dataset1$ID, dataset2$ID))
  if(length(missing_ID) > 0)
    stop(paste0("These keepID values are neither in dataset1 nor in dataset2: ",
                paste(missing_ID, collapse = ", ")), call. = FALSE)
  if (length(drop_trtIA) != 1 || !is.numeric(drop_trtIA) || !is.finite(drop_trtIA) ||
      drop_trtIA != round(drop_trtIA) || drop_trtIA < 0)
    stop("drop_trtIA must be a non-negative integer", call. = FALSE)
  is_na_or_finite <- function(x)
    length(x) == 1 && (is.na(x) || (is.numeric(x) && is.finite(x)))
  for (nm in c("stat_unadj_obs", "stat_adj_obs")) {
    if (!is_na_or_finite(get(nm)))
      stop(sprintf("'%s' must be either NA or a finite number", nm), call. = FALSE)
  }
  if (!is.character(cov_list) || length(cov_list) == 0 || anyNA(cov_list) || 
      anyDuplicated(cov_list) > 0)
    stop("cov_list must be a non-empty character vector without duplicates", call. = FALSE)
  missing_cov1 <- setdiff(cov_list, colnames(dataset1))
  if (length(missing_cov1) > 0)
    stop(paste0("These cov_list entries are not columns of dataset1: ",
                paste(missing_cov1, collapse = ", ")), call. = FALSE)
  missing_cov2 <- setdiff(cov_list, colnames(dataset2))
  if (length(missing_cov2) > 0)
    stop(paste0("These cov_list entries are not columns of dataset2: ",
                paste(missing_cov2, collapse = ", ")), call. = FALSE)
  
  
  S_unadj_L <- rep(NA_real_,L)
  I_unadj_L <- rep(NA_real_,L)
  S_adj_L <- rep(NA_real_,L)
  I_adj_L <- rep(NA_real_,L)
  
  l <- 1  # valid permutations counter
  attempts <- 0 # number of attempts counter
  
  conv_issues_rand <- 0 
  na_issues_rand <- 0
  
  while(l <= L){
    attempts <- attempts + 1
    
    if(attempts > max_atts){
      warning(paste("Only ", l-1, " permutations obtained out of ", L, " permutations."))
      S_unadj_L <- S_unadj_L[seq_len(l-1)]
      I_unadj_L <- I_unadj_L[seq_len(l-1)]
      S_adj_L <- S_adj_L[seq_len(l-1)]
      I_adj_L <- I_adj_L[seq_len(l-1)]
      break
    }
    
    N1 <- nrow(dataset1)
    
    covmat1 <- as.matrix(dataset1[, cov_list])
    
    res1 <- numeric(N1)
    res1[1] <- sample(trt_labels, prob=probs, size=1)
    for(j in 2:N1) {
      res1[j] <- Minirand(
        covmat = covmat1,
        j = j,
        covwt = rep(1,ncol(covmat1)),
        ratio = probs,
        ntrt = length(trt_labels),
        trtseq = trt_labels,
        method = "Range",
        result = res1,
        p = p_min
      )
    }
    
    new_dataset1 <- dataset1
    new_dataset1 <- new_dataset1 %>%
      arrange(t0)
    new_dataset1$trt_rand <- res1
    
    # Interim analysis
    new_Zscore_IA <- numeric(M)
    for (new_s_assign in 1:M) {
      new_sub_df <- subset(new_dataset1, stratum_assign == new_s_assign)
      new_sub_df$covX_f <- factor(new_sub_df$covX)
      new_sub_df <- new_sub_df %>%
        arrange(t0)
      new_sub_df$trt_rand_bin <- new_sub_df$trt_rand - 1
      new_cox_IA <- coxph(Surv(t_eventIA, eventIA) ~ trt_rand_bin + strata(covX_f), 
                          data=new_sub_df)
      new_logHR_IA <- as.numeric(coef(new_cox_IA)["trt_rand_bin"])
      new_SE_logHR_IA <- sqrt(as.numeric(vcov(new_cox_IA)["trt_rand_bin","trt_rand_bin"]))
      new_Zscore_IA[new_s_assign] <- new_logHR_IA / new_SE_logHR_IA
    }
    new_worst_stratum <- which.max(new_Zscore_IA)
    if (new_Zscore_IA[new_worst_stratum] > z_c) {
      new_drop_trt <- new_worst_stratum
    }
    else{
      new_drop_trt <- 0
    }
    
    
    if(identical(drop_trtIA,new_drop_trt)){ # conditioning on the interim decision
      N2 <- nrow(dataset2)
      
      covmat2 <- as.matrix(dataset2[, cov_list])
      covmat_tot <- as.matrix(rbind(covmat1,covmat2))
      
      res_tot <- c(res1,numeric(N2))
      for(j in (N1+1):(N1+N2)) {
        res_tot[j] <- Minirand(
          covmat = covmat_tot,
          j = j,
          covwt = rep(1,ncol(covmat_tot)),
          ratio = probs,
          ntrt = length(trt_labels),
          trtseq = trt_labels,
          method = "Range",
          result = res_tot,
          p = p_min
        )
      }
      res2 <- res_tot[(N1 + 1):(N1 + N2)]
      
      new_dataset2 <- dataset2
      new_dataset2 <- new_dataset2 %>%
        arrange(t0)
      new_dataset2$trt_rand <- res2
      
      new_dataset1$t_eventIA <- NULL
      new_dataset1$eventIA <- NULL
      new_dataset1$t_real <- NULL
      new_dataset1$t_drop <- NULL
      new_dataset_tot <- rbind(new_dataset1,new_dataset2)
      data_analysis <- new_dataset_tot %>% filter(ID %in% keepID)
      data_analysis$trt_rand_bin <- data_analysis$trt_rand - 1
      data_analysis$stratum_assign_f <- factor(data_analysis$stratum_assign)
      data_analysis$stratum_assign_f <- droplevels(data_analysis$stratum_assign_f)
      data_analysis$covX_f <- factor(data_analysis$covX)
      data_analysis <- data_analysis %>%
        arrange(stage, t0)
      
      # unadjusted score test (or log-rank test)
      # computation of the test statistic at this iteration
      S_unadj_L[l] <- survdiff(Surv(t_event,event) ~ trt_rand_bin, data = data_analysis)$chisq
      # checking the condition |Sl|>=|Sobs| at this iteration
      I_unadj_L[l] <- ifelse(abs(S_unadj_L[l])>=abs(stat_unadj_obs),1,0) 
      
      # adjusted score test
      if(!is.na(stat_adj_obs)){
        new_strata_levels <- levels(droplevels(data_analysis$stratum_assign_f))
        new_K <- length(new_strata_levels) 
        new_ref_stratum_assign <- new_strata_levels[new_K]
        new_weights <- prop.table(table(data_analysis$stratum_assign_f)) # or prop.table(table(data_analysis$stratum_assign_f[data_analysis$event == 1]))
        new_w_ref <- as.numeric(new_weights[new_ref_stratum_assign])
        # new_int_names <- c()
        new_z_names <- c()
        new_int_ref <- as.numeric(data_analysis$stratum_assign_f == new_ref_stratum_assign) * data_analysis$trt_rand_bin
        data_analysis[[paste0("int", new_K)]] <- new_int_ref
        # new_int_names <- c(new_int_names, paste0("int", new_K))
        new_int_ref_name <- paste0("int", new_K)
        for (new_k in 1:(new_K - 1)) {
          new_s_k <- new_strata_levels[new_k]
          new_w_k <- as.numeric(new_weights[new_s_k])
          new_int_k_name <- paste0("int", new_k)
          new_int_k <- as.numeric(data_analysis$stratum_assign_f == new_s_k) * data_analysis$trt_rand_bin
          data_analysis[[new_int_k_name]] <- new_int_k
          # new_int_names <- c(new_int_names, new_int_k_name)
          new_z_k_name <- paste0("Z", new_k)
          data_analysis[[new_z_k_name]] <- new_int_k - (new_w_k / new_w_ref) * new_int_ref
          new_z_names <- c(new_z_names, new_z_k_name)
        }
        # new_int_names <- sort(new_int_names)
        new_formule_red <- as.formula(
          paste("Surv(t_event, event) ~ stratum_assign_f +", 
                paste(new_z_names, collapse = " + "), "+ covX_f")
        )
        new_formule_full <- as.formula(
          paste("Surv(t_event, event) ~ stratum_assign_f +", 
                paste(new_z_names, collapse = " + "), "+ covX_f +", new_int_ref_name)
        )
        
        conv_perm <- FALSE
        na_perm   <- FALSE
        
        new_model_red <- withCallingHandlers(
          tryCatch({
            coxph(new_formule_red, data = data_analysis)
          }, error = function(e) NULL),
          warning = function(w) {
            if (grepl("Loglik converged before variable", w$message) || 
                grepl("coefficient may be infinite", w$message)) {
              has_warning_rand <<- TRUE
              invokeRestart("muffleWarning")
            }
          }
        )
        new_model_full <- NULL
        if (is.null(new_model_red)) {
          conv_perm <- TRUE                          
        } else if (anyNA(coef(new_model_red))) {
          na_perm <- TRUE                          
        } else if (!conv_perm) {
          new_model_full <- tryCatch(
            suppressWarnings(
              coxph(new_formule_full,
                    data = data_analysis,
                    init = c(coef(new_model_red), 0),
                    control = coxph.control(iter.max = 0))
            ),
            error = function(e) NULL
          )
          if (is.null(new_model_full)) {
            conv_perm <- TRUE                        
          } else if (anyNA(coef(new_model_full)) || is.na(new_model_full$score)) {
            na_perm <- TRUE                          
          }
        }
        
        if (conv_perm) conv_issues_rand <- conv_issues_rand + 1
        if (na_perm)   na_issues_rand   <- na_issues_rand + 1
        
        if (!conv_perm && !na_perm) {
          # computation of the test statistic at this iteration
          S_adj_L[l] <- unname(new_model_full$score)
          # checking the condition |Sl|>=|Sobs| at this iteration
          I_adj_L[l] <- ifelse(abs(S_adj_L[l]) >= abs(stat_adj_obs), 1, 0)
        }
      }
      
      l <- l+1
    }
  }
  
  # computation two-sided p-value
  pvalue_unadj <- (1/sum(!is.na(I_unadj_L)))*sum(I_unadj_L, na.rm = T)
  
  valid_permutations <- sum(!is.na(I_adj_L))
  if (valid_permutations > 0) {
    pvalue_adj <- sum(I_adj_L, na.rm = TRUE) / valid_permutations
  } else {
    pvalue_adj <- NA_real_
  }
  
  tibble::tibble(pvalue_unadj,pvalue_adj,tot_attempts = attempts,conv_issues_rand,na_issues_rand) 
  
}

## First simulation run: data generation and hypothesis testing
## Arguments:
## * 'type_outcomes': type of generated outcomes ('continuous' or 'survival');
## * 'type_rand': type of randomization procedure ('block' for stratified block randomization
## or 'minim' for minimization); 
## the other parameters are described in the previous functions.
simulate_one <- function(type_outcome = "continuous",
                         type_rand = "block",
                         type_adapt = "EXT",
                         stageII_approach = "DROP",
                         N = 200, N1 = 100, M = 4, prevalence = c(1/4,1/4,1/4,1/4), trt_labels = c(1,2), 
                         probs = c(0.5,0.5), bs = 4, p_min = NULL, adapt_time = NULL, 
                         Tr = NULL, Ts = NULL, delta = NULL, alpha = 1, lambda0 = NULL, 
                         sigma = c(0.5,0.5,0.5,0.5), coefX = NULL, probX = NULL,
                         tau = matrix(c(0,0,0,0,
                                        0,0,0,0,
                                        0,0,0,0,
                                        0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE), 
                         drop_trt = 4, 
                         z_c = NULL, probs_orient = NULL,
                         L = 1000, max_atts = NULL){
  
  # controls
  if (length(N) != 1 || !is.numeric(N) || !is.finite(N) || N != round(N) || N <= 0)
    stop("N must be a finite and positive integer", call. = FALSE)
  if (length(M) != 1 || !is.numeric(M) || !is.finite(M) || M != round(M) || M <= 0)
    stop("M must be a finite and positive integer", call. = FALSE)
  if (length(L) != 1 || !is.numeric(L) || !is.finite(L) || L != round(L) || L <= 0)
    stop("L must be a finite and positive integer", call. = FALSE)
  
  if(length(type_outcome) != 1 || !type_outcome %in% c("continuous", "survival")) # extensions required for binary outcomes
    stop("type_outcome must be 'continuous' or 'survival'", call. = FALSE) 
  if (length(type_rand) != 1 || !type_rand %in% c("block", "minim")) # extensions required to account for other randomization procedures
    stop("type_rand must be 'block' or 'minim'", call. = FALSE)
  
  if(type_outcome == "continuous" && type_rand != "block") # extensions required to account for other randomization procedures
    stop("If type_outcome is 'continuous', type_rand must be 'block'", call. = FALSE) 
  if(type_outcome == "continuous"){
    if(length(alpha) != 1 || !is.numeric(alpha) || !is.finite(alpha))
      stop("alpha must be a finite number", call. = FALSE)
    if(!is.null(adapt_time) || !is.null(Tr) || !is.null(Ts) || !is.null(delta) || 
       !is.null(lambda0))
      stop("If type_outcome = 'continuous', adapt_time, Tr, Ts, delta, and lambda0 
           must be NULL", call. = FALSE)
  }
  if(type_rand == "block"){
    if(length(bs) != 1 || !is.numeric(bs) || !is.finite(bs) || bs != round(bs) ||
       bs < 1)
      stop("bs must be a positive integer", call. = FALSE)
    expected_counts <- bs * probs
    if (!all(abs(expected_counts - round(expected_counts)) < 1e-3)) {
      stop("bs must be proportional to the treatment ratio derived by probs.", call. = FALSE)
    }
    if(!is.null(p_min) || !is.null(coefX) || !is.null(probX))
      stop("If type_rand = 'block', p_min, coefX, and probX must be NULL", call. = FALSE)
  }
  
  if(type_outcome == "survival" && type_rand != "minim") # extensions required to account for other randomization procedures
    stop("If type_outcome is 'survival', type_rand must be 'minim'", call. = FALSE) 
  if(type_outcome == "survival"){
    if(!is.null(alpha))
      stop("If type_outcome = 'survival', alpha must be NULL", call. = FALSE)
    if (length(Ts) != 1 || !is.numeric(Ts) || !is.finite(Ts) || Ts <= 0)
      stop("Ts must be a finite and positive number", call. = FALSE)
    if (length(Tr) != 1 || !is.numeric(Tr) || !is.finite(Tr) || Tr <= 0 || Tr > Ts)
      stop("Tr must be a positive number lower or equal than Ts", call. = FALSE)
    if (length(adapt_time) != 1 || !is.numeric(adapt_time) || !is.finite(adapt_time) ||
        adapt_time <= 0 || adapt_time >= Tr)
      stop("adapt_time must be a positive number strictly between 0 and Tr", call. = FALSE)
    if (length(delta) != 1 || !is.numeric(delta) || !is.finite(delta))
      stop("delta must be a finite number", call. = FALSE)
    if (length(lambda0) != 1 || !is.numeric(lambda0) || !is.finite(lambda0) || lambda0 <= 0)
      stop("lambda0 must be a finite and positive number", call. = FALSE)
  }
  if(type_rand == "minim"){
    if(!is.null(bs))
      stop("If type_rand = 'minim', bs must be NULL", call. = FALSE)
    if(length(p_min) != 1 || !is.numeric(p_min) || !is.finite(p_min) ||
       p_min <= 0.5 || p_min >= 1)
      stop("p_min must be a number (probability) strictly between 0.5 and 1", call. = FALSE)
    if(length(probX) != 1 || !is.numeric(probX) || !is.finite(probX) ||
       probX <= 0 || probX >= 1)
      stop("probX must be a number (probability) strictly between 0 and 1", call. = FALSE)
    if(is.na(coefX) || length(p_min) != 1 || !is.numeric(coefX) || !is.finite(coefX)) # extensions required for more minimization factors
      stop("coefX must be a finite number", call. = FALSE)
  }
  
  if(length(type_adapt) != 1 || !type_adapt %in% c("EXT", "IA")) # extensions required for other/more adaptations
    stop("type_adapt must be 'EXT' or 'IA'", call. = FALSE) 
  if(type_adapt == "EXT"){
    if(length(drop_trt) != 1 || !is.numeric(drop_trt) || !is.finite(drop_trt) ||
       drop_trt != round(drop_trt) || drop_trt < 1 || drop_trt > M)
      stop(paste0("If type_adapt = 'EXT', drop_trt must be an integer between 1 and M = ", M), 
           call. = FALSE)
    if(!is.null(z_c))
      stop("If type_adapt = 'EXT', z_c must be NULL", call. = FALSE)
    if(!is.null(max_atts))
      stop("If type_adapt = 'EXT', max_atts must be NULL", call. = FALSE)
  }
  if(type_adapt == "IA"){
    if(!is.null(drop_trt))
      stop("If type_adapt = 'IA', drop_trt must be NULL", call. = FALSE)
    if(length(z_c) != 1 || !is.numeric(z_c) || !is.finite(z_c))
      stop("If type_adapt = 'IA', z_c must be a number", call. = FALSE)
    if(length(max_atts) != 1 || !is.numeric(max_atts) || !is.finite(max_atts) ||
       max_atts != round(max_atts) || max_atts < L)
      stop(paste0("If type_adapt = 'IA', max_atts must be a positive integer equal or greater than L = ", L), 
           call. = FALSE)
  }
  
  if(length(stageII_approach) != 1 || !stageII_approach %in% c("DROP", "ORIENT"))
    stop("stageII_approach must be 'DROP' or 'ORIENT'", call. = FALSE) 
  if(stageII_approach == "DROP"){
    if(!is.null(probs_orient))
      stop("If stageII_approach = 'DROP', probs_orient must be NULL", call. = FALSE)
  }
  if(stageII_approach == "ORIENT"){
    if(!is.numeric(probs_orient) || length(probs_orient) != M - 1 ||
       any(!is.finite(probs_orient)) || any(probs_orient < 0) || any(probs_orient > 1))
      stop(paste0("probs_orient must be a numeric vector of M-1 = ", M - 1,
                  " probabilities"), call. = FALSE)
    if(abs(sum(probs_orient) - 1) > 1e-8)
      stop(paste0("probs_orient must sum to 1 (current sum: ", sum(probs_orient), ")"),
           call. = FALSE)
  }
  
  if(!is.numeric(prevalence) || length(prevalence) != M ||
     any(!is.finite(prevalence)) || any(prevalence < 0) || any(prevalence > 1))
    stop(paste0("prevalence must be a numeric vector of M = ", M,
                " probabilities"), call. = FALSE)
  if (abs(sum(prevalence) - 1) > 1e-8)
    stop(paste0("prevalence must sum to 1 (current sum: ", sum(prevalence), ")"),
         call. = FALSE)
  if(length(trt_labels) != 2 || anyNA(trt_labels) || anyDuplicated(trt_labels) > 0)
    stop("trt_labels must contain exactly 2 distinct, non-missing labels", call. = FALSE)
  if (!isTRUE(all(trt_labels == c(1, 2)))) { # trt_labels transformation
    warning("trt_labels has been transformed to c(1, 2): ",
            trt_labels[1], " -> 1 (control), ",
            trt_labels[2], " -> 2 (experimental)", call. = FALSE)
    trt_labels_orig <- trt_labels
    trt_labels <- c(1, 2)
  }
  if(!is.numeric(probs) || length(probs) != 2 ||
     any(!is.finite(probs)) || any(probs < 0) || any(probs > 1))
    stop("probs must be a numeric vector of 2 probabilities", call. = FALSE)
  if(abs(sum(probs) - 1) > 1e-8)
    stop(paste0("probs must sum to 1 (current sum: ", sum(probs), ")"), call. = FALSE)
  if(!is.numeric(sigma) || length(sigma) != M || any(!is.finite(sigma)))
    stop(paste0("sigma must be a numeric vector of M = ", M, " finite values"), call. = FALSE)
  if (!is.matrix(tau) || !is.numeric(tau) || nrow(tau) != M || ncol(tau) != M ||
      any(!is.finite(tau)))
    stop(paste0("tau must be a ", M, "x", M, " numeric matrix with finite entries"),
         call. = FALSE)
  
  
  if(type_outcome=="continuous"){
    
    # DATA GENERATION
    
    df <- data_generation_cont(type_adapt, stageII_approach, 
                               N, N1, M, prevalence, trt_labels, probs, bs, 
                               alpha, sigma, tau, drop_trt, z_c, probs_orient)
    data_strPBR1 <- df$data_strPBR1
    data_strPBR2 <- df$data_strPBR2
    data_strPBR <- df$data_strPBR
    drop_trtIA <- df$drop_trtIA
    
    
    # TESTING
    
    # t-test (NAIVE)
    outNAIVET <- t.test(y ~ factor(trt), data = data_strPBR)
    pvaluesNAIVET <- outNAIVET$p.value
    rejectsNAIVET <- ifelse(pvaluesNAIVET<=0.05,1,0)
    
    # F-test (F-specification 1)
    F1_analysis <- Fanalysis(dataset = data_strPBR, analysis = 1)
    F1_data_strPBR_dummy <- F1_analysis$db
    F1_formula_str <- F1_analysis$formula
    F1_hypothesis <- F1_analysis$hypoth
    
    F1_modelWALD <- lm(as.formula(paste("y", F1_formula_str)), data = F1_data_strPBR_dummy)
    if(any(is.na(coef(F1_modelWALD)))) {
      F1_S_obsWALD <- NA
      F1_pvaluesWALD <- NA
      F1_rejectsWALD <- NA
      F1_nullWALD <- 1
      F1_pvaluesPERM <- NA
      F1_rejectsPERM <- NA
      F1_nullsPERM <- NA
    } 
    else {
      F1_outWALD <- linearHypothesis(F1_modelWALD, F1_hypothesis)
      F1_S_obsWALD <- F1_outWALD$F[2] # observed statistics
      F1_pvaluesWALD <- F1_outWALD$`Pr(>F)`[2]
      F1_rejectsWALD <- ifelse(F1_pvaluesWALD<=0.05,1,0)
      F1_nullWALD <- NA
    }
    
    # F-test (F-specification 2)
    F2_analysis <- Fanalysis(dataset = data_strPBR, analysis = 2)
    F2_data_strPBR_dummy <- F2_analysis$db
    F2_formula_str <- F2_analysis$formula
    F2_hypothesis <- F2_analysis$hypoth
    
    F2_modelWALD <- lm(as.formula(paste("y", F2_formula_str)), data = F2_data_strPBR_dummy)
    if(any(is.na(coef(F2_modelWALD)))) {
      F2_S_obsWALD <- NA
      F2_pvaluesWALD <- NA
      F2_rejectsWALD <- NA
      F2_nullWALD <- 1
      F2_pvaluesPERM <- NA
      F2_rejectsPERM <- NA
      F2_nullsPERM <- NA
    } 
    else {
      F2_outWALD <- linearHypothesis(F2_modelWALD, F2_hypothesis)
      F2_S_obsWALD <- F2_outWALD$F[2] # observed statistics 
      F2_pvaluesWALD <- F2_outWALD$`Pr(>F)`[2]
      F2_rejectsWALD <- ifelse(F2_pvaluesWALD<=0.05,1,0)
      F2_nullWALD <- NA
    }
    
    # F-test (F-specification 3)
    F3_analysis <- Fanalysis(dataset = data_strPBR, analysis = 3)
    F3_data_strPBR_dummy <- F3_analysis$db
    F3_formula_str <- F3_analysis$formula
    F3_hypothesis <- F3_analysis$hypoth
    
    F3_modelWALD <- lm(as.formula(paste("y", F3_formula_str)), data = F3_data_strPBR_dummy)
    if(any(is.na(coef(F3_modelWALD)))) {
      F3_S_obsWALD <- NA
      F3_pvaluesWALD <- NA
      F3_rejectsWALD <- NA
      F3_nullWALD <- 1
      F3_pvaluesPERM <- NA
      F3_rejectsPERM <- NA
      F3_nullsPERM <- NA
    } 
    else {
      F3_outWALD <- linearHypothesis(F3_modelWALD, F3_hypothesis)
      F3_S_obsWALD <- F3_outWALD$F[2] # observed statistics
      F3_pvaluesWALD <- F3_outWALD$`Pr(>F)`[2]
      F3_rejectsWALD <- ifelse(F3_pvaluesWALD<=0.05,1,0)
      F3_nullWALD <- NA
    }
    
    if(stageII_approach=="DROP"){
      # Permutation test 
      if(type_adapt=="EXT"){
        out_perm <- perm_test_pval_EXT(dataset = data_strPBR, stat_F1_obs = F1_S_obsWALD, 
                                       stat_F2_obs = NA, stat_F3_obs = NA, L)
        Nattempts <- NA
      }
      else if(type_adapt=="IA"){
        out_perm <- perm_test_pval_IA(dataset1 = data_strPBR1, dataset2 = data_strPBR2, 
                                      z_c, drop_trtIA, stat_F1_obs = F1_S_obsWALD, 
                                      stat_F2_obs = NA, stat_F3_obs = NA, L, max_atts)
        Nattempts <- out_perm$tot_attempts
      }
      
      F1_pvaluesPERM <- out_perm$F1_pvalue
      F1_rejectsPERM <- ifelse(F1_pvaluesPERM<=0.05,1,0)
      F1_nullsPERM <- out_perm$F1_nulls
      
      tibble::tibble(pvaluesNAIVET, rejectsNAIVET,
                     F1_pvaluesWALD, F1_rejectsWALD, F1_nullWALD, F1_pvaluesPERM, 
                     F1_rejectsPERM, F1_nullsPERM, Nattempts, drop_trtIA
      )
    }
    else if(stageII_approach=="ORIENT"){
      # Permutation test 
      if(type_adapt=="EXT"){
        out_perm <- perm_test_pval_EXT(dataset = data_strPBR, stat_F1_obs = F1_S_obsWALD,
                                       stat_F2_obs = F2_S_obsWALD, stat_F3_obs = F3_S_obsWALD, L)
        Nattempts <- NA
      }
      else if(type_adapt=="IA"){
        out_perm <- perm_test_pval_IA(dataset1 = data_strPBR1, dataset2 = data_strPBR2, 
                                      z_c, drop_trtIA, stat_F1_obs = F1_S_obsWALD, 
                                      stat_F2_obs = F2_S_obsWALD, stat_F3_obs = F3_S_obsWALD, 
                                      L, max_atts)
        Nattempts <- out_perm$tot_attempts
      }
      
      F1_pvaluesPERM <- out_perm$F1_pvalue
      F1_rejectsPERM <- ifelse(F1_pvaluesPERM<=0.05,1,0)
      F1_nullsPERM <- out_perm$F1_nulls
      
      F2_pvaluesPERM <- out_perm$F2_pvalue
      F2_rejectsPERM <- ifelse(F2_pvaluesPERM<=0.05,1,0)
      F2_nullsPERM <- out_perm$F2_nulls
      
      F3_pvaluesPERM <- out_perm$F3_pvalue
      F3_rejectsPERM <- ifelse(F3_pvaluesPERM<=0.05,1,0)
      F3_nullsPERM <- out_perm$F3_nulls
      
      tibble::tibble(pvaluesNAIVET, rejectsNAIVET,
                     F1_pvaluesWALD, F1_rejectsWALD, F1_nullWALD, 
                     F1_pvaluesPERM, F1_rejectsPERM, F1_nullsPERM,
                     F2_pvaluesWALD, F2_rejectsWALD, F2_nullWALD, 
                     F2_pvaluesPERM, F2_rejectsPERM, F2_nullsPERM,
                     F3_pvaluesWALD, F3_rejectsWALD, F3_nullWALD, 
                     F3_pvaluesPERM, F3_rejectsPERM, F3_nullsPERM, Nattempts, drop_trtIA
      )
    }
  }
  else if (type_outcome=="survival"){
    
    # DATA GENERATION
    
    df <- data_generation_surv(type_adapt, stageII_approach, 
                               N, N1, M, prevalence, trt_labels, probs, p_min, adapt_time, 
                               Tr, Ts, delta, lambda0, sigma, coefX, probX, 
                               tau, drop_trt, z_c, probs_orient)
    data1 <- df$data1
    data2 <- df$data2
    data_tot <- df$data_tot
    keepID <- data_tot$ID
    drop_trtIA <- df$drop_trtIA
    
    
    # TESTING
    
    # unadjusted score test (or log-rank test)
    data_tot$trt_bin <- data_tot$trt - 1
    out_asympt <- survdiff(Surv(t_event, event) ~ trt_bin, data=data_tot)
    S_unadj_asympt <- out_asympt$chisq
    p_unadj_asympt <- out_asympt$pvalue
    rejects_unadj_asympt <- ifelse(p_unadj_asympt<=0.05,1,0)
    
    # adjusted score test 
    data_tot$stratum_assign_f <- factor(data_tot$stratum_assign)
    data_tot$covX_f <- factor(data_tot$covX)
    
    strata_levels <- levels(droplevels(data_tot$stratum_assign_f))
    K <- length(strata_levels) 
    ref_stratum_assign <- strata_levels[K]
    weights <- prop.table(table(data_tot$stratum_assign_f)) # or prop.table(table(data_tot$stratum_assign_f[data_tot$event == 1]))
    w_ref <- as.numeric(weights[ref_stratum_assign])
    # int_names <- c()
    z_names <- c()
    int_ref <- as.numeric(data_tot$stratum_assign_f == ref_stratum_assign) * data_tot$trt_bin
    data_tot[[paste0("int", K)]] <- int_ref
    # int_names <- c(int_names, paste0("int", K))
    int_ref_name <- paste0("int", K)
    for (k in 1:(K - 1)) {
      s_k <- strata_levels[k]
      w_k <- as.numeric(weights[s_k])
      int_k_name <- paste0("int", k)
      int_k <- as.numeric(data_tot$stratum_assign_f == s_k) * data_tot$trt_bin
      data_tot[[int_k_name]] <- int_k
      # int_names <- c(int_names, int_k_name)
      z_k_name <- paste0("Z", k)
      data_tot[[z_k_name]] <- int_k - (w_k / w_ref) * int_ref
      z_names <- c(z_names, z_k_name)
    }
    # int_names <- sort(int_names)
    formule_red <- as.formula(
      paste("Surv(t_event, event) ~ stratum_assign_f +", 
            paste(z_names, collapse = " + "), "+ covX_f")
    )
    formule_full <- as.formula(
      paste("Surv(t_event, event) ~ stratum_assign_f +", 
            paste(z_names, collapse = " + "), "+ covX_f +", int_ref_name)
    )
    
    conv_issue_asympt <- 0
    na_issue_asympt   <- 0
    
    model_red <- withCallingHandlers(
      tryCatch({
        coxph(formule_red, data = data_tot)
      }, error = function(e) NULL),
      warning = function(w) {
        if (grepl("Loglik converged before variable", w$message) ||
            grepl("coefficient may be infinite", w$message)) {
          conv_issue_asympt <<- 1
          invokeRestart("muffleWarning")
        }
      }
    )
    
    model_full <- NULL
    if (is.null(model_red)) {
      conv_issue_asympt <- 1                      
    } else if (anyNA(coef(model_red))) {
      na_issue_asympt <- 1                        
    } else if (conv_issue_asympt == 0) {
      model_full <- tryCatch(
        suppressWarnings(
          coxph(formule_full,
                data = data_tot,
                init = c(coef(model_red), 0),
                control = coxph.control(iter.max = 0))
        ),
        error = function(e) NULL
      )
      if (is.null(model_full)) {
        conv_issue_asympt <- 1                   
      } else if (anyNA(coef(model_full)) || is.na(model_full$score)) {
        na_issue_asympt <- 1                    
      }
    }
    
    if (conv_issue_asympt == 1 || na_issue_asympt == 1) {
      # Randomization
      if(type_adapt=="EXT"){
        data1$t_eventIA <- NULL
        data1$eventIA <- NULL
        data1$t_real <- NULL
        data1$t_drop <- NULL
        data_full <- rbind(data1,data2)
        data_full <- data_full %>%
          arrange(stage, t0)
        out_rand <- rand_test_pval_EXT(data_full, keepID, stat_unadj_obs = S_unadj_asympt, 
                                       stat_adj_obs = NA, L, cov_list = c("stratum_orig","covX"),
                                       trt_labels, probs, p_min)
        Nattempts <- NA
      }
      else if(type_adapt=="IA"){
        out_rand <- rand_test_pval_IA(M, data1, data2, keepID, z_c, drop_trtIA, 
                                      stat_unadj_obs = S_unadj_asympt, stat_adj_obs = NA,
                                      L, max_atts, cov_list = c("stratum_orig","covX"),
                                      trt_labels, probs, p_min)
        Nattempts <- out_rand$tot_attempts
      }
      # unadjusted score test (or log-rank test)
      p_unadj_rand <- out_rand$pvalue_unadj
      rejects_unadj_rand <- ifelse(p_unadj_rand<=0.05,1,0)
      
      return(tibble::tibble(
        p_unadj_asympt = p_unadj_asympt,
        rejects_unadj_asympt = rejects_unadj_asympt,
        Nattempts = Nattempts,
        p_unadj_rand = p_unadj_rand,
        rejects_unadj_rand = rejects_unadj_rand,
        p_adj_asympt = NA_real_,
        rejects_adj_asympt = NA_real_,
        conv_issue_asympt = conv_issue_asympt,
        na_issue_asympt = na_issue_asympt,
        p_adj_rand = NA_real_,
        rejects_adj_rand = NA_real_,
        conv_issues_rand = NA_real_,
        na_issues_rand = NA_real_,
        drop_trtIA = drop_trtIA
      ))
    }
    
    S_adj_asympt <- unname(model_full$score)
    p_adj_asympt <- pchisq(S_adj_asympt, df = 1, lower.tail = FALSE)
    rejects_adj_asympt <- ifelse(p_adj_asympt<=0.05,1,0)
    
    # Randomization
    if(type_adapt=="EXT"){
      data1$t_eventIA <- NULL
      data1$eventIA <- NULL
      data1$t_real <- NULL
      data1$t_drop <- NULL
      data_full <- rbind(data1,data2)
      data_full <- data_full %>%
        arrange(stage, t0)
      out_rand <- rand_test_pval_EXT(data_full, keepID, stat_unadj_obs = S_unadj_asympt, 
                                     stat_adj_obs = S_adj_asympt, L, cov_list = c("stratum_orig","covX"),
                                     trt_labels, probs, p_min)
      Nattempts <- NA
    }
    if(type_adapt=="IA"){
      out_rand <- rand_test_pval_IA(M, data1, data2, keepID, z_c = z_c, drop_trt = drop_trtIA, 
                                    stat_unadj_obs = S_unadj_asympt, stat_adj_obs = S_adj_asympt,
                                    L, max_atts, cov_list = c("stratum_orig","covX"),
                                    trt_labels, probs, p_min)
      Nattempts <- out_rand$tot_attempts
    }
    # unadjusted score test (or log-rank test)
    p_unadj_rand <- out_rand$pvalue_unadj
    rejects_unadj_rand <- ifelse(p_unadj_rand<=0.05,1,0)
    
    # adjusted score test
    p_adj_rand <- out_rand$pvalue_adj
    rejects_adj_rand <- ifelse(p_adj_rand<=0.05,1,0)
    conv_issues_rand <- out_rand$conv_issues_rand
    na_issues_rand <- out_rand$na_issues_rand
    
    
    tibble::tibble(p_unadj_asympt,rejects_unadj_asympt,Nattempts,p_unadj_rand,rejects_unadj_rand,
                   p_adj_asympt,rejects_adj_asympt,conv_issue_asympt,na_issue_asympt,
                   p_adj_rand,rejects_adj_rand,conv_issues_rand,na_issues_rand,drop_trtIA
    )
  }
}
