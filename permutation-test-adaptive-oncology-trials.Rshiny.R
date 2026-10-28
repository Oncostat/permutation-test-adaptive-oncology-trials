# install.packages("shiny")
# install.packages("bslib")
# install.packages("DT")         # editable tables
# install.packages("car")        
# install.packages("tibble")     
# install.packages("future")     
# install.packages("furrr")      
# install.packages("survival")   
# install.packages("dplyr")      
# install.packages("Minirand") 
library(shiny)
library(bslib)
library(DT)         
library(car)        
library(tibble)     
library(future)     
library(furrr)      
library(survival)   
library(dplyr)      
library(Minirand)   

# Auxiliary functions saved in aux_functions.R file.
# Preferred location: aux_functions.R in the same folder as this file.
aux_file <- "//nas-01/SBE_ETUDES/SAFIR02/Translation_Projects/Christina/PhD project/abstracts and papers/Paper axis 1/github code/aux_functions.R"
source(aux_file)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

# Workers must be strictly lower than the number of available cores
n_cores     <- as.integer(future::availableCores())
max_workers <- max(1L, n_cores - 1L)

# Maximum and minimum number of rows (molecular alterations) in the table
max_M <- 20L
min_M <- 1L

# Packages needed on the parallel workers
worker_pkgs <- function(p) {
  if (identical(p$type_outcome, "continuous")) c("car", "tibble")
  else c("car", "tibble", "survival", "dplyr", "Minirand")
}

# Wilson confidence interval for a proportion x/n
wilson_interval <- function(x, n, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2)
  p <- x / n
  
  center <- (p + z^2/(2*n)) / (1 + z^2/n)
  margin <- (z / (1 + z^2/n)) * sqrt(p*(1-p)/n + z^2/(4*n^2))
  
  c(lower = center - margin, upper = center + margin)
}

# Type I error and power of each test, in %, and Wilson interval centred on the
# nominal level (5%): it tells whether the estimated type I error is compatible with 5%.
#  rej0 = 1/0 indicator of rejection under H0 (NA = test not computable)
#  rej1 = 1/0 indicator of rejection under H1 (NA = test not computable)
# The two scenarios are summarised separately, so they can have different numbers
# of valid replications.
summ_test <- function(rej0, rej1, name, alpha = 0.05) {
  v0  <- !is.na(rej0)
  v1  <- !is.na(rej1)
  n0  <- sum(v0)
  n1  <- sum(v1)
  est <- 100 * mean(rej0[v0])                                  # estimated type I error (%)
  ci  <- wilson_interval(alpha * n0, n0, conf = 0.95) * 100    # interval around the nominal level
  data.frame(test = name,
             type_I_error_pct = est,
             nominal_ci_lower_pct = unname(ci["lower"]),
             nominal_ci_upper_pct = unname(ci["upper"]),
             # within_nominal_ci = est >= ci["lower"] & est <= ci["upper"],
             power_pct = 100 * mean(rej1[v1]),
             n_valid_H0 = n0,
             n_valid_H1 = n1)
}

# Column of a data frame, or NA if the column does not exist
getcol <- function(d, col) if (col %in% names(d)) d[[col]] else rep(NA_real_, nrow(d))

# Tests reported for each scenario: label shown in the app and rejection column
# returned by simulate_one()
test_defs <- function(outcome, post) {
  if (identical(outcome, "continuous")) {
    d <- list(c("Naive t-test", "rejectsNAIVET"),
              c("F-test (specification 1)", "F1_rejectsWALD"),
              c("Permutation test (F1)", "F1_rejectsPERM"))
    if (identical(post, "ORIENT"))
      d <- c(d, list(c("F-test (specification 2)", "F2_rejectsWALD"),
                     c("Permutation test (F2)", "F2_rejectsPERM"),
                     c("F-test (specification 3)", "F3_rejectsWALD"),
                     c("Permutation test (F3)", "F3_rejectsPERM")))
  } else {
    d <- list(c("Unadjusted score test (asymptotic)", "rejects_unadj_asympt"),
              c("Unadjusted score test (randomization)", "rejects_unadj_rand"),
              c("Adjusted score test (asymptotic)", "rejects_adj_asympt"),
              c("Adjusted score test (randomization)", "rejects_adj_rand"))
  }
  data.frame(label  = vapply(d, function(x) x[1], character(1)),
             column = vapply(d, function(x) x[2], character(1)))
}

# Parse an allocation ratio written "a:b" (treated:control) into the number a/b (NA if invalid)
parse_ratio <- function(txt) {
  x <- strsplit(gsub("[[:space:]]", "", txt %||% ""), ":")[[1]]
  if (length(x) != 2) return(NA_real_)
  v <- suppressWarnings(as.numeric(x))
  if (anyNA(v) || any(v <= 0)) return(NA_real_)
  v[1] / v[2]
}

# Parse a number or a simple arithmetic expression such as "log(2)/4" (NA if invalid).
# Only digits, + - * / ^ ( ) and log/exp/sqrt are accepted.
parse_hazard <- function(txt) {
  txt <- trimws(txt %||% "")
  if (!nzchar(txt)) return(NA_real_)
  rest <- gsub("log|exp|sqrt", "", txt)
  if (!grepl("^[0-9eE.+*/^() -]*$", rest)) return(NA_real_)
  v <- tryCatch(eval(parse(text = txt), envir = baseenv()), error = function(e) NA_real_)
  if (length(v) != 1 || !is.numeric(v)) return(NA_real_)
  as.numeric(v)
}

# Parse "0.5, 0.3 0.2" into c(0.5, 0.3, 0.2) (NA if a value is not a number)
parse_probs <- function(txt) {
  x <- strsplit(trimws(gsub("[,;]", " ", txt)), "[[:space:]]+")[[1]]
  x <- x[nzchar(x)]
  suppressWarnings(as.numeric(x))
}

# Off-diagonal treatment effects of the orientation table, by default the effect of the
# original alteration (row i = eff[i]); the diagonal is always taken from the main table
default_tau_off <- function(eff) matrix(eff, nrow = length(eff), ncol = length(eff))

# Orientation table shown in the app: tau[original alteration, assigned stratum]
tau_df <- function(eff, m) {
  M <- length(eff)
  diag(m) <- eff
  df <- data.frame(paste("Alteration", seq_len(M)), m, check.names = FALSE)
  names(df) <- c("Original alteration (M)", paste("Assigned to stratum", seq_len(M)))
  df
}

# Orientation probabilities: M x M matrix, 0 on the diagonal, equal probabilities elsewhere
# (each row sums to 1; for M = 1 the matrix is just 0)
default_probs_orient <- function(M) {
  m <- matrix(if (M > 1) 1 / (M - 1) else 0, nrow = M, ncol = M)
  diag(m) <- 0
  m
}

# Data frame shown in the app: probs[original alteration, assigned stratum]
probs_df <- function(m) {
  M <- nrow(m)
  df <- data.frame(paste("Alteration", seq_len(M)), m, check.names = FALSE)
  names(df) <- c("Original alteration (M)", paste("Assigned to stratum", seq_len(M)))
  df
}

# Matrix of alteration x assigned-stratum treatment effects (tau[original, assigned]).
# DROP: patients stay in their own stratum, so only the diagonal is used.
# ORIENT: `off` (M x M) holds the effect of each possible orientation; its diagonal is replaced
# by the treatment effects of the main table.
make_tau <- function(eff, M, off = NULL) {
  tau <- if (is.null(off)) matrix(0, nrow = M, ncol = M) else off
  diag(tau) <- eff
  tau
}

# Arguments of simulate_one() for the scenario chosen in the app (list `p`).
# NULL elements are kept on purpose: simulate_one() requires some arguments to be NULL.
# build_args <- function(p) {
#   cont   <- identical(p$type_outcome, "continuous")
#   ext    <- identical(p$type_adapt, "EXT")
#   orient <- identical(p$stageII_approach, "ORIENT")
#   list(
#     type_outcome     = p$type_outcome,
#     type_rand        = if (cont) "block" else "minim",
#     type_adapt       = p$type_adapt,
#     stageII_approach = p$stageII_approach,
#     N = p$N, N1 = p$N1, M = p$M, prevalence = p$prevalence,
#     trt_labels = c(1, 2),
#     probs = c(1, p$probs) / (1 + p$probs),
#     bs           = if (cont) p$bs else NULL,
#     p_min        = if (cont) NULL else p$p_min,
#     adapt_time = if (cont) NULL else p$adapt_time,   # adapt_time in the app
#     Tr           = if (cont) NULL else p$Tr,
#     Ts           = if (cont) NULL else p$Ts,
#     delta        = if (cont) NULL else p$delta,
#     alpha        = if (cont) p$alpha else NULL,
#     lambda0      = if (cont) NULL else p$lambda0,    # baseline hazard (exponential: median = log(2) / lambda0)
#     sigma = p$sigma,
#     coefX        = if (cont) NULL else p$coefX,
#     probX        = if (cont) NULL else p$probX,
#     drop_trt     = if (ext) p$drop_trt else NULL,
#     z_c          = if (ext) NULL else p$z_c,
#     probs_orient = if (orient) p$probs_orient else NULL,
#     L = p$L,
#     max_atts     = if (ext) NULL else p$max_atts
#   )
# }

simulate_args <- c("type_outcome", "type_rand", "type_adapt", "stageII_approach",
                   "N", "N1", "M", "prevalence", "probs", "bs", "p_min", "adapt_time",
                   "Tr", "Ts", "delta", "alpha", "lambda0", "sigma", "coefX", "probX",
                   "drop_trt", "z_c", "probs_orient", "L", "max_atts")

build_args <- function(p) {
  a <- lapply(simulate_args, function(n) p[[n]])   # NULL when absent (required by simulate_one)
  names(a) <- simulate_args
  a$trt_labels <- c(1, 2)
  a                                                # tau is added in run_simulation (H0 / H1)
}

# ------------------------------------------------------------------
# 1. R CODE SHOWN IN THE APP
#    Builds a stand-alone script: the content of aux_functions.R followed by the
#    code that runs simulate_one() under H0 and H1 with the parameters of a given
#    run, and summarises the results exactly as the app does.
# ------------------------------------------------------------------
build_code <- function(p, tau_H1, functions_file = aux_file) {
  dp <- function(x) paste(deparse(x), collapse = "")
  M    <- p$M
  args <- build_args(p)
  args_txt <- paste0("  ", names(args), " = ", vapply(args, dp, character(1)), collapse = ",\n")
  pkgs     <- worker_pkgs(p)
  pkgs_txt <- paste0("c(", paste0('"', pkgs, '"', collapse = ", "), ")")
  td       <- test_defs(p$type_outcome, p$stageII_approach)
  tests_txt <- paste0(
    sprintf('  summ_test(rslt_H0$%s, rslt_H1$%s, "%s")', td$column, td$column, td$label),
    c(rep(",", nrow(td) - 1), ""))
  src <- readLines(functions_file, warn = FALSE)
  
  paste(c(
    "# R code reproducing the results shown in the app.",
    "# It contains the functions of aux_functions.R, followed by the parameters chosen",
    "# in the app and the code that runs the simulation under H0 and H1.",
    paste0("library(", c(pkgs, "future", "furrr"), ")"),
    "",
    "## ---- Functions (aux_functions.R) ----",
    src,
    "",
    "## ---- Wilson interval (centred on the nominal 5%) and performance measures of one test (in %) ----",
    paste0("wilson_interval <- ", paste(deparse(wilson_interval), collapse = "\n")),
    "",
    paste0("summ_test <- ", paste(deparse(summ_test), collapse = "\n")),
    "",
    "## ---- Parameters chosen in the app (arguments of simulate_one()) ----",
    paste0("Nsim <- ", dp(p$Nsim)),
    paste0("seed <- ", dp(p$seed)),
    paste0("M <- ", dp(M)),
    "base_args <- list(",
    args_txt,
    ")",
    "tau_H0 <- matrix(0, nrow = M, ncol = M)                 # no treatment effect",
    paste0("tau_H1 <- matrix(", dp(as.vector(tau_H1)), ", nrow = M, ncol = M)   # tau[original, assigned]"),
    "",
    "## ---- Running the simulation: Nsim replications under H0 and under H1 ----",
    "# The same seed is used for both scenarios, so replication i uses the same random",
    "# numbers under H0 and H1: the two scenarios differ only in the treatment effect.",
    "run_scenario <- function(tau) {",
    "  res <- future_map(seq_len(Nsim), function(i) {",
    "    tryCatch(do.call(simulate_one, c(base_args, list(tau = tau))),",
    "             error = function(e) NULL)                   # failed replications are discarded",
    paste0("  }, .options = furrr_options(seed = seed, packages = ", pkgs_txt, "))"),
    "  do.call(rbind, Filter(is.data.frame, res))",
    "}",
    "",
    if (p$workers > 1) paste0("plan(multisession, workers = ", p$workers, ")") else "plan(sequential)",
    "rslt_H0 <- run_scenario(tau_H0)   # type I error",
    "rslt_H1 <- run_scenario(tau_H1)   # power",
    "plan(sequential)",
    "",
    "## ---- Type I error and power, in %; Wilson 95% interval centred on the nominal 5% ----",
    "summary_df <- rbind(",
    tests_txt,
    ")",
    "print(summary_df)",
    "",
    "## ---- % of simulations in which each stratum is dropped (or none); each column sums to 100 ----",
    "drop_pct <- function(x) 100 * prop.table(table(factor(x,",
    "                                          levels = c(seq_len(M), 0),",
    '                                          labels = c(paste("Alteration", seq_len(M)), "No drop"))))',
    "drop_df <- rbind(H0 = drop_pct(rslt_H0$drop_trtIA), H1 = drop_pct(rslt_H1$drop_trtIA))",
    "print(round(drop_df, 2))",
    'barplot(drop_df, beside = TRUE, ylim = c(0, 105), ylab = "Simulations (%)", main = "Dropped stratum",',
    '        legend.text = TRUE, args.legend = list(x = "topright", bty = "n"))'
  ), collapse = "\n")
}

# ------------------------------------------------------------------
# 2. METHODS TEXT (shown in the app and saved in the .rds file)
# ------------------------------------------------------------------
methods_content <- function(p) {
  cont   <- identical(p$type_outcome, "continuous")
  ext    <- identical(p$type_adapt, "EXT")
  orient <- identical(p$stageII_approach, "ORIENT")
  
  list(
    # ---------------------------------------------------------------
    list(title = "Design, randomization, adaptation and treatment of new subjects",
         text = c(
           "In the simulation study, 'N' patients are planned to be enrolled. Patients are initially grouped into different strata, each corresponding to a specific molecular alteration ('M'). Each alteration has a certain prevalence ('f'), prognosis ('sigma') and treatment effect ('tau').",
           
           if (cont) {
             "Patients are randomized to the matching experimental treatment or to the strandad-of-care through a stratified permuted block randomization within each stratum, taking the allocation ratio 'treated:control' into account, with blocks of size 'b'."
           } else {
             "Patients are randomized to the matching experimental treatment or to the strandad-of-care through a minimization-based randomization (range method) on the original stratum and on one additional binary covariate 'X' (of effect on the log-hazard equal to 'coefX' and probability of X=1 equal to 'probX'), taking the allocation ratio 'treated:control' into account, with probability 'p_min' of assigning the arm that minimizes the imbalance."
           },
           
           if (ext) {
             "After the enrollment of 'N1' patients, one fixed treatment (and the corresponding stratum) is excluded due to an external event (EXT)."
           } else if (cont) {
             "After the enrollment of 'N1' patients, an interim futility analysis (IA) is performed to exclude at least one futile treatment (and the corresponding stratum). Under continuous outcomes, in each stratum the standardized difference in treatment effect is computed; the stratum with the smallest value is dropped if it is at most 'z_c'."
           } else {
             "After the enrollment of 'N1' patients, an interim futility analysis (IA) is performed to exclude at least one futile treatment (and the corresponding stratum). Under time-to-event outcomes, in each stratum a Cox model stratified by the covariate 'X' is fitted and the Z-score for the log-hazard ratio is computed; the stratum with the largest Z-score (worst effect) is dropped if it is above z_c."
           },
           
           if (!orient) {
             "The adaptation determines two stages of the trial. In case of adaptation, stage I patients belonging to the dropped strata are removed from the analysis (new stage I sample size 'N1_kept') and 'N2' = N-N1_kept new patients are enolled at stage II from the remaining alterations only (DROP), with prevalence rescaled to sum to 1. Treatment effects 'tau' are unchanged."
           } else {
             "The adaptation determines two stages of the trial. In case of adaptation, stage I patients belonging to the dropped strata are removed from the analysis (new stage I sample size 'N1_kept') and 'N2' = N-N1_kept new patients are enolled at stage II from all the original alterations (ORIENT) with the original prevalence; those belonging to the dropped alteration are oriented to one of the remaining strata with probability given by the row of the matrix 'probs_orient' corresponding to the dropped alteration (the diagonal is 0) by block randomization. The treatment effect 'tau' of oriented patients changes according to the stratum they are oriented into."
             },
           
           "At stage II, new patients are randomized to the matching experimental treatment or to the standard-of-care with the same randomization procedure used in stage I."
         )),
    
    # ---------------------------------------------------------------
    list(title = "Data generation",
         text = c(
           if (cont) {
             "Continuous outcome are generated through a regression model taking into account the control effect 'alpha', the original prognostic effects 'sigma', and the matrix of treatment effects 'tau', where each row corresponds to the original molecular alteration 'M' and each column to the assigned strata."
           } else {
             "Survival outcome are generated by considering exponential event times with baseline hazard 'lambda0', and the effects on the log-hazard scale of the covariate X ('coefX'), the original prognosis ('sigma'), and the matrix of treatments 'tau', where each row corresponds to the original molecular alteration 'M' and each column to the assignment strata. Patients enter uniformly over (0, 'adapt_time') in stage I and over ('adapt_time', 'Tr') in stage II, where 'Tr' is the last recruitment time. Follow-up continues util the time of total study duration 'Ts'. If patients are still alive at 'adapt_time' and Ts', they are right-censored (administrative censoring). Random censorings are exponentially distributed with rate 'delta'."
           }
         )),
    
    # ---------------------------------------------------------------
    list(title = "Tests",
         text = c(
           "All tests are two-sided at the 5% nominal level. The tested global null hypothesis 'H0' corresponds to the absence of an overall treatment effect, calculated as the weighted sum of stratum-specific treatment effects.",
           
           if (cont) {
             c(
               "Under continuous outcomes the following test statistics are evaluated:",
               "- Naive unstratified t-test, ignoring the strata;",
               if (!orient) {
                 "- F-test for linear contrasts (specification 1) under both asymptotic and permutation frameworks."
               }
               else {
                 c(
                   "- F-test for linear contrasts (specification 1), which ignores possible shifts in the prognostic and treatment effects due to the orientation of patients;",
                   "- F-test for linear contrasts (specifications 2), which accounts only for possible shifts in the prognostic effects due to the orientation of patients;",
                   "- F-test for linear contrasts (specifications 3), which accounts for both possible shifts in the prognostic and treatment effects due to the orientation of patients.",
                   "All these F-tests are evaluated under both asymptotic and permutation frameworks."
                 )
               },
               "Permutation tests are performed by randomly permuting the treatment assignment labels within the randomization blocks.",
               if (!ext) {
                 "Only the permutations leading to the same interim decision as the observed data are kept, and permutations are repeated until 'L' valid ones are obtained or 'max_atts' attempts are reached."
               }
             )
           }
           else {
             c(
               "Under time-to-event outcomes the following test statistics are evaluated:",
               "- Unadjusted score test, ignoring the minimization factors (equivalent to the log-rank test), under both asymptotic and randomization frameworks;",
               "- Score test adjusted for the minimization factors, under both asymptotic and randomization frameworks.",
               "Randomization tests are performed by re-assigning treatment labels through the minimization procedure, thereby preserving the accrual times.",
               if (!ext) {
                 "Only the permutations leading to the same interim decision as the observed data are kept, and permutations are repeated until 'L' valid ones are obtained or 'max_atts' attempts are reached."
               }
             )
           }
         )),
    
    # ---------------------------------------------------------------
    list(title = "Performance measures",
         text = c(
           if (cont) {
             c(
               "P-values for the asymptotic and permutation tests are computed ad each simulation step, and the following measures are compared:",
               "- Type I error: percentage of the valid simlations in which the test rejects H0 when is true. The 95% Wilson interval centred on the nominal level (5%) is reported: if the estimate of the type I error falls inside the interval, it is compatible with the nominal 5% given the Monte Carlo error; if it falls outside, the test is conservative (below the interval) or liberal (above it).", 
               "- Statistical power: percentage of valid simlations in which the test rejects H0 when is false.",
               "The p-value for the permutation test is estimated through the bilateral Monte Carlo estimator, as the proportion of permutations yielding a test statistic which is at least as extreme as the observed one, in absolute value."
             )
           }
           else {
             c(
               "P-values for the asymptotic and randomization tests are computed ad each simulation step, and the following measures are compared:",
               "- Type I error: percentage of the valid simulations in which the test rejects H0 when is true. The 95% Wilson interval centred on the nominal level (5%) is reported: if the estimate of the type I error falls inside the interval, it is compatible with the nominal 5% given the Monte Carlo error; if it falls outside, the test is conservative (below the interval) or liberal (above it).", 
               "- Statistical power: percentage of valid simulations in which the test rejects H0 when is false.",
               "The p-value for the randomization test is estimated through the bilateral Monte Carlo estimator, as the proportion of randomizations yielding a test statistic which is at least as extreme as the observed one, in absolute value."
             )
           }
         ))
  )
}

methods_ui <- function(p) {
  tagList(lapply(methods_content(p), function(s) {
    tagList(h5(s$title), lapply(s$text, function(t) tags$p(t)))
  }))
}

methods_text <- function(p) {
  paste(unlist(lapply(methods_content(p), function(s) c(toupper(s$title), s$text, ""))), collapse = "\n")
}

# ------------------------------------------------------------------
# 3. SIMULATION WRAPPER
#    Maps the UI parameters (list `p`) onto simulate_one(), repeats it
#    Nsim times in parallel under H0 and H1 and returns type I error / power
#    of each test, the % of simulations dropping each stratum, the raw
#    replications, the parameters used and the R code that reproduces the results.
# ------------------------------------------------------------------
run_simulation <- function(p) {
  
  M      <- p$M
  orient <- identical(p$stageII_approach, "ORIENT")
  
  p$prevalence <- p$prevalence / sum(p$prevalence)
  # if (orient) p$probs_orient <- p$probs_orient / sum(p$probs_orient)
  # if (orient && p$M >= 2) p$probs_orient <- prop.table(p$probs_orient, 1)
  
  tau_H0 <- matrix(0, nrow = M, ncol = M)
  # tau_H1 <- make_tau(p$tau, M, if (orient) p$tau_off else NULL)
  tau_H1 <- p$tau
  
  base_args <- build_args(p)
  
  # One replication; returns list(result = tibble or error message, warnings = character vector)
  one_rep <- function(i, tau) {
    warns <- character(0)
    # warnings shown in the app (the others only go to the console)
    shown <- c("blocksize is greater than the number of subjects",
               "permutations obtained out of")
    out <- withCallingHandlers(
      tryCatch(do.call(simulate_one, c(base_args, list(tau = tau))),
               error = function(e) conditionMessage(e)),
      warning = function(w) {
        msg <- conditionMessage(w)
        if (any(vapply(shown, function(s) grepl(s, msg, fixed = TRUE), logical(1))))
          warns <<- c(warns, msg)
        # no muffleWarning: the warning keeps going to the console
      })
    list(result = out, warnings = unique(warns))
  }
  
  # Parallel execution setup
  old_plan <- if (p$workers > 1) future::plan(future::multisession, workers = p$workers) else future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)
  
  # Same seed under H0 and H1: replication i uses the same random numbers in both
  # scenarios, which differ only in tau
  opts <- furrr::furrr_options(seed = p$seed, packages = worker_pkgs(p))
  res_h0 <- furrr::future_map(seq_len(p$Nsim), one_rep, tau = tau_H0, .options = opts)
  res_h1 <- furrr::future_map(seq_len(p$Nsim), one_rep, tau = tau_H1, .options = opts)
  
  # Separate the warnings from the results
  w_h0   <- lapply(res_h0, `[[`, "warnings")
  w_h1   <- lapply(res_h1, `[[`, "warnings")
  res_h0 <- lapply(res_h0, `[[`, "result")
  res_h1 <- lapply(res_h1, `[[`, "result")
  
  ok_h0 <- vapply(res_h0, is.data.frame, logical(1))
  ok_h1 <- vapply(res_h1, is.data.frame, logical(1))
  
  if (!any(ok_h0)) stop("All H0 replications failed. First error: ", res_h0[[1]])
  if (!any(ok_h1)) stop("All H1 replications failed. First error: ", res_h1[[1]])
  
  r_h0 <- as.data.frame(do.call(rbind, res_h0[ok_h0]))
  r_h1 <- as.data.frame(do.call(rbind, res_h1[ok_h1]))
  
  # Type I error from the H0 replications, power from the H1 replications
  td <- test_defs(p$type_outcome, p$stageII_approach)
  summary_df <- do.call(rbind, lapply(seq_len(nrow(td)), function(i)
    summ_test(getcol(r_h0, td$column[i]), getcol(r_h1, td$column[i]), td$label[i])))
  
  # drop_trtIA = stratum dropped in the replication (0 = none)
  drop_pct <- function(x) {
    cnt <- table(factor(x, levels = c(seq_len(M), 0),
                        labels = c(paste("Alteration", seq_len(M)), "No drop")))
    100 * as.numeric(cnt) / sum(cnt)
  }
  drop_df <- data.frame(category = c(paste("Alteration", seq_len(M)), "No drop"),
                        pct_H0 = drop_pct(r_h0$drop_trtIA),
                        pct_H1 = drop_pct(r_h1$drop_trtIA))
  
  # IA: mean number of attempts needed to obtain L valid permutations / randomizations
  attempts <- if (identical(p$type_adapt, "IA"))
    c(H0 = mean(r_h0$Nattempts, na.rm = TRUE), H1 = mean(r_h1$Nattempts, na.rm = TRUE)) else NULL
  
  # Warnings: number of replications in which each message was raised
  norm_warn <- function(x) gsub("Only +[0-9]+ +permutations", "Only k permutations", x)
  warn_tab <- function(w, scen) {
    u <- norm_warn(unlist(w))
    if (!length(u)) return(NULL)
    tb <- table(u)
    data.frame(scenario = scen, message = names(tb), n = as.integer(tb),
               stringsAsFactors = FALSE)
  }
  warnings_df <- rbind(warn_tab(w_h0, "H0"), warn_tab(w_h1, "H1"))
  
  list(summary = summary_df, drop = drop_df, attempts = attempts,
       replications = list(H0 = r_h0, H1 = r_h1),
       params = p, tau_H1 = tau_H1, code = build_code(p, tau_H1),
       n_failed = sum(!ok_h0) + sum(!ok_h1),
       warnings = warnings_df)
}

# Grouped bar chart of the % of simulations dropping each stratum (or none)
plot_drop <- function(d) {
  op <- par(mar = c(7, 5, 4, 1))
  on.exit(par(op))
  m <- rbind(H0 = d$pct_H0, H1 = d$pct_H1)
  b <- barplot(m, beside = TRUE, names.arg = d$category, ylim = c(0, 112),
               col = c("steelblue", "orange"), las = 2, ylab = "Simulations (%)",
               main = "Dropped stratum across simulations",
               legend.text = c("Under H0", "Under H1"),
               args.legend = list(x = "topright", bty = "n"))
  text(b, m, sprintf("%.1f", m), pos = 3, cex = 0.7)
}

# ------------------------------------------------------------------
# 4. PARAMETER CHECKS
#    Returns a character vector with one message per problem found
#    (empty vector = everything is fine).
# ------------------------------------------------------------------
validate_params <- function(p) {
  e <- character(0)
  add <- function(msg) if (!is.null(msg)) e <<- c(e, msg)
  
  num1 <- function(x) length(x) == 1 && is.numeric(x) && !is.na(x) && is.finite(x)
  
  # message if x is not a single positive integer (NULL otherwise)
  pos_int_msg <- function(x, name) {
    if (!num1(x))         paste0(name, " must be a number.")
    else if (x <= 0)      paste0(name, " must be positive.")
    else if (x != round(x)) paste0(name, " must be an integer.")
    else NULL
  }
  
  cont <- identical(p$type_outcome, "continuous")
  surv <- identical(p$type_outcome, "survival")
  ext  <- identical(p$type_adapt, "EXT")
  
  # ---- randomization method ----
  if ((cont && !identical(p$type_rand, "block")) || (surv && !identical(p$type_rand, "minim")))
    add("Only stratified block randomization (continuous outcome) and minimization (survival outcome) are available.")
  probs_ok <- length(p$probs) == 2 && all(is.finite(p$probs)) && all(p$probs > 0)
  if (!probs_ok)
    add("The allocation ratio must be written as a:b with two positive numbers (e.g. 1:1 or 2:1).")
  
  # ---- N and N1 ----
  m <- pos_int_msg(p$N, "N")
  add(m)
  if (is.null(m) && p$N < 2) add("N must be at least 2.")
  m1 <- pos_int_msg(p$N1, "N1")
  add(m1)
  if (is.null(m1) && num1(p$N) && p$N1 >= p$N) add("N1 must be smaller than N.")
  if (is.null(m1) && p$N1 < 2) add("N1 must be at least 2.")
  
  # ---- M (= number of rows of the table) and prevalence ----
  M_ok <- num1(p$M) && p$M >= min_M
  if (!M_ok) add(sprintf("The table of alterations needs at least %d rows.", min_M))
  if (M_ok) {
    f <- p$prevalence
    if (!is.numeric(f) || length(f) != p$M || anyNA(f)) {
      add("All prevalence must be numbers.")
    } else {
      if (any(f <= 0 | f >= 1)) add("Each prevalence must be strictly between 0 and 1.")
      if (abs(sum(f) - 1) > 1e-3)
        add(sprintf("prevalence must sum to 1 (currently %.4f).", sum(f)))
    }
    if (!is.numeric(p$sigma) || length(p$sigma) != p$M || anyNA(p$sigma))
      add("All prognosis values (sigma) in the table must be numbers.")
    # if (!is.numeric(p$tau) || length(p$tau) != p$M || anyNA(p$tau))
    #   add("All treatment effects (under H1) in the table must be numbers.")
    if (!is.matrix(p$tau) || any(dim(p$tau) != p$M) || anyNA(p$tau))
      add("All treatment effects in the tables must be numbers.")
  }
  
  # ---- randomization ----
  if (cont) {
    m <- pos_int_msg(p$bs, "Block size")
    add(m)
    if (is.null(m) && probs_ok) {
      expected <- p$bs * p$probs        # probs is already c(control, treated)
      if (!all(abs(expected - round(expected)) < 1e-3))
        add("Block size must be proportional to the allocation ratio (e.g. a multiple of 2 for 1:1, of 3 for 2:1).")
    }
  }
  if (surv) {
    if (!num1(p$p_min) || p$p_min <= 0.5 || p$p_min >= 1)
      add("The probability of assignment (minimization) must be strictly between 0.5 and 1.")
    if (!num1(p$coefX))                       add("The effect of the covariate X must be a number.")
    if (!num1(p$probX) || p$probX <= 0 || p$probX >= 1)
      add("The probability of X = 1 must be strictly between 0 and 1.")
  }
  
  # ---- outcome-specific parameters ----
  if (cont && !num1(p$alpha)) add("The control effect (alpha) must be a number.")
  if (surv) {
    if (!num1(p$adapt_time) || p$adapt_time <= 0) add("The adaptation time must be positive.")
    if (!num1(p$Tr) || p$Tr <= 0)   add("The accrual period must be positive.")
    if (!num1(p$Ts) || p$Ts <= 0) add("The study duration must be positive.")
    if (num1(p$adapt_time) && num1(p$Tr) && p$adapt_time >= p$Tr)
      add("The adaptation time must be before the end of accrual.")
    if (num1(p$Tr) && num1(p$Ts) && p$Tr > p$Ts)
      add("The accrual period cannot be longer than the study duration.")
    if (!num1(p$lambda0) || p$lambda0 <= 0)         add("The baseline hazard must be a positive number or expression (e.g. log(2)/4).")
    if (!num1(p$delta) || p$delta < 0)        add("The random censoring rate (delta) must be >= 0.")
  }
  
  # ---- adaptation ----
  if (ext && M_ok) {
    if (!num1(p$drop_trt) || p$drop_trt != round(p$drop_trt) || p$drop_trt < 1 || p$drop_trt > p$M)
      add("The category to exclude must be an integer between 1 and M.")
  }
  if (!ext) {
    if (!num1(p$z_c)) add("The futility threshold (z_c) must be a number.")
    m <- pos_int_msg(p$max_atts, "max_atts")
    add(m)
    if (is.null(m) && num1(p$L) && p$max_atts < p$L) add("max_atts must be >= L.")
  }
  
  # ---- new subjects after adaptation ----
  if (identical(p$stageII_approach, "ORIENT") && M_ok) {
    po <- p$probs_orient
    if (!is.matrix(po) || !is.numeric(po) || any(dim(po) != p$M) || anyNA(po)) {
      add("All orientation probabilities in the table must be numbers.")
    } else {
      if (any(diag(po) != 0)) add("The diagonal of the orientation probabilities table must be 0.")
      off <- po[row(po) != col(po)]
      if (any(off < 0 | off > 1)) add("Orientation probabilities must be between 0 and 1.")
      # if (p$M >= 2) {
      #   rs <- rowSums(po)
      #   bad <- which(abs(rs - 1) > 1e-3)
      #   if (length(bad))
      #     add(sprintf("Each row of the orientation probabilities table must sum to 1 (check: %s).",
      #                 paste("Alteration", bad, collapse = ", ")))
      # }
    }
  }
  
  # ---- simulation size, seed, workers ----
  add(pos_int_msg(p$Nsim, "Nsim"))
  add(pos_int_msg(p$L, "L"))
  if (!num1(p$seed) || p$seed != round(p$seed)) add("The seed must be an integer.")
  m <- pos_int_msg(p$workers, "Workers")
  add(m)
  if (is.null(m) && p$workers > max_workers)
    add(sprintf("Workers must be lower than the number of available cores (%d): the maximum is %d.",
                n_cores, max_workers))
  
  e
}

# ------------------------------------------------------------------
# 5. UI
# ------------------------------------------------------------------
warn_box <- function(...) div(class = "alert alert-warning py-1 px-2 small", ...)

ui <- page_sidebar(
  title = "Adaptive stratified trial: simulation study",
  theme = bs_theme(version = 5),
  
  sidebar = sidebar(
    width = 440,
    accordion(
      open = c("1. Outcome", "2. Design"),
      multiple = TRUE,
      
      # ---- 1. Outcome ------------------------------------------------
      accordion_panel(
        "1. Outcome",
        radioButtons("outcome", "Outcome type",
                     choices = c("Continuous" = "continuous", "Survival" = "survival"),
                     inline = TRUE),
        warn_box("Warning: only the outcomes types listed above are currently available."),
      ),
      
      # ---- 2. Design -------------------------------------------------
      accordion_panel(
        "2. Design",
        numericInput("N", "Total number of subjects (N)", value = 200, min = 2, step = 10),
        strong("Molecular alterations (one row each)"),
        DT::DTOutput("mtab"),
        div(class = "d-flex gap-2 mt-2 mb-1",
            actionButton("add_row", "Add row", class = "btn-sm"),
            actionButton("del_row", "Remove last row", class = "btn-sm"),
            actionButton("eq_prev", "Equal prevalence", class = "btn-sm")),
        # uiOutput("freq_sum"),
        helpText("M is the number of rows. Double-click a cell to edit it. ",
                 "prevalence must sum to 1. Treatment effects under H0 are set to 0."),
        helpText("Continuous outcome: prognosis = shift of the mean outcome, treatment effect = mean difference. ",
                 "Survival outcome: prognosis = log-hazard shift, treatment effect = log-hazard ratio ",
                 "(negative = benefit). Changing the outcome type resets the treatment effects to a default value.")
      ),
      
      # ---- 3. Randomization ------------------------------------------
      accordion_panel(
        "3. Randomization",
        radioButtons("rand", "Randomization method",
                     choices = c("Stratified block randomization" = "block")),   # updated by the outcome
        warn_box("Warning: only the methods listed above are currently available."),
        textInput("ratio", "Allocation ratio (treated : control)", value = "1:1",
                  placeholder = "e.g. 1:1 or 2:1"),
        conditionalPanel(
          "input.rand == 'block'",
          numericInput("bs", "Block size (b)", value = 4, min = 2, step = 2)
        ),
        conditionalPanel(
          "input.rand == 'minim'",
          numericInput("p_min", "Probability of assignment (p_min)", min = 0.51, max = 0.99, value = 0.9, step = 0.05),
          warn_box("Choose 1 additional minimization factor (besides the molecular alteration): ",
                   "only one additional binary factor X is currently available."),
          numericInput("probX", "Probability of X = 1", min = 0.01, max = 0.99,
                       value = 0.4, step = 0.1),
          numericInput("coefX", "Effect on the log-hazard (coefX)", 
                       value = 0.3, step = 0.1)
        )
      ),
      
      # ---- 4. Outcome-specific parameters ----------------------------
      accordion_panel(
        "4. Outcome-specific parameters",
        conditionalPanel(
          "input.outcome == 'continuous'",
          numericInput("alpha_ctrl", "Control effect (alpha)", value = 1, step = 0.1),
        ),
        conditionalPanel(
          "input.outcome == 'survival'",
          textInput("lambda0", "Baseline hazard (log(2)/median)", value = "log(2)/4",
                    placeholder = "e.g. log(2)/4 or 0.1733"),
          uiOutput("lambda0_note"),
          numericInput("accrual", "Last recruitment time (Tr)", value = 72, min = 0.1, step = 1),
          numericInput("followup", "Total study duration (Ts)", value = 84, min = 0.1, step = 1),
          numericInput("delta", "Random censoring rate (delta)", value = 0.01, min = 0, step = 0.01)
        )
      ),
      
      # ---- 5. Adaptation ---------------------------------------------
      accordion_panel(
        "5. Adaptation",
        warn_box("Warning: only one adaptation is currently available."),
        numericInput("N1", "Subjects in the first stage (N1)", value = 100, min = 2, step = 10),
        conditionalPanel(
          "input.outcome == 'survival'",
          numericInput("adapt_time", "Adaptation time (adapt_time)", value = 36, min = 0.1, step = 1),
          helpText("Same time unit as the survival parameters (e.g. months). SHould be consistent with N1.")
        ),
        radioButtons("adapt", "Type of adaptation",
                     choices = c("EXT (external event)" = "EXT",
                                 "IA (interim futility analysis)" = "IA")),
        conditionalPanel(
          "input.adapt == 'EXT'",
          uiOutput("ext_ui")
        ),
        conditionalPanel(
          "input.adapt == 'IA'",
          numericInput("fut_thr", "Futility threshold (z_c)", value = 0.9, step = 0.05),
          helpText("Continuous: the stratum with the smallest standardized effect is dropped if it is <= z_c. ",
                   "Survival: the stratum with the largest Z-score for the logHR is dropped if it is > z_c.")
        )
      ),
      
      # ---- 6. New subjects after adaptation --------------------------
      accordion_panel(
        "6. New subjects after adaptation",
        radioButtons("post", "What to do with new subjects",
                     choices = c("DROP" = "DROP", "ORIENT" = "ORIENT"), inline = TRUE),
        helpText("DROP: new subjects are enrolled only from the remaining alterations. Treatment effects are supposed to remain the same.",
                 "ORIENT: new subjects of the dropped alteration are enrolled and oriented to the remaining strata. Treatment effects of oriented patients are supposed to change."),
        conditionalPanel(
          "input.post == 'ORIENT'",
          strong("Orientation probabilities (probs_orient)"),
          DT::DTOutput("probs_tab"),
          div(class = "mt-2 mb-3",
              actionButton("probs_reset", "Equal probabilities", class = "btn-sm")),
          helpText("Rows: original alteration of the patient. Columns: stratum the patient is oriented to. ",
                   "Each row gives the orientation probabilities used when that alteration is dropped: ",
                   "each row must sum to 1 (or to 0 if the alteration is never dropped). ",
                   "The diagonal is always 0 and cannot be edited."),
          strong("Treatment effect of each possible orientation under H1 (tau)"),
          DT::DTOutput("tau_tab"),
          div(class = "mt-2 mb-1",
              actionButton("tau_reset", "Fill with the effect of the original alteration", class = "btn-sm")),
          helpText("Rows: original alteration of the patient. Columns: stratum the patient is assigned to. ",
                   "The diagonal reports the treatment effect of the table of alterations (section 2) and cannot be ",
                   "edited here. Double-click the other cells to edit them.")
        )
      ),
      
      # ---- 7. Simulation size ----------------------------------------
      accordion_panel(
        "7. Simulation size",
        numericInput("Nsim", "Number of simulations (Nsim)", value = 10000, min = 1, step = 100),
        numericInput("L", "Number of permutations / randomizations (L)", value = 1000, min = 1, step = 100),
        helpText("The permutation / randomization tests are run in every replication, under H0 and H1: ",
                 "computing time grows roughly with Nsim x L, and it is much longer for the survival outcome and under the IA scenario. ",
                 "Try small values first."),
        conditionalPanel(
          "input.adapt == 'IA'",
          numericInput("max_atts", "Maximum number of attempts (max_atts)", value = 100000, min = 1, step = 500),
          helpText("Under the IA scenario, only permutations giving the same interim decision as the observed data are kept.")
        )
      ),
      
      # ---- 8. Reproducibility / computing ----------------------------
      accordion_panel(
        "8. Reproducibility and computing",
        numericInput("seed", "Seed (for reproducibility)", value = 123, step = 1),
        numericInput("workers", "Workers (for parallelization)",
                     value = min(max_workers, max(1L, n_cores - 2L)), min = 1,
                     max = max_workers, step = 1),
        helpText(sprintf("Available cores: %d. Workers must be lower than this number. The largest the number of workers, the fastest the code runs.", n_cores))
      )
    ),
    
    actionButton("run", "Run simulation", class = "btn-primary w-100"),
    uiOutput("dl_ui")      # download button, shown only once results exist
  ),
  
  navset_card_tab(
    selected = "Results",
    nav_panel("Methods", uiOutput("methods_out")),
    nav_panel("Parameter summary", verbatimTextOutput("params_txt")),
    nav_panel(
      "Results",
      uiOutput("errors"),
      uiOutput("scenario"),
      uiOutput("warns"),
      h5("Type I error and statistical power (nominal level alpha = 0.05)"),
      tableOutput("tab"),
      uiOutput("attempts_note")
    ),
    nav_panel(
      "Plot",
      plotOutput("drop_plot", width = "auto", height = "auto")   
    ),
    nav_panel(
      "Code",
      tags$button("Copy code", class = "btn btn-sm btn-outline-secondary mb-2",
                  onclick = "
                    var btn = this;
                    var txt = document.getElementById('code_text').textContent;
                    function done(ok) {
                      btn.innerText = ok ? 'Copied!' : 'Copy failed: select the code and press Ctrl+C';
                      setTimeout(function() { btn.innerText = 'Copy code'; }, 2500);
                    }
                    function legacy() {
                      var ta = document.createElement('textarea');
                      ta.value = txt;
                      ta.style.position = 'fixed';
                      ta.style.top = '0';
                      ta.style.left = '0';
                      ta.style.opacity = '0';
                      document.body.appendChild(ta);
                      ta.focus();
                      ta.select();
                      var ok = false;
                      try { ok = document.execCommand('copy'); } catch (e) {}
                      document.body.removeChild(ta);
                      return ok;
                    }
                    if (legacy()) {
                      done(true);
                    } else if (navigator.clipboard) {
                      navigator.clipboard.writeText(txt).then(function() { done(true); },
                                                              function() { done(false); });
                    } else {
                      done(false);
                    }
                  "),
      div(style = "max-height: 75vh; overflow: auto;", verbatimTextOutput("code_text"))
    )
  )
)

# ------------------------------------------------------------------
# 6. SERVER
# ------------------------------------------------------------------
server <- function(input, output, session) {
  
  # ---- Table of the molecular alterations: one row = one alteration, so M = nrow ----
  tab <- reactiveVal(data.frame(
    alteration = paste("Alteration", 1:4),
    prevalence = rep(0.25, 4),
    prognosis  = rep(0.5, 4),
    effect     = rep(0.4, 4)
  ))
  n_alt <- reactiveVal(4L)     # M: changes only when a row is added or removed
  
  # Orientation table: M x M treatment effects tau[original, assigned]; the diagonal is
  # taken from the main table, so only the off-diagonal cells are stored here
  tau_off <- reactiveVal(default_tau_off(rep(0.4, 4)))
  
  # Orientation probabilities: M x M, diagonal = 0 (not editable)
  probs_or <- reactiveVal(default_probs_orient(4))
  
  output$mtab <- DT::renderDT({
    DT::datatable(
      isolate(tab()), rownames = FALSE, selection = "none", class = "compact stripe",
      colnames = c("Alteration (M)", "Prevalence (f)", "Prognosis (sigma)", "Treatment effect under H1 (tau)"),
      editable = list(target = "cell", disable = list(columns = 0)),
      options = list(dom = "t", paging = FALSE, ordering = FALSE)
    )
  })
  proxy <- DT::dataTableProxy("mtab")
  
  # The orientation table is redrawn when M changes (n_alt) and updated in place otherwise
  output$tau_tab <- DT::renderDT({
    n_alt()
    DT::datatable(
      isolate(tau_df(tab()$effect, tau_off())), rownames = FALSE, selection = "none",
      class = "compact stripe",
      editable = list(target = "cell", disable = list(columns = 0)),
      options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE)
    )
  })
  outputOptions(output, "tau_tab", suspendWhenHidden = FALSE)
  proxy_tau <- DT::dataTableProxy("tau_tab")
  
  refresh_tau <- function() {
    eff <- tab()$effect
    m   <- tau_off()
    if (length(eff) == nrow(m))
      DT::replaceData(proxy_tau, tau_df(eff, m), resetPaging = FALSE, rownames = FALSE)
  }
  
  output$probs_tab <- DT::renderDT({
    n_alt()
    DT::datatable(
      isolate(probs_df(probs_or())), rownames = FALSE, selection = "none",
      class = "compact stripe",
      editable = list(target = "cell", disable = list(columns = 0)),
      options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE)
    )
  })
  outputOptions(output, "probs_tab", suspendWhenHidden = FALSE)
  proxy_probs <- DT::dataTableProxy("probs_tab")
  
  refresh_probs <- function() {
    DT::replaceData(proxy_probs, probs_df(probs_or()), resetPaging = FALSE, rownames = FALSE)
  }
  
  # Replace the table (labels are renumbered, M is updated)
  set_tab <- function(d) {
    d$alteration <- paste("Alteration", seq_len(nrow(d)))
    tab(d)
    n_alt(nrow(d))
    DT::replaceData(proxy, d, resetPaging = FALSE, rownames = FALSE)
  }
  
  # Edit of a cell (non-numeric text becomes NA and is reported as an error)
  observeEvent(input$mtab_cell_edit, {
    info <- input$mtab_cell_edit
    j <- info$col + 1                    # rownames = FALSE: column index is 0-based
    if (j < 2) return()                  # the label column cannot be edited
    d <- tab()
    d[info$row, j] <- suppressWarnings(as.numeric(info$value))
    tab(d)
    DT::replaceData(proxy, d, resetPaging = FALSE, rownames = FALSE)
    if (j == 4) refresh_tau()            # the diagonal of the orientation table follows the effects
  })
  
  # Edit of a cell of the orientation table (the diagonal cannot be edited)
  observeEvent(input$tau_tab_cell_edit, {
    info <- input$tau_tab_cell_edit
    i <- info$row
    j <- info$col                        # rownames = FALSE: 0 = label column, 1..M = assigned stratum
    m <- tau_off()
    if (j >= 1 && i != j && i <= nrow(m) && j <= ncol(m)) {
      m[i, j] <- suppressWarnings(as.numeric(info$value))
      tau_off(m)
    }
    refresh_tau()
  })
  
  observeEvent(input$tau_reset, {
    tau_off(default_tau_off(tab()$effect))
    refresh_tau()
  })
  
  # Edit of a cell of the orientation-probabilities table:
  # the diagonal cannot be edited; the value must be a number in [0, 1]
  observeEvent(input$probs_tab_cell_edit, {
    info <- input$probs_tab_cell_edit
    i <- info$row
    j <- info$col                        # rownames = FALSE: 0 = label column, 1..M = assigned stratum
    m <- probs_or()
    if (j >= 1 && i != j && i <= nrow(m) && j <= ncol(m)) {
      v <- suppressWarnings(as.numeric(info$value))
      if (is.na(v) || v < 0 || v > 1) {
        showNotification("Orientation probabilities must be numbers between 0 and 1 (inclusive).",
                         type = "warning")
      } else {
        m[i, j] <- v
        probs_or(m)
      }
    }
    refresh_probs()
  })
  
  observeEvent(input$probs_reset, {
    probs_or(default_probs_orient(n_alt()))
    refresh_probs()
  })
  
  observeEvent(input$add_row, {
    d <- tab()
    if (nrow(d) >= max_M) {
      showNotification(sprintf("At most %d alterations.", max_M), type = "warning")
      return()
    }
    last <- d[nrow(d), ]
    # new row: same prognosis and effect as the last one; the prevalence has to be set
    d <- rbind(d, data.frame(alteration = NA_character_, prevalence = 0,
                             prognosis = last$prognosis, effect = last$effect))
    # orientation table: keep the existing cells, new row and column filled with the default
    M1 <- nrow(d)
    m_new <- default_tau_off(d$effect)
    m_new[seq_len(M1 - 1), seq_len(M1 - 1)] <- tau_off()
    tau_off(m_new)
    probs_or(default_probs_orient(nrow(d)))
    set_tab(d)
  })
  
  observeEvent(input$del_row, {
    d <- tab()
    if (nrow(d) <= min_M) {
      showNotification(sprintf("At least %d alterations are needed.", min_M), type = "warning")
      return()
    }
    tau_off(tau_off()[-nrow(d), -nrow(d), drop = FALSE])
    probs_or(default_probs_orient(nrow(d) - 1L))
    set_tab(d[-nrow(d), ])
  })
  
  # Rounded to 4 decimals for display; the simulation rescales them to sum exactly to 1
  observeEvent(input$eq_prev, {
    d <- tab()
    d$prevalence <- round(1 / nrow(d), 4)
    set_tab(d)
  })
  
  # The treatment effect is a mean difference (continuous) or a log-hazard ratio (survival):
  # reset it to a sensible default when the outcome type changes
  observeEvent(input$outcome, ignoreInit = TRUE, {
    d <- tab()
    d$effect <- rep(if (identical(input$outcome, "survival")) -0.4 else 0.4, nrow(d))
    tau_off(default_tau_off(d$effect))
    set_tab(d)
    refresh_tau()
  })
  
  # The randomization methods offered follow the outcome
  # (only block for the continuous outcome and minimization for the survival outcome so far)
  observeEvent(input$outcome, {
    ch <- if (identical(input$outcome, "continuous")) c("Stratified block randomization" = "block")
    else c("Minimization (Range method)" = "minim")
    updateRadioButtons(session, "rand", choices = ch, selected = unname(ch)[1])
  })
  
  # # M and live sum of the prevalence
  # output$freq_sum <- renderUI({
  #   d <- tab()
  #   s <- sum(d$prevalence)
  #   ok <- !is.na(s) && abs(s - 1) <= 1e-3
  #   span(style = paste0("font-size: 0.85em; color: ", if (ok) "#2e7d32" else "#c62828", ";"),
  #        sprintf("M = %d alterations. ", nrow(d)),
  #        if (is.na(s)) "Sum of prevalence: NA (check the values)"
  #        else sprintf("Sum of prevalence: %.4f %s", s, if (ok) "(OK)" else "(must be 1)"))
  # })
  
  # Category to exclude
  output$ext_ui <- renderUI({
    M <- n_alt()
    sel <- suppressWarnings(as.integer(isolate(input$ext_cat)))
    sel <- if (length(sel) == 1 && !is.na(sel) && sel >= 1 && sel <= M) sel else M
    selectInput("ext_cat", "Alteration M to exclude",
                choices = setNames(seq_len(M), paste("Alteration", seq_len(M))),
                selected = sel)
  })
  
  # Collect all parameters in one list
  params <- reactive({
    d <- tab()
    M <- nrow(d)
    cont <- identical(input$outcome, "continuous")
    ratio  <- parse_ratio(input$ratio)
    orient <- identical(input$post, "ORIENT")
    off    <- if (orient) tau_off() else NULL
    tau    <- if (is.null(off) || nrow(off) == M) make_tau(d$effect, M, off)
    else matrix(NA_real_, M, M)       # transient state while M is changing
    p <- list(
      type_outcome = input$outcome,
      type_rand = input$rand,
      type_adapt = input$adapt,
      stageII_approach = input$post,
      N = input$N, 
      N1 = input$N1, 
      M = M,
      prevalence = d$prevalence,
      probs = c(1, ratio) / (1 + ratio),   
      sigma = d$prognosis,
      tau = tau,                                      # treatment effects under H1 
      Nsim = input$Nsim, 
      L = input$L,
      seed = input$seed, 
      workers = input$workers
    )
    if (cont) {
      p <- c(p, list(alpha = input$alpha_ctrl, 
                     bs = input$bs))
    } else {
      p <- c(p, list(p_min = input$p_min, 
                     coefX = input$coefX,
                     probX = input$probX,
                     adapt_time = input$adapt_time, 
                     lambda0 = parse_hazard(input$lambda0),
                     Tr = input$accrual, 
                     Ts = input$followup, 
                     delta = input$delta))
    }
    if (input$adapt == "EXT") {
      p$drop_trt <- suppressWarnings(as.integer(input$ext_cat %||% M))
    } else {
      p$z_c <- input$fut_thr
      p$max_atts <- input$max_atts
    }
    if (input$post == "ORIENT") {
      pm <- probs_or()
      p$probs_orient <- if (nrow(pm) == M) pm else matrix(NA_real_, M, M)   # transient state while M is changing
    }
    p
  })
  
  output$params_txt <- renderPrint(str(params()))
  
  output$methods_out <- renderUI(methods_ui(params()))
  
  # Value read from the baseline hazard field, with the corresponding median survival
  output$lambda0_note <- renderUI({
    v <- parse_hazard(input$lambda0)
    helpText(if (is.na(v) || v <= 0) "Enter a positive number or an expression such as log(2)/4."
             else sprintf("Baseline hazard = %.4f (median survival = %.2f).", v, log(2) / v))
  })
  
  # Live parameter checks
  errors <- reactive(validate_params(params()))
  
  output$errors <- renderUI({
    e <- errors()
    if (!length(e)) return(NULL)
    div(class = "alert alert-danger",
        strong("Please check the parameters:"),
        tags$ul(lapply(e, tags$li)))
  })
  
  # Results
  res <- reactiveVal(NULL)
  
  observeEvent(input$run, {
    e <- errors()
    if (length(e)) {
      showNotification("Please fix the parameter errors before running the simulation.",
                       type = "error", duration = 8)
      return()
    }
    p <- params()
    out <- tryCatch(
      withProgress(message = "Running simulation...", detail = "this can take several minutes",
                   value = 0.5, run_simulation(p)),
      error = function(err) {
        showNotification(paste("Simulation failed:", conditionMessage(err)),
                         type = "error", duration = NULL)
        NULL
      }
    )
    if (is.null(out)) return()
    if (out$n_failed > 0)
      showNotification(sprintf("%d of %d replications (H0 and H1) failed and were discarded.",
                               out$n_failed, 2L * p$Nsim),
                       type = "warning", duration = 10)
    res(out)
  })
  
  # Scenario of the displayed results
  output$scenario <- renderUI({
    r <- res()
    if (is.null(r)) return(NULL)
    p <- r$params
    tags$p(class = "text-muted",
           sprintf("Scenario: %s outcome (%s), adaptation %s, new subjects %s; N = %d, N1 = %d, M = %d, Nsim = %d, L = %d.",
                   p$type_outcome, if (identical(p$type_outcome, "continuous")) "block randomization" else "minimization",
                   p$type_adapt, p$stageII_approach, as.integer(p$N), as.integer(p$N1), p$M,
                   as.integer(p$Nsim), as.integer(p$L)))
  })
  
  # Warnings raised by simulate_one() during the simulation
  output$warns <- renderUI({
    r <- res()
    if (is.null(r) || is.null(r$warnings) || !nrow(r$warnings)) return(NULL)
    w <- r$warnings
    div(class = "alert alert-warning",
        strong("Warnings raised during the simulation:"),
        tags$ul(lapply(seq_len(nrow(w)), function(i)
          tags$li(sprintf("[%s] %s (in %d of %d replications)",
                          w$scenario[i], w$message[i], w$n[i], r$params$Nsim)))))
  })
  
  # Type I error (with Wilson interval) and power, in %
  output$tab <- renderTable({
    validate(need(!is.null(res()), "Set the parameters and press 'Run simulation'."))
    s <- res()$summary
    data.frame(
      "Test" = s$test,
      "Type I error" = sprintf("%.2f%%", s$type_I_error_pct),
      "Wilson 95% interval around 5%" = sprintf("[%.2f%%, %.2f%%]",
                                                s$nominal_ci_lower_pct, s$nominal_ci_upper_pct),
      "Power" = sprintf("%.2f%%", s$power_pct),
      "Valid replications (H0 / H1)" = sprintf("%d / %d", s$n_valid_H0, s$n_valid_H1),
      check.names = FALSE
    )
  }, align = "lcccc")
  
  # IA: attempts needed to obtain L valid permutations / randomizations
  output$attempts_note <- renderUI({
    r <- res()
    if (is.null(r) || is.null(r$attempts)) return(NULL)
    helpText(sprintf("Mean number of attempts to obtain L = %d valid permutations / randomizations (max_atts = %d): %.0f under H0, %.0f under H1.",
                     as.integer(r$params$L), as.integer(r$params$max_atts),
                     r$attempts[["H0"]], r$attempts[["H1"]]))
  })
  
  # ---- Tab "Plot": % of simulations dropping each stratum (fixed size) ----
  output$drop_plot <- renderPlot({
    validate(need(!is.null(res()), "Run the simulation first."))
    plot_drop(res()$drop)
  }, width  = function() { r <- res(); if (is.null(r)) 820 else max(820, 60 * nrow(r$drop)) },
  height = 460
  )
  
  # ---- Tab "Code": R code that gives the results ----
  output$code_text <- renderText({
    validate(need(!is.null(res()), "Run the simulation first."))
    res()$code
  })
  
  # (The "Methods" tab is static text, defined in the UI.)
  
  # ---- Download button: only available once a simulation has been run ----
  output$dl_ui <- renderUI({
    if (is.null(res())) return(NULL)
    downloadButton("dl", "Download results (RDS)", class = "w-100")
  })
  
  # ---- Everything in one .rds file ----
  output$dl <- downloadHandler(
    filename = function() sprintf("simulation_results_%s.rds", format(Sys.time(), "%Y%m%d_%H%M")),
    content = function(file) {
      r <- res()
      req(r)
      saveRDS(list(
        created = Sys.time(),
        r_version = R.version.string,
        parameters = r$params,
        tau_H1 = r$tau_H1,                   # matrix tau[original, assigned] used under H1
        summary = r$summary,                 # type I error, power (%), Wilson interval around 5%
        dropped_strata = r$drop,             # % of simulations dropping each stratum, under H0 and H1
        mean_attempts = r$attempts,          # IA only
        replications = r$replications,       # raw output of each replication: list(H0 = ..., H1 = ...)
        n_failed = r$n_failed,
        code = r$code,                       # R code reproducing the results
        methods = methods_text(r$params)     # description of the methods for the scenario that was run
      ), file)
    }
  )
}

shinyApp(ui, server)