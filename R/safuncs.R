# This is an R package containing useful functions for my work.

# Available functions with documentation:
# 1. Con_Simul() -- simulates contingency tables based on the multinomial distribution.
# 1b. Con_Simul_PR() -- calculates positive rates for statistical tests on contingency tables.
# 6. Surv_Gen() -- generate rows of survivors given a starting number of fish per tank and data containing morts and sampled fish.
# 7. Surv_Plots() -- generate Kaplan-Meier survival curve and hazard curve from survival data.

# Available functions without documentation:
# 2. Simul_Con_MULT.FISH.ORD() -- simulates ordinal-distributed data across treatments and lesions with inter-fish variation in the PO.
# 3. Simul_Surv() -- simulate survival data based on a reference hazard function, the specified hazard ratio(s), and inter-tank variation.
# 4. theme_Publication() -- ggplot theme for generating publication-ready plots.
# 5. Surv_Pred() -- predict future survival rate(s) for ongoing experiment based on a reference hazard function from older data.

###################################################################################################################################
################################################## Function 1 - Con_Simul() ######################################################

#' @title Simulate Contingency Table
#'
#' @description Simulate a contingency table with fish counts distributed across \emph{n} lesion categories and \emph{n} treatment groups. Probability values for generating counts in each cell (i.e. each factor level combination) can be assigned using the \code{probs} argument. This function is designed for use in power and/or false positive rate calculations; for details, see \code{Con_Simul_PR()}.
#'
#' @details Counts are simulated from a multinomial distribution using \code{rmultinom()}. Counts may be assumed to have a fixed total in the marginals (e.g. per treatment group) or no fixed total in row or column marginals.
#'
#' For further discussion into the types of marginals in contingency tables, refer to: \url{https://www.uvm.edu/~statdhtx/StatPages/More_Stuff/Chi-square/Contingency-Tables.pdf} and the comments on \bold{Arguments}.
#'
#' @param probs Matrix of probability values created using \code{matrix()}. Each row in the matrix should represent a treatment group and each column a lesion category. All probability values in the matrix should sum to 1. Default = equal probability across all cells.
#' @param total_count Total number of counts in the contingency table. Defaults to 750.
#' @param n_lesion Number of lesion categories. Ignored if \code{probs} specified. Defaults to 3.
#' @param n_Trt. Number of treatment groups. Ignored if \code{probs} specified. Defaults to 5.
#' @param margin_fixed_Trt. Whether margins are fixed per treatment group (i.e. fixed number of fish per treatment). Default = FALSE. See \bold{Details} for further information on marginals.
#' @param verbose Whether to print the parameters and probability matrix used. Default = TRUE.
#'
#' @return Returns a list containing:\tabular{lll}{
#'  \code{sim_tab} \tab \tab The simulated contingency table containing counts across different treatment groups (rows) and lesion categories (columns) \cr
#'  \code{params} \tab \tab The simulation parameters as a vector \cr
#'  \code{probs} \tab \tab The probability matrix used for simulation \cr
#' }
#'
#' @export
#'
#' @examples
#' Con_Tab = Con_Simul(total_count = 750, n_lesion = 3, n_Trt. = 5)
#' Con_Tab = Con_Simul(probs = matrix(nrow = 2, ncol = 3, c(1/6, 3/12, 1/6, 1/6, 1/6, 1/12)))
Con_Simul = function(probs = "equal",
                      total_count = 750,
                      n_lesion = 3,
                      n_Trt. = 5,
                      margin_fixed_Trt. = FALSE,
                      verbose = TRUE) { #default conditions

  if (probs[1] == "equal") {
    probs = matrix(nrow = n_Trt., ncol = n_lesion, 1/(n_lesion * n_Trt.))
  }

  if (is.matrix(probs)) {
    if(abs(1 - sum(probs)) > 0.01) {stop("Will not run. Sum of all probabilities must be near equal to 1 (+/- 0.01)")}
    n_Trt. = nrow(probs)
    n_lesion = ncol(probs)
  }


  #Generate vector of counts
  if(margin_fixed_Trt. == FALSE) {
    Count_Mat = rmultinom(n = 1, size = total_count, p = probs)
  } else {
    row_prob = apply(probs, MARGIN = 1, FUN = sum)
    Count_Mat = t(apply(probs, MARGIN = 1, FUN = rmultinom, n = 1, size = 750 * row_prob))
  }

  #Create contingency table
  Con_Tab = matrix(Count_Mat, nrow = n_Trt., ncol = n_lesion,
                   dimnames = list(Trt. = LETTERS[1:n_Trt.],
                                   Lesion_Category = 1:n_lesion))

  if(verbose == FALSE) {
    return(Con_Tab)
  } else {
    return(list(sim_tab = Con_Tab,
                params = c(total_count = total_count, n_lesion = n_lesion, n_Trt. = n_Trt., margin_fixed_Trt. = margin_fixed_Trt.),
                probs = probs))
  }
}

###################################################################################################################################
################################################## Function 1b - Con_Simul_PR() #################################################

#' @title Calculate Positive Rates for Contingency Table
#'
#' @description Computes statistical power and optionally false positive rates for tests applied to contingency tables based on Monte Carlo simulations. Specify the simulation process using \code{Con_Simul()}, which serves as input. Positive rates are computed for the Chi-square test and optionally for Fisher's exact test and the Wald test applied to an ordinal regression model.
#'
#' @details Power is defined as the percentage of tests yielding positive results (p-value < 0.05) on a set of contingency tables simulated based on the data-generating process and probability matrix specified in \code{Con_Simul()}. The specified probability matrix should represent the parameters of the population where there is a desire to detect a significant effect in the sample. The simulated contingency tables then reflect the sample outcomes of the specified population parameter.
#'
#' False Positive Rate is defined as the percentage of tests yielding positive results (p-value < 0.05) on contingency tables simulated based on a probability matrix without treatment effect. The odd ratios between lesion scores for all treatment groups are assumed to be the same as that for control (first row in the probability matrix).
#'
#' P-values for the Chi-square test are computed using \code{stats::chisq.test()} with default parameters. P-values for Fisher's exact test are computed using \code{stats::fisher.test()} with \code{simulate.p.value} set to TRUE, alongside default parameters. P-values for ordinal regression are computed from \code{stats::anova()} applied to the output of \code{ordinal::clm()}.
#'
#' @param Con_Simul_Object Output from \code{Con_Simul()} with the argument \code{verbose} set to TRUE.
#' @param add_fisher_exact Whether to compute positive rates for Fisher's Exact test. May add >1 min of calculation time. Defaults to FALSE.
#' @param add_ord Whether to compute positive rates for Wald test on a fitted ordinal regression model. May add >1 min of calculation time. Defaults to FALSE.
#' @param sample_sizes A vector of sample sizes over which false positive rates are to be calculated. A sample size is defined as the total number of counts in a contingency table. Defaults to total count received by \code{Con_Simul_Object}.
#' @param n_sim Number of contingency tables simulated for each positive rate calculation. Defaults to 1000.
#' @param FPR Whether to calculate false positive rate in addition to power. Defaults to TRUE.
#'
#' @return Returns a list containing:\tabular{lll}{
#'  \code{pos_rate} \tab \tab A dataframe containing positive rates of various tests at specified sample sizes \cr
#'  \code{eff_mat} \tab \tab The probability matrix used to calculate power  \cr
#'  \code{null_mat} \tab \tab The probability matrix used to calculate false positive rate \cr
#'  \code{plot} \tab \tab A plot of positive rates of various tests at specified sample sizes \cr
#' }
#'
#' @export
#'
#' @examples
#' ## Below I show how we can perform a simple power calculation using this tool.
#' ## Suppose I want to calculate power for Treatment B which halves the lesions in category 2 and 3.
#' ## I then specify the following probability matrix and feed it into Con_Simul():
#' probs_mat = matrix(nrow = 2, ncol = 3, data = c(1/6, 1/3, 1/6, 1/12, 1/6, 1/12))
#' sim_tab = Con_Simul(probs_mat)
#'
#' ## Next, I feed the output into Con_Simul_PR():
#' Con_Simul_PR(sim_tab, sample_sizes = c(50, 100, 150))
#' ## Results: Power is ~55, 86, and 97% for the Chi-square test using total counts of 50, 100, and 150, respectively.
#'
#' ## The same power for Chi-square test can be calculated using Cohen's omega (w) method which is faster but has its own limitations;
#' ## e.g. assumes one data generating process for the contingency table (the no fixed marginals).
#' library(pwr)
#' pwr.chisq.test(w = ES.w2(probs_mat), df = 2, sig.level = 0.05, N = 100)
#' ## Results: Power is 85.6% for the Chi-square test at the total count of 100.
Con_Simul_PR = function(Con_Simul_Object = Con_Simul(),
                        add_fisher_exact = FALSE,
                        add_ord = FALSE,
                        sample_sizes = NULL,
                        n_sim = 1000,
                        FPR = TRUE) {

  Params = Con_Simul_Object$params
  PR_DB = data.frame()
  if(is.null(sample_sizes)) {sample_sizes <- Params[1]}

  for(tot_count in as.vector(sample_sizes)) {

    P_Chisq1_Eff = c()
    P_Chisq1_Null = c()
    P_Ord_Eff = c()
    P_Ord_Null = c()
    P_Fish_Eff = c()
    P_Fish_Null = c()

    for(iter in 1:n_sim) {
      Sim_Tab_Eff = Con_Simul(total_count = tot_count,
                               n_lesion = Params[2],
                               n_Trt. = Params[3],
                               margin_fixed_Trt. = Params[4],
                               probs = Con_Simul_Object$probs,
                               verbose = FALSE)

      P_Chisq1_Eff = append(P_Chisq1_Eff, suppressWarnings(chisq.test(x = Sim_Tab_Eff)$p.value))
      if(add_ord == TRUE) {
        DB_Ord_Eff = as.data.frame(as.table(Sim_Tab_Eff))
        DB_Ord_Eff = data.frame(lapply(DB_Ord_Eff, rep, DB_Ord_Eff$Freq))
        Ord_Eff = ordinal::clm(data = DB_Ord_Eff, factor(Lesion_Category, ordered = TRUE) ~ Trt.)
        P_Ord_Eff = append(P_Ord_Eff, anova(Ord_Eff)$'Pr(>Chisq)')
      }
      if(add_fisher_exact == TRUE) {
        P_Fish_Eff <- append(P_Fish_Eff, fisher.test(x = Sim_Tab_Eff, simulate.p.value = TRUE)$p.value)
      }
    }

    if(FPR == TRUE) {
      probs_null_mat = matrix(nrow = nrow(Con_Simul_Object$probs),
                              ncol = ncol(Con_Simul_Object$probs),
                              data = rep(Con_Simul_Object$probs[1,],
                                         each = nrow(Con_Simul_Object$probs)))
      probs_null_mat = probs_null_mat * rowSums(Con_Simul_Object$probs)
      probs_null_mat = probs_null_mat / sum(probs_null_mat)

      for(iter in 1:n_sim) {

        Sim_Tab_Null = Con_Simul(total_count = tot_count,
                                  n_lesion = Params[2],
                                  n_Trt. = Params[3],
                                  margin_fixed_Trt. = Params[4],
                                  probs = probs_null_mat,
                                  verbose = FALSE)


        P_Chisq1_Null = append(P_Chisq1_Null, suppressWarnings(chisq.test(x = Sim_Tab_Null)$p.value))
        if(add_ord == TRUE) {
          DB_Ord_Null = as.data.frame(as.table(Sim_Tab_Null))
          DB_Ord_Null = data.frame(lapply(DB_Ord_Null, rep, DB_Ord_Null$Freq))
          Ord_Null = ordinal::clm(data = DB_Ord_Null, factor(Lesion_Category, ordered = TRUE) ~ Trt.)
          P_Ord_Null = append(P_Ord_Null, anova(Ord_Null)$'Pr(>Chisq)')
        }
        if(add_fisher_exact == TRUE) {
          P_Fish_Null <- append(P_Fish_Null, fisher.test(x = Sim_Tab_Null, simulate.p.value = TRUE)$p.value)
        }
      }
    }

    PR_DB = rbind(PR_DB, data.frame(Total_count = tot_count,
                                    Power_Chisquare = sum(P_Chisq1_Eff < 0.05) / length(P_Chisq1_Eff),
                                    Power_OrdinalRegression = sum(P_Ord_Eff < 0.05) / length(P_Ord_Eff),
                                    Power_FishersExact = sum(P_Fish_Eff < 0.05) / length(P_Fish_Eff),
                                    FPR_Chisquare = sum(P_Chisq1_Null < 0.05) / length(P_Chisq1_Null),
                                    FPR_OrdinalRegression = sum(P_Ord_Null < 0.05) / length(P_Ord_Null),
                                    FPR_FishersExact = sum(P_Fish_Null < 0.05) / length(P_Fish_Null)))
  }

  PR_DB_stacked = data.frame(PR_DB$Total_count,
                             stack(x = PR_DB, select = c("Power_Chisquare", "Power_OrdinalRegression", "Power_FishersExact",
                                                         "FPR_Chisquare", "FPR_OrdinalRegression", "FPR_FishersExact")))
  PR_DB_stacked = na.omit(PR_DB_stacked)
  split_txt = strsplit(x = as.character(PR_DB_stacked$ind), split = "_")
  PR_DB_stacked$Statistical_Test = sapply(split_txt, '[', 2)
  PR_DB_stacked$Class = sapply(split_txt, '[', 1)
  PR_DB_stacked = PR_DB_stacked[,-3]
  colnames(PR_DB_stacked) <- c("Total_Count", "Percent_of_Significant_Results", "Statistical_Test", "Class")
  PR_DB_stacked$Percent_of_Significant_Results = 100 * PR_DB_stacked$Percent_of_Significant_Results

  #remove na columns
  PR_DB = PR_DB[, colSums(is.na(PR_DB)) < nrow(PR_DB)]

  plot1 = ggplot2::ggplot(data = PR_DB_stacked, ggplot2::aes(x = Total_Count, y = Percent_of_Significant_Results, color = Statistical_Test)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ Class) +
    ggplot2::xlab("Total Counts in Contingency Table") +
    ggplot2::ylab("Percent of Significant Results (%)") +
    ggplot2::scale_y_continuous(breaks = seq(0, 100, 5), limits = c(0, 100)) +
    ggplot2::labs(color = "Statistical Test")

  return(list(pos_rate = PR_DB_stacked,
              eff_mat = Con_Simul_Object$probs,
              null_mat = probs_null_mat,
              plot = plot1))
}


###################################################################################################################################
############################################# Function 2 - Simul_Con_MULT.FISH.ORD() ##############################################

Simul_Con_MULT.FISH.ORD = function(total_count = 15000,
                                   n_Trt. = 5,
                                   n_Fish_total = 750,
                                   base_OR_vec2 = c(1/2, 2/1),
                                   fish_LOR_var = 0,
                                   fish_susceptibility_var = 0) { #default conditions

  n_lesion = 3 #only number of lesion scores supported due to difficulty in simulating log odds variation of greater complexity
  odds_vec = c()

  for(fish_num in 1:n_Fish_total) {

    #Generate variable log-odds between lesion scores
    base_LOR_vec = log(base_OR_vec2) #LOR = Log_odds_ratio.
    LOR_var = rnorm(n = 1, mean = 0, sd = fish_LOR_var)
    new_LOR_vec = base_LOR_vec + LOR_var

    #Calculate odds for each lesion score. Got these by solving simultaneous equations (work not shown)
    odds_LS3 = (1 + (1 / exp(new_LOR_vec[1]))) / (exp(new_LOR_vec[2]) + 1)
    odds_LS2 = (exp(new_LOR_vec[2]) * odds_LS3) - 1
    odds_LS1 = 1

    odds_vec = append(odds_vec, c(odds_LS1, odds_LS2, odds_LS3) / sum(c(odds_LS1, odds_LS2, odds_LS3)))
  }

  #Generate vector of counts
  fish_susc_var = rep(rnorm(mean = log(1/(n_lesion * n_Fish_total)),
                            sd = fish_susceptibility_var, n = n_Fish_total), each = n_lesion) #variation in fish susceptibility in the log-scale

  fish_ls_var = exp(fish_susc_var) * odds_vec #variation in p between lesions, returned to absolute scale
  Count_Vec1 = rmultinom(n = 1, size = total_count, p = fish_ls_var)

  #Create Count Dataframe
  Count_DF = data.frame(Trt. = rep(LETTERS[1:n_Trt.], each = n_Fish_total * n_lesion / n_Trt.),
                        LS = rep(1:n_lesion, times = n_Fish_total),
                        Fish.ID = rep(1:n_Fish_total, each = n_lesion),
                        Fish_Susc_Var = exp(fish_susc_var),
                        Fish_LS_Var = fish_ls_var,
                        Counts = Count_Vec1)

  return(Count_DF)
}

#' @export

###################################################################################################################################
################################################### Function 3 - Simul_SURV() #####################################################

Simul_SURV = function(Haz_DB = Simul_SURV_Haz_DB, #object from bshazard() as reference hazard function. See data(Simul_SURV_Haz_DB)
                      fish_num_per_tank = 100,
                      tank_num_per_trt = 4,
                      Treatments_HR = c(1, 1, 1, 1), #HR for treatment groups, starting from control (1)
                      logHR_sd_intertank = 0.35, #inter-tank variation (sd) in HR in the log-scale
                      round_digits_days = 0) { #Time is rounded to Days (i.e. fish cannot die at Day 1.4 or any other decimal)
  require(devtools)
  CDF_Yval = c()

  for(Treatment_Term in Treatments_HR) {

    for(Tank_num in 1:tank_num_per_trt) {

      Tank_eff = rnorm(n = 1, mean = 0, sd = logHR_sd_intertank)

      U = runif(n = fish_num_per_tank, min = 0, max = 1)

      CDF_Yval_temp = -log(U) * exp(-(log(Treatment_Term) + Tank_eff))
      CDF_Yval = append(CDF_Yval, CDF_Yval_temp)
    }
  }

  #Get Time to Event
  TTE = approx(x = cumsum(Haz_DB$hazard), y = Haz_DB$time, xout = CDF_Yval, method = "linear")$y
  TTE = round(TTE, digits = round_digits_days)

  #Turn NA (from out of bound CDF_Yval) to the last follow up time
  TTE = ifelse(is.na(TTE), max(Haz_DB$time), TTE)

  #Label Status (1 - dead, or 0 - survived) given TTE
  Surv_simul_DB = data.frame(TimetoEvent = TTE,
                             CDF_Yval = CDF_Yval,
                             Status = ifelse(TTE == max(Haz_DB$time), 0, 1),
                             Treatment_G = as.factor(rep(LETTERS[1:length(Treatments_HR)], each = fish_num_per_tank * tank_num_per_trt)),
                             Tank_num = as.factor(rep(1:(length(Treatments_HR) * tank_num_per_trt), each = fish_num_per_tank)))

  return(Surv_simul_DB)
}

###################################################################################################################################
############################################### Function 4 - theme_Publication() ##################################################

theme_Publication = function(base_size = 14, base_family = "helvetica") {
  require(grid)
  require(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.6),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_blank(),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.4, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(size=12),
            plot.margin=unit(c(10,15,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

###################################################################################################################################
################################################## Function 5 - Surv_Pred() #######################################################

#' Predict Survival Rate
#' @description Predict survival rate for a given survival dataset provided a reference survival database used to estimate a reference hazard curve. Prediction done seperately by treatment group.
#'
#' @details P
#'
#' @param pred_db Placeholder
#' @param ref_db Placeholder
#' @param predsr_tte Placeholder
#' @param method Placeholder
#' @param coxph_mod Placeholder
#' @param lambda_pred Placeholder
#'
#' @return Placeholder
#' @export
#'
#' @examples Placeholder
Surv_Pred = function(pred_db, #Data from ongoing study, with SR to be predicted. See \bold{Details} for specifics.
                     ref_db, #Reference survival data from to create the reference hazard function.
                     predsr_tte, #The day at which SR is to be predicted. Minimum is Day 5 post challenge.
                     method = 2, #SR prediction method. See \bold{Details} for more info.
                     coxph_mod = "GLMM", #Model used to estimate HR. Can be either "GLMM" or "GEE". See \bold{Details} for more info.
                     lambda_pred = NULL, #Lambda parameter for the bshazard curve of the predicted dataset.
                     phi_pred = NULL)
{

  pred_db = pred_db[pred_db$TTE > 0, ] #ensure positive TTE
  ref_db = ref_db[ref_db$TTE > 0, ] #ensure positive TTE

  #Generate reference level bshazard curve
  ref_id = levels(as.factor(ref_db$Trt.ID))
  ref_db$Trt.ID = paste("ref", ref_id)
  ref_bshaz = bshazard::bshazard(data = ref_db, survival::Surv(TTE, Status) ~ Tank.ID,
                                 verbose = FALSE)

  #initialize dataframes
  pred_SR_DB = data.frame()
  haz_db = data.frame()

  #loop for every treatment
  for(pred_trt in levels(as.factor(pred_db$Trt.ID))) {

    comb_db = rbind(ref_db, pred_db[pred_db$Trt.ID == pred_trt,])
    comb_db$Trt.ID = relevel(as.factor(comb_db$Trt.ID), ref = paste("ref", ref_id))

    #loop for every predictable day
    for(SR_Day in 5:(predsr_tte-1)) {
      comb_db2 = survival::survSplit(comb_db, cut = SR_Day, end = "TTE", event = "Status", episode = "Obs")
      comb_db2 = comb_db2[comb_db2$Obs == 1, ]

      #Get HR
      if(coxph_mod == "GLMM"){cox_comp = coxme::coxme(survival::Surv(TTE, Status) ~ Trt.ID + (1|Tank.ID),
                                                      data = comb_db2)}
      if(coxph_mod == "GEE"){cox_comp = survival::coxph(survival::Surv(TTE, Status) ~ Trt.ID, cluster = Tank.ID,
                                                        data = comb_db2)}
      pred_HR = exp(coef(cox_comp))

      #Get SR
      ref_bshaz_t = data.frame(hazard = ref_bshaz$hazard[ref_bshaz$time < predsr_tte],
                               time = ref_bshaz$time[ref_bshaz$time < predsr_tte])

      if(method == 1) {
        cumhaz = DescTools::AUC(x = c(ref_bshaz_t$time),
                                y = c(ref_bshaz_t$hazard) * pred_HR)
        pred_SR = 100 * exp(-cumhaz)
      }

      if(method == 2) {
        pred_db2 = survival::survSplit(pred_db, cut = SR_Day, end = "TTE", event = "Status", episode = "Obs")
        pred_db2 = pred_db2[pred_db2$Obs == 1, ]
        pred_bshaz = bshazard::bshazard(data = droplevels(pred_db2[pred_db2$Trt.ID == pred_trt,]),
                                        survival::Surv(TTE, Status) ~ Tank.ID, verbose = FALSE, lambda = lambda_pred)

        cumhaz_precut = DescTools::AUC(x = c(pred_bshaz$time),
                                       y = c(pred_bshaz$hazard))
        ref_bshaz_t2 = ref_bshaz_t[ref_bshaz_t$time > SR_Day,]
        cumhaz_postcut = DescTools::AUC(x = c(ref_bshaz_t2$time),
                                        y = c(ref_bshaz_t2$hazard) * pred_HR)

        if(is.na(cumhaz_postcut)) {cumhaz_postcut <- 0}
        pred_SR = 100 * exp(-(cumhaz_precut + cumhaz_postcut))
      }
      #compile survival prediction database
      pred_SR_DB = rbind(pred_SR_DB, data.frame(Trt.ID = pred_trt, Observable_SR_Day = SR_Day, pred_SR, pred_HR))
    }
  }
  row.names(pred_SR_DB) = NULL

  #create survival prediction plot
  library(ggplot2)
  Pred_SR_Plot = ggplot(data = pred_SR_DB, aes(x = Observable_SR_Day, y = pred_SR, color = Trt.ID)) +
    geom_point() +
    geom_line() +
    facet_wrap(~Trt.ID) +
    scale_y_continuous(name = paste("Predicted Survival (%) on TTE = ", predsr_tte), breaks = seq(0, 100, 10), limits = c(0, 100)) +
    scale_x_continuous(name = "Observable TTEs used in Prediction", breaks = seq(0, 100, 1))

  Pred_HR_Plot = ggplot(data = pred_SR_DB, aes(x = Observable_SR_Day, y = pred_HR, color = Trt.ID)) +
    geom_point() +
    geom_line() +
    facet_wrap(~Trt.ID) +
    scale_y_continuous(name = paste("Predicted HR ", predsr_tte)) +
    scale_x_continuous(name = "Observable TTEs used in Prediction", breaks = seq(0, 100, 1))

  return(list(pred_SR_DB, Pred_SR_Plot, Pred_HR_Plot))
}

###################################################################################################################################
################################################## Function 6 - Surv_Gen() ########################################################

#' @title Generate Survivor Data
#'
#' @description Returns a survival dataframe that includes rows representing every surviving fish based on the starting number of fish provided and the mortality count provided in a dataframe. To generate survivor data for tanks absent in the input dataframe, specify the arguments \code{tank_without_mort} and \code{trt_without_mort}.
#'
#' @details The mort dataframe supplied as input should consist of the following 4 columns at minimum:
#' * "Trt.ID" = Labels for treatment groups in the study.
#' * "Tank.ID" = Labels for tanks in the study (each tank must have a unique label).
#' * "TTE" = Time to Event. Event could be fish death or being sampled and removed depending on "Status".
#' * "Status" = Value indicating what happened at TTE. 1 for dead fish, 0 for those sampled and removed.
#'
#' Each row should represent one fish.
#'
#' For an example dataset, view \code{data(mort_db_ex)}.
#' @md
#'
#' @param mort_db A mort dataframe as described in \bold{Details}.
#' @param starting_fish_count Value representing the starting number of fish per tank.
#' @param today_tte Value representing the day or time-to-event the fish survived to. Value assigned to every row of survivor data generated.
#' @param tank_without_mort A vector of strings specifying the tanks absent from \code{mort_db} for which to generate survivor data.
#' @param trt_without_mort A vector of strings corresponding to \code{tank_without_mort}. Keep the order the same.
#'
#' @return A dataframe produced by combining the input mort data and generated rows of survivor data.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' Surv_Gen(mort_db = mort_db_ex,
#'          today_tte = 54,
#'          starting_fish_count = 100,
#'          tank_without_mort = c("C99", "C100"),
#'          trt_without_mort = c("A", "B"))
#'
Surv_Gen = function(mort_db,
                    starting_fish_count,
                    today_tte,
                    tank_without_mort = NULL,
                    trt_without_mort = NULL) {

  DB_Mort_Gensum = data.frame(mort_db %>%
                                dplyr::group_by(Trt.ID, Tank.ID) %>%
                                dplyr::summarise(Num_dead = n()))

  if(!is.null(tank_without_mort) && !is.null(trt_without_mort)) {
    WM_DB = data.frame(Trt.ID = trt_without_mort,
                       Tank.ID = tank_without_mort,
                       Num_dead = 0)
    DB_Mort_Gensum = rbind(DB_Mort_Gensum, WM_DB)
  }

  if(is.data.frame(starting_fish_count)) {
    DB_Mort_Gensum = merge(DB_Mort_Gensum, starting_fish_count)
  } else {DB_Mort_Gensum$starting_fish_count = starting_fish_count}

  DB_Mort_Gensum$Num_alive = DB_Mort_Gensum$starting_fish_count - DB_Mort_Gensum$Num_dead
  DB_Mort_Genalive = data.frame(lapply(DB_Mort_Gensum, rep, DB_Mort_Gensum$Num_alive))
  DB_Mort_Genalive$Status = 0
  DB_Mort_Genalive$TTE = today_tte
  DB_Mort_Gencomb = plyr::rbind.fill(mort_db, DB_Mort_Genalive[, -c(3:5)])

  return(DB_Mort_Gencomb)
}

###################################################################################################################################
################################################# Function 7 - Surv_Plots() #######################################################

#' Generate Survival Plots
#'
#' @description Uses survival data to produce a Kaplan-Meier Survival Plot and/or Hazard Time Plot. Each plot contains multiple curves for the different treatment groups. Plots saved automatically to working directory.
#'
#' @details The survival dataset should be a dataframe containing at least 4 different columns:
#' * "Trt.ID" = Labels for treatment groups in the study.
#' * "Tank.ID" = Labels for tanks in the study (each tank must have a unique label).
#' * "TTE" = Time to Event. Event depends on "Status".
#' * "Status" = Value indicating what happened at TTE. 1 for dead fish, 0 for survivors or those sampled and removed.
#'
#' Each row should represent one fish. For an example dataset, view \code{data(surv_db_ex)}.
#'
#' For details on the statistical methodology used by \code{bshazard()}, refer to: \url{https://www.researchgate.net/publication/287338889_bshazard_A_Flexible_Tool_for_Nonparametric_Smoothing_of_the_Hazard_Function}.
#'
#' General concept: the author first considers h(t) the hazard function in an event-count model: count(t) = h(t) * P(t) where P(t) is the number alive as a function of time. h(t) is modeled over time using basis splines. By assuming the basis spline curvatures have a normal distribution with mean 0 (a random effect), the author found that the curvature's variance can be estimated as a function of the degree of over-dispersion (phi) of counts. The author determined that the variance of curvatures (smoothness) is equal to phi divided by a denominator lambda representing the smoothness parameter. The user can specify phi or lambda which affects the variance of curvatures (smoothness) or can allow them to be estimated from the data.
#' @md
#'
#' @param surv_db A survival dataframe as described in \bold{Details}.
#' @param plot_prefix A string specifying the prefix for the filename of the saved plots.
#' @param xlim A vector specifying the plots x-axis lower and upper limits, respectively.
#' @param ylim A vector specifying the Survival Plot y-axis lower and upper limits, respectively.
#' @param xlab A string specifying the plot x-axis label.
#' @param lambda Smoothing value for the hazard curve. Higher lambda produces greater smoothing. Defaults to NULL where \code{bshazard()} uses the provided survival data to estimate lambda; recommended for large sample size situations which usually occurs on our full-scale studies with many mortalities and tank-replication. At low sample sizes, the lambda estimate can be unreliable. Choosing a lambda of 10 (or anywhere between 1-100) probably produces the most accurate hazard curve for these situations. In place of choosing lambda, choosing \code{phi} is recommended; see below.
#' @param phi Dispersion parameter for the count model used in hazard curve estimation. Defaults to NULL where \code{bshazard()} uses the provided survival data to estimate phi; recommended for large sample size situations. At low sample sizes, the phi estimate can be unreliable. Choosing a phi value of 1 for low sample sizes is recommended. This value of 1 (or close) seems to be that estimated in past Tenaci data (phi ~ 0.8-1.4) in Onda where there are large sample sizes with tank-replication. The value of 1 indicates counts (deaths) vary as "expected" (Poisson) according to the different hazard values along the hazard curve, and are not overdispersed (phi > 1).
#' @param dailybin Whether to set time bins at daily (1 TTE) intervals. Refer to the \code{bshazard()} documentation for an understanding on the role of bins to hazard curve estimation. Please set to TRUE at low sample sizes and set to FALSE at large sample sizes with tank-replication. Defaults to TRUE.
#' @param plot Which plot to output. Use "surv" for the Kaplan-Meier Survival Curve, "haz" for the Hazard Curve, or "both" for both. Defaults to "both".
#' @param colours Vector of color codes for the different treatment groups in the plot. Defaults to ggplot2 default palette.
#'
#' @return If \code{plot == "surv"}, returns a ggplot2 object reflecting the Kaplan-Meier Survival Curve.
#' If \code{plot == "haz"}, returns a ggplot2 object reflecting the Hazard Curve.
#' If \code{plot == "both"}, returns both ggplot2 objects in a list.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' Surv_Plots(surv_db = surv_db_ex,
#'            plot_prefix = "QCATC777",
#'            xlim = c(0, 50),
#'            ylim = c(0, 1),
#'            xlab = "TTE",
#'            phi = 1,
#'            plot = "both")
#'
Surv_Plots = function(surv_db,
                      plot_prefix = "plot_prefix",
                      xlim = NULL,
                      ylim = c(0, 1),
                      xlab = "Days Post Challenge",
                      lambda = NULL,
                      phi = NULL,
                      dailybin = TRUE,
                      plot = "both",
                      colours = NULL) {

  if(is.null(xlim)) {xlim <- c(0, max(surv_db$TTE))}

  if(plot == "surv" | plot == "both") {
  surv_obj = survminer::surv_fit(survival::Surv(TTE, Status) ~ Trt.ID, data = surv_db)
  attributes(surv_obj$strata)$names <- levels(as.factor(surv_db$Trt.ID))

  surv_plot = survminer::ggsurvplot(surv_obj,
                                    conf.int = FALSE,
                                    ggtheme = theme(plot.background = element_rect(fill = "white")),
                                    break.y.by = 0.1,
                                    break.x.by = max(round(max(xlim) / 15), 1),
                                    xlim = xlim,
                                    ylim = ylim,
                                    xlab = xlab,
                                    surv.scale = "percent")
  Survival_Plot = surv_plot$plot + theme(legend.position = "right") + guides(color = guide_legend("Trt.ID"))

  if(!is.null(colours)) {Survival_Plot = Survival_Plot + scale_color_manual(values = colours)}
  ggsave(paste(plot_prefix, "Survival Curve.tiff"), dpi = 300, width = 6, height = 4, plot = Survival_Plot)

  }

  if(dailybin == TRUE) {dbin <- max(surv_db$TTE)}
  if(dailybin == FALSE) {dbin <- NULL}

  #Create Haz_list

  if(plot == "haz" | plot == "both") {
  Haz_list = list()
  for(Haz_Trt in levels(as.factor(surv_db$Trt.ID))) {
    surv_db_trt = surv_db[surv_db$Trt.ID == Haz_Trt,]
    if(length(levels(as.factor(surv_db_trt$Tank.ID))) > 1) {
      Haz_bs = bshazard::bshazard(nbin = dbin,
                                  data = surv_db_trt,
                                  survival::Surv(TTE, Status) ~ Tank.ID,
                                  verbose = FALSE,
                                  lambda = lambda,
                                  phi = phi)
    } else {
      Haz_bs = bshazard::bshazard(nbin = dbin,
                                  data = surv_db_trt,
                                  survival::Surv(TTE, Status) ~ 1,
                                  verbose = FALSE,
                                  lambda = lambda,
                                  phi = phi)
    }
    Haz_list[[Haz_Trt]] = data.frame(Hazard = Haz_bs$hazard,
                                     Time = Haz_bs$time)
  }

  haz_db = dplyr::bind_rows(Haz_list, .id = "Trt.ID")
  Hazard_Plot = ggplot(data = haz_db, aes(x = Time, y = Hazard, color = Trt.ID)) +
    geom_line(linewidth = 1) +
    geom_point() +
    xlab(xlab) +
    scale_x_continuous(breaks = seq(from = min(xlim),
                                    to = max(xlim),
                                    by = max(round(max(xlim) / 15), 1)),
                       limits = xlim)

  if(!is.null(colours)) {Hazard_Plot = Hazard_Plot + scale_color_manual(values = colours)}
  ggsave(paste(plot_prefix, "Hazard Curve.tiff"), dpi = 300, width = 6, height = 4, plot = Hazard_Plot)

  }

  if(plot == "surv") {return(Survival_Plot)}
  if(plot == "haz") {return(Hazard_Plot)}
  if(plot == "both") {return(list(Survival_Plot, Hazard_Plot))}
}
