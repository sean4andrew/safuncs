# This is an R package containing useful functions for my work.

# Available functions with documentation:
# 1. Simul_Mult() -- simulates contingency tables based on the multinomial distribution.
# 1b. Pow_Simul_Mult() -- calculates positive rates for statistical tests on contingency tables.

# Available functions without documentation:
# 2. Simul_Con_MULT.FISH.ORD() -- simulates ordinal-distributed data across treatments and lesions with inter-fish variation in the PO. WITHOUT HELP PAGE/DOCUMENTATION.
# 3. Simul_Surv() -- simulate survival data based on a reference hazard function, the specified hazard ratio(s), and inter-tank variation.
# 4. theme_Publication() -- ggplot theme for generating publication-ready plots.
# 5. Predict_SR() -- predict future survival rate(s) for ongoing experiment based on a reference hazard function from older data.
# 6. Surv_Gen0() -- generate rows of survivors given a starting number of fish per tank and data containing morts and sampled fish.

###################################################################################################################################
################################################## Function 1 - Simul_Mult() ######################################################

#' @title Simulate a Contingency Table
#'
#' @description Simulate a contingency table with fish counts distributed across \emph{n} lesion categories and \emph{n} treatment groups. Probability values for generating counts in each cell (i.e. each factor level combination) can be assigned using the \code{probs} argument. This function is designed for use in power and/or false positive rate calculations; for details, see \code{Pow_Simul_Mult()}.
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
#' Con_Tab = Simul_Mult(total_count = 750, n_lesion = 3, n_Trt. = 5)
#' Con_Tab = Simul_Mult(probs = matrix(nrow = 2, ncol = 3, c(1/6, 3/12, 1/6, 1/6, 1/6, 1/12)))
Simul_Mult = function(probs = "equal",
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
################################################## Function 1b - Pow_Simul_Mult() #################################################

#' @title Positive Rates for Contingency Tables
#'
#' @description Computes statistical power and optionally false positive rates for tests applied to contingency tables based on Monte Carlo simulations. Specify the simulation process using \code{Simul_Mult()}, which serves as input. Positive rates are computed for the Chi-square test and optionally for Fisher's exact test and the Wald test applied to an ordinal regression model.
#'
#' @details Power is defined as the percentage of tests yielding positive results (p-value < 0.05) on a set of contingency tables simulated based on the data-generating process and probability matrix specified in \code{Simul_Mult()}. The specified probability matrix should represent the parameters of the population where there is a desire to detect a significant effect in the sample. The simulated contingency tables then reflect the sample outcomes of the specified population parameter.
#'
#' False Positive Rate is defined as the percentage of tests yielding positive results (p-value < 0.05) on contingency tables simulated based on a probability matrix without treatment effect. The lesion odds ratios for the control group (assumed to be the first row in the probability matrix) is assumed to be true across all treatments.
#'
#' P-values for the Chi-square test are computed using \code{stats::chisq.test()} with default parameters. P-values for Fisher's exact test are computed using \code{stats::fisher.test()} with \code{simulate.p.value} set to TRUE, alongside default parameters. P-values for ordinal regression are computed from \code{stats::anova()} applied to the output of \code{ordinal::clm()}.
#'
#' @param Simul_Mult_Object Output from \code{Simul_Mult()} with the argument \code{verbose} set to TRUE.
#' @param add_fisher_exact Whether to compute positive rates for Fisher's Exact test. May add >1 min of calculation time. Defaults to FALSE.
#' @param add_ord Whether to compute positive rates for Wald test on a fitted ordinal regression model. May add >1 min of calculation time. Defaults to FALSE.
#' @param sample_sizes A vector of sample sizes over which false positive rates are to be calculated. A sample size is defined as the total number of counts in a contingency table. Defaults to total count received by \code{Simul_Mult_Object}.
#' @param n_sim Number of contingency tables simulated for each positive rate calculation. Defaults to 1000.
#' @param FPR Whether to calculate false positive rate in addition to power. Defaults to TRUE.
#'
#' @return Returns a list containing:\tabular{lll}{
#'  \code{pos_rate} \tab \tab A data frame containing positive rates of various tests at specified sample sizes \cr
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
#' ## I then specify the following probability matrix and feed it into Simul_Mult():
#' probs_mat = matrix(nrow = 2, ncol = 3, data = c(1/6, 1/3, 1/6, 1/12, 1/6, 1/12))
#' sim_tab = Simul_Mult(probs_mat)
#'
#' ## Next, I feed the output into Pow_Simul_Mult():
#' Pow_Simul_Mult(sim_tab, sample_sizes = c(50, 100, 150))
#' ## Results: Power is ~55, 86, and 97% for the Chi-square test using total counts of 50, 100, and 150, respectively.
#'
#' ## The same power for Chi-square test can be calculated using Cohen's omega (w) method which is faster but has its own limitations:
#' library(pwr)
#' pwr.chisq.test(w = ES.w2(probs_mat), df = 2, sig.level = 0.05, N = 100)
#' ## Results: Power is 85.6% for the Chi-square test at the total count of 100.
Pow_Simul_Mult = function(Simul_Mult_Object = Simul_Mult(),
                          add_fisher_exact = FALSE,
                          add_ord = FALSE,
                          sample_sizes = NULL,
                          n_sim = 1000,
                          FPR = TRUE) {

  Params = Simul_Mult_Object$params
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
      Sim_Tab_Eff = Simul_Mult(total_count = tot_count,
                               n_lesion = Params[2],
                               n_Trt. = Params[3],
                               margin_fixed_Trt. = Params[4],
                               probs = Simul_Mult_Object$probs,
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
      probs_null_mat = matrix(nrow = nrow(Simul_Mult_Object$probs),
                              ncol = ncol(Simul_Mult_Object$probs),
                              data = rep(Simul_Mult_Object$probs[1,],
                                         each = nrow(Simul_Mult_Object$probs)))
      probs_null_mat = probs_null_mat * rowSums(Simul_Mult_Object$probs)
      probs_null_mat = probs_null_mat / sum(probs_null_mat)

      for(iter in 1:n_sim) {

        Sim_Tab_Null = Simul_Mult(total_count = tot_count,
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
              eff_mat = Simul_Mult_Object$probs,
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

  #Create Count Data Frame
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
################################################## Function 5 - Predict_SR() ######################################################

Predict_SR = function(New_DB = Predict_SR_New_DB, #Data from ongoing study, with SR to be predicted.
                      Ref_DB = Predict_SR_Ref_DB, #Reference survival data to create the reference hazard function.
                      End_Day = 18,  #The end date at which SR is to be predicted
                      Method = 2,      #SR prediction method, minor differences between Method 1-2 (choose any should be OK)
                      PH_Mod = "GLMM") #Model used to estimate HR. Can be either "GLMM" or "GEE". GLMM recommended.
  {

  require(bshazard)
  require(coxme)
  require(DescTools)
  require(survival)
  require(dplyr)
  require(SimDesign)
  require(devtools)

  SR_Days = 5:End_Day
  New_DB = New_DB[New_DB$Std_Time > 0, ]

  Ref_ID = levels(as.factor(Ref_DB$Trt.ID))
  Ref_bshaz = quiet(bshazard(data = Ref_DB, Surv(Std_Time, Status) ~ Tank.ID,
                             verbose = FALSE, nbin = max(Ref_DB$Std_Time)))

  Comb_DB = rbind(Ref_DB, New_DB)
  Comb_DB$Trt.ID = relevel(as.factor(Comb_DB$Trt.ID), ref = Ref_ID)

  pred_SR_DB = data.frame()

  for(Comp_HR in levels(as.factor(New_DB$Trt.ID))) {

    for(SR_Day in SR_Days) {
      Comb_DB2 = survSplit(Comb_DB, cut = SR_Day, end = "Std_Time", event = "Status", episode = "Obs")
      Comb_DB2 = Comb_DB2[Comb_DB2$Obs == 1, ]

      #Get HR
      if(PH_Mod == "GLMM"){
        cox_comp = coxme(Surv(Std_Time, Status) ~ Trt.ID + (1|Tank.ID),
                         data = droplevels(Comb_DB2[Comb_DB2$Trt.ID %in% c(Ref_ID, Comp_HR),]))
      }
      if(PH_Mod == "GEE"){
        cox_comp = coxph(Surv(Std_Time, Status) ~ Trt.ID, cluster = Tank.ID,
                         data = droplevels(Comb_DB2[Comb_DB2$Trt.ID %in% c(Ref_ID, Comp_HR),]))
      }
      pred_HR = exp(coef(cox_comp))

      #Get SR
      if(Method == 1) {
        Cumhaz1 = -(AUC(x = Ref_bshaz$time[Ref_bshaz$time < max(SR_Days)],
                        y = Ref_bshaz$hazard[Ref_bshaz$time < max(SR_Days)] * pred_HR) +
                      last(Ref_bshaz$hazard[Ref_bshaz$time < max(SR_Days)]) * 0.5 * pred_HR +
                      first(Ref_bshaz$hazard[Ref_bshaz$time < max(SR_Days)]) * 0.5 * pred_HR)
        pred_SR = exp(Cumhaz1)
      }

      if(Method == 2) {
        New_DB2 = survSplit(New_DB, cut = SR_Day, end = "Std_Time", event = "Status", episode = "Obs")
        New_DB2 = New_DB2[New_DB2$Obs == 1, ]
        Comp_bshaz = quiet(bshazard(data = droplevels(New_DB2[New_DB2$Trt.ID == Comp_HR,]),
                                    Surv(Std_Time, Status) ~ Tank.ID, verbose = FALSE, nbin = SR_Day))

        Cumhaz1 = -(AUC(x = Comp_bshaz$time,
                        y = Comp_bshaz$hazard) +
                      last(Comp_bshaz$hazard) * 0.5 +
                      first(Comp_bshaz$hazard) * 0.5)

        Cumhaz2 = ifelse(SR_Day < max(SR_Days) - 1,
                         -(AUC(x = Ref_bshaz$time[Ref_bshaz$time < max(SR_Days) & Ref_bshaz$time > SR_Day],
                               y = Ref_bshaz$hazard[Ref_bshaz$time < max(SR_Days) & Ref_bshaz$time > SR_Day] * pred_HR) +
                             last(Ref_bshaz$hazard[Ref_bshaz$time < max(SR_Days) & Ref_bshaz$time > SR_Day]) * 0.5 * pred_HR +
                             first(Ref_bshaz$hazard[Ref_bshaz$time < max(SR_Days) & Ref_bshaz$time > SR_Day]) * 0.5 * pred_HR),
                         -(last(Ref_bshaz$hazard[Ref_bshaz$time < max(SR_Days) & Ref_bshaz$time > SR_Day]) * 0.5 * pred_HR +
                             first(Ref_bshaz$hazard[Ref_bshaz$time < max(SR_Days) & Ref_bshaz$time > SR_Day]) * 0.5 * pred_HR))
        ifelse(is.na(Cumhaz2), Cumhaz2 <- 0, NA)
        pred_SR = exp(Cumhaz1 + Cumhaz2)
      }
      pred_SR_DB = rbind(pred_SR_DB, data.frame(Trt.ID = Comp_HR, Observable_SR_Day = SR_Day, pred_SR, pred_HR))
    }
  }
  row.names(pred_SR_DB) = NULL
  return(pred_SR_DB)
}

###################################################################################################################################
################################################## Function 6 - Surv_Gen() ########################################################

#' Title
#'
#' @param DB_Mort
#' @param Starting_Number_of_Fish_per_Tank
#'
#' @return
#' @export
#'
#' @examples
#'
Surv_Gen = function(mort_db = db_mort_ex,
                    today_tte,
                    tank_without_mort,
                    trt_without_mort,
                    starting_fish_count) {

  library(dplyr)
  DB_Mort_Gensum = data.frame(mort_db %>%
                                dplyr::group_by(Trt.ID, Tank.ID) %>%
                                dplyr::summarise(Num_dead = n()))

  WM_DB = data.frame(Trt.ID = trt_without_mort,
                     Tank.ID = tank_without_mort,
                     Num_dead = 0)
  DB_Mort_Gensum = rbind(DB_Mort_Gensum, WM_DB)

  DB_Mort_Gensum$Num_alive = starting_fish_count - DB_Mort_Gensum$Num_dead
  DB_Mort_Genalive = data.frame(lapply(DB_Mort_Gensum, rep, DB_Mort_Gensum$Num_alive))
  DB_Mort_Genalive$Status = 0
  DB_Mort_Genalive$TTE = today_tte
  DB_Mort_Gencomb = plyr::rbind.fill(mort_db, DB_Mort_Genalive[,-c(3:4)])

  return(DB_Mort_Gencomb)
}

###################################################################################################################################
################################################# Function 7 - Surv_Plots() #######################################################

#' Title
#'
#' @param surv_db
#' @param figure_name_prefix
#' @param x_axis_limits
#' @param y_axis_limits
#' @param x_lab
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
#'
Surv_Plots = function(surv_db,
                      figure_name_prefix = "figure_name_prefix",
                      x_axis_limits = c(0, max(surv_db$TTE)),
                      y_axis_limits = c(0, 1),
                      x_lab = "Days Post Challenge",
                      lambda = NULL) {

  library(ggplot2)

  surv_obj = survival::survfit(survival::Surv(TTE, Status) ~ Trt.ID, data = surv_db)
  attributes(surv_obj$strata)$names <- levels(as.factor(surv_db$Trt.ID))

  surv_plot = survminer::ggsurvplot(surv_obj,
                                    conf.int = FALSE,
                                    ggtheme = theme(plot.background = element_rect(fill = "white")),
                                    break.y.by = 0.1,
                                    break.x.by = min(round(max(x_axis_limits) / 15), 1),
                                    xlim = x_axis_limits,
                                    ylim = y_axis_limits,
                                    xlab = x_lab,
                                    surv.scale = "percent")
  plot_a = surv_plot$plot + ggplot2::theme(legend.position = "right") + ggplot2::guides(color = guide_legend("Trt."))
  ggplot2::ggsave(paste(figure_name_prefix, "Survival Curve.tiff"), dpi = 300, width = 6, height = 4, plot = plot_a)

  Haz_list = list()
  for(Haz_Trt in levels(as.factor(surv_db$Trt.ID))) {
    surv_db_trt = surv_db[surv_db$Trt.ID == Haz_Trt,]
    if(length(levels(as.factor(surv_db_trt$Tank.ID))) > 1) {
      Haz_bs = bshazard::bshazard(nbin = max(surv_db$TTE),
                                  data = surv_db_trt,
                                  survival::Surv(TTE, Status) ~ Tank.ID,
                                  verbose = FALSE,
                                  lambda = lambda)
    } else {
      Haz_bs = bshazard::bshazard(nbin = max(surv_db$TTE),
                                  data = surv_db_trt,
                                  survival::Surv(TTE, Status) ~ 1,
                                  verbose = FALSE,
                                  lambda = lambda)
    }

    Haz_DB = data.frame(Hazard = Haz_bs$hazard,
                        Time = Haz_bs$time)
    Haz_list[[Haz_Trt]] = data.frame(Hazard = Haz_bs$hazard,
                                     Time = Haz_bs$time)
  }

  Haz_DB = dplyr::bind_rows(Haz_list, .id = "Trt.ID")
  plot_b = ggplot(data = Haz_DB, aes(x = Time, y = Hazard, color = Trt.ID)) +
    geom_line(linewidth = 1) +
    geom_point() +
    xlab("Days Post Challenge") +
    scale_x_continuous(breaks = seq(from = min(x_axis_limits),
                                    to = max(x_axis_limits),
                                    by = min(round(max(x_axis_limits) / 15), 1)),
                       limits = x_axis_limits)

  ggplot2::ggsave(paste(figure_name_prefix, "Hazard Curve.tiff"), dpi = 300, width = 6, height = 4, plot = plot_b)

  return(list(plot_a, plot_b))
}
