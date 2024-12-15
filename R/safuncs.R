# This is an R package containing useful functions for my work.

# Available functions with documentation:
# 1. Con_Simul() -- simulates contingency tables based on the multinomial distribution.
# 1b. Con_Simul_PR() -- calculates positive rates for statistical tests on contingency tables.
# 3. Simul_Surv() -- simulate survival data based on a reference hazard function, the specified hazard ratio(s), and inter-tank variation.
# 4. theme_Publication() -- ggplot theme for generating publication-ready plots.
# 6. Surv_Gen() -- generate rows of survivors given a starting number of fish per tank and data containing morts and sampled fish.
# 7. Surv_Plots() -- generate Kaplan-Meier survival curve and hazard curve from survival data.
# 8. GG_Color_Hue() -- returns the default colour codes assigned by ggplot to a given number of categorical groups (n)

# Available functions without documentation:
# 2. Simul_Con_MULT.FISH.ORD() -- simulates ordinal-distributed data across treatments and lesions with inter-fish variation in the PO.
# 5. Surv_Pred() -- predict future survival rate(s) for ongoing experiment based on a reference hazard function from older data.
# 9. Label_Gen() -- generate combination of strings to use as labels.

########################################################## Function 1 - Con_Simul() ########################################################

#' @title Simulate Contingency Table
#'
#' @description Simulate a contingency table with fish counts distributed across \emph{n} lesion categories and \emph{n} treatment groups. Probability values for generating counts in each cell (i.e. each factor level combination) can be assigned using the \code{probs} argument. This function is designed for use in power and/or false positive rate calculations; for details, see \code{Con_Simul_PR()}.
#'
#' @details Counts are simulated from a multinomial distribution using \code{rmultinom()}. Counts may be assumed to have a fixed total in the marginals (e.g. per treatment group) or no fixed total in row or column marginals.
#'
#' For further discussion into the types of marginals in contingency tables, refer to: \href{https://www.uvm.edu/~statdhtx/StatPages/More_Stuff/Chi-square/Contingency-Tables.pdf}{here} and the comments on \bold{Arguments}.
#'
#' @param probs Matrix of probability values created using \code{matrix()}. Each row in the matrix should represent a treatment group and each column a lesion category. All probability values in the matrix should sum to 1. Default = equal probability across all cells.
#' @param total_count Total number of counts in the contingency table. Defaults to 750.
#' @param n_lesion Number of lesion categories. Ignored if \code{probs} specified. Defaults to 3.
#' @param n_Trt. Number of treatment groups. Ignored if \code{probs} specified. Defaults to 5.
#' @param margin_fixed_Trt. Whether margins are fixed per treatment group (i.e. fixed number of fish per treatment). Default = FALSE. See \bold{Details} for further information on marginals.
#' @param verbose Whether to print the parameters and probability matrix used. Default = TRUE.
#'
#' @seealso \href{https://sean4andrew.github.io/safuncs/reference/Con_Simul.html}{Link} for executed \bold{Examples} which includes any figure outputs.
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
#' # Simulate table with uniform probabilities across cells
#' Con_Simul(total_count = 750, n_lesion = 3, n_Trt. = 5)
#'
#' # Simulate table with specified probabilities across cells
#' Con_Simul(probs = matrix(nrow = 2, ncol = 3, c(1/6, 3/12, 1/6, 1/6, 1/6, 1/12)))
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

################################################## Function 1b - Con_Simul_PR() #################################################

#' @title Calculate Positive Rates for Contingency Table
#'
#' @description Computes statistical power and optionally false positive rates for tests applied to contingency tables based on simulations. Specify the simulation process using \code{Con_Simul()}, which serves as input. Positive rates are computed for the Chi-square test and optionally for Fisher's exact test and the Wald test applied to an ordinal regression model.
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
#' @seealso \href{https://sean4andrew.github.io/safuncs/reference/Con_Simul_PR.html}{Link} for executed \bold{Examples} which includes any figure outputs.
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
#' # Below I show how we can perform a simple power calculation using this tool.
#' # Suppose I want to calculate power for Treatment B which halves the lesions in
#' # category 2 and 3. I then specify the following probability matrix and feed it into
#' # Con_Simul():
#' probs_mat = matrix(nrow = 2, ncol = 3, data = c(1/6, 1/3, 1/6, 1/12, 1/6, 1/12))
#' sim_tab = Con_Simul(probs_mat)
#'
#' # Next, I feed the output into Con_Simul_PR():
#' Con_Simul_PR(sim_tab, sample_sizes = c(50, 100, 150))
#' # Results: Power is ~55, 86, and 97% for the Chi-square test using total counts of
#' # 50, 100, and 150, respectively.
#'
#' # The same power for Chi-square test can be calculated using Cohen's omega (w) method
#' # which is faster but has its own limitations; e.g. assumes one data generating
#' # process for the contingency table (the no fixed marginals).
#' library(pwr)
#' pwr::pwr.chisq.test(w = ES.w2(probs_mat), df = 2, sig.level = 0.05, N = 100)
#' # Results: Power is 85.6% for the Chi-square test at the total count of 100.
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

################################################### Function 3 - Surv_Simul() #####################################################

#' @title Simulate Survival Data
#'
#' @description Simulates survival data based on a set of user-specified experimental parameters and a reference hazard curve (e.g. hazard curve from the control group). Able to simulate data with inter-cluster (e.g. tank) variation which is added based on the framework of the mixed cox proportional hazards model (\code{coxme::coxme}). Able to simulate right censored data (e.g. sampled fish) using \code{sampling_specs} argument. Able to simulate treatment- and tank- specific fish numbers. Optionally produces a plot illustrating the characteristics of the simulated data and that of the population / truth from which the data (sample) is simulated.
#'
#' @details Simulations are based on uniform-probability draws (\emph{U} ~ (0, 1)) from a set of events which can be expressed as a function of time using the cumulative density function of failures (\emph{F(t)}, i.e. cumulative mort. curve). Because the cumulative mort curve (\emph{F(t)}) can be expressed in terms of the cumulative hazard function \emph{H(t)}, the relationship between \emph{H(t)} and \emph{U} draws is known (for derivation and equation, see \href{https://epub.ub.uni-muenchen.de/1716/1/paper_338.pdf}{Bender et al. (2003)}. Because \emph{H(t)} is related (as the integral) to the hazard function \emph{h(t)}, and since \emph{h(t)} is related to effects (e.g. treatment or tank) based on the cox proportional hazards model, such effects can now be incorporated into the simulation process as they interact with \emph{U}. The simulation process is as follows:
#'
#' \enumerate{
#' \item \code{Surv_Simul()} takes a random sample from \emph{U} (e.g. 0.7).
#'
#' \item U is then transformed into \emph{H} as they are related as discussed. The equation relating \emph{U}, \emph{H}, and treatment effects is shown below (obtain from \href{https://epub.ub.uni-muenchen.de/1716/1/paper_338.pdf}{Bender et al. 2003} which also shows the derivation of the equation):
#'
#' \emph{H = -log(U) ⋅ exp(-log(β))}; \emph{β} representing treatment or tank effects.
#'
#' \item The function \emph{H(t)} inverse (known from the supplied reference hazard curve) is applied to \emph{H} to obtain \emph{t} (time to event) which represents the survival data.
#'
#' \item Data with \emph{t} beyond the last follow-up period represent survivors (Status set to 0), and below it, represents mortalities (Status set to 1).
#' }
#'
#' To verify the correct "randomness" is produced in the simulated survival data, given that adding "randomness" is the whole point of simulations (to me), 5 different validation checks have been performed (documented in a pdf to be uploaded to github). Those checks showed that the HR estimated by fitting two curves sampled from the same population, converges to a mean of 1 (as should be) over many simulations, and across simulations the HR varies as expected (SD of simulated HRs = SE of HR as supposed by the cox model). Those checks also showed that the p-value obtained by applying log-rank test to null (no-effect) simulated datasets, has a distribution that is uniform (as should be), with a false positive rate of 0.05 given the alpha used was 0.05 (as should be). Additionally, power calculated from the simulations equal to the power calculated from an \href{https://homepage.univie.ac.at/robin.ristl/samplesize.php?test=logrank}{online calculator}. Last, the checks showed that the variations in a simulated survival curve is similar to that observed in curves simulated using a different, more limited, method (bootstraping / re-sampling with replacement).
#'
#' @param haz_db A dataframe representing the reference hazard curve; can be generated from \code{bshazard::bshazard()} or \code{Surv_Plots()}.
#' @param fish_num_per_tank The number of fish to simulate per tank, defaults to 100. If this differs by treatment, specify a vector of numbers ordered according to \code{treatments_hr}. When there is a need to compare experiments with different setups (fish numbers), specify the different setups as elements in a list (see \bold{Examples}). This is useful for comparing power between experimental setups (for calculations see \code{Surv_Power()}). Only 1 input parameter in \code{Surv_Simul()} can be specified as a list.
#' @param tank_num_per_trt The number of tanks to simulate per treatment group, defaults to 4. If this differs by treatment, specify a vector of numbers ordered according to \code{treatments_hr}. Input can be specified as elements in a list, with each element representing different experimental setups as described for \code{fish_num_per_tank}.
#' @param treatments_hr A vector representing the hazard ratios of the treatment groups starting with the reference/control (HR = 1), defaults to \code{c(1, 1, 1, 1)}. Length of the vector represents the number of treatment groups. Input can be specified as elements in a list, with each element representing different experimental setups as described for \code{fish_num_per_tank}.
#' @param logHR_sd_intertank The standard deviation of inter-tank variation (which contributes to overall data variation) in the log-HR scale according to the \code{coxme} framework. Defaults to 0 (no tank effect) which has been and quite oftenly, the estimate for injected Trojan fish data. For reference 0.1 reflects a low tank effect, while 0.35 is fairly high but can and has occurred in some immersion challenged fish datasets. Input can be specified as elements in a list, with each element representing different experimental setups as described for \code{fish_num_per_tank}.
#' @param sampling_specs A dataframe containing at least 2 columns; "Amount" representing the number of right censored data (e.g. sampled fish) per tank; "TTE" representing the time the sampling occurred; optionally a "Trt.ID" column to account for different sampling conditions per tank per treatment. Trt.IDs must start with "Control", then capitalized letters (see \bold{Examples}). Defaults to NULL (no sampling). Input can be specified as elements in a list, with each element representing different experimental setups as described for \code{fish_num_per_tank}.
#' @param exp_design A string specifying the type of experimental design. Can be "between-tank" which indicates each tank has a unique treatment hence the treatment effect occurs "between-tanks". Or, "within-tank" where each tank contains fish exposed to various treatments.
#' @param prog_notes Whether to print the number of simulations completed. Defaults to TRUE.
#' @param n_sim Number of survival dataset to simulate. Defaults to 1.
#' @param plot_out Whether to output the information plot (further details in \bold{Value}). Defaults to TRUE.
#' @param pop_out Whether to output a dataframe containing the survival probability values for the population. Defaults to TRUE.
#' @param theme A string specifying the graphics theme for the plots. Theme "ggplot2" and "prism" currently available. Defaults to "ggplot2".
#' @param plot_save Whether to save plot as .tiff in the working directory. Defaults to TRUE.
#'
#' @return Returns a list that, at minimum, contains the simulated survival dataframe which has 5 columns: TTE (Time to Event), Status (0 / 1), Trt.ID, Tank.ID, and n_sim which represents the simulation number for the data subsets.
#'
#' If \code{plot_out = TRUE}, the list additionally contains a Kaplan-Meier survival plot. The plot illustrates the survival curve with end survival rates for the simulated dataset as well as the population. If the number of simulated dataset is greater than 1, multiple curves are drawn representing each and a statement is provided regarding the power to detect the effect of Treatment -- specifically, the percent positive (p < 0.05) from a global log-rank test using \code{survival::survdiff()}.
#'
#' If \code{pop_out = TRUE}, the list additionally contains a dataframe representing the survival probability values for the population / truth from which the sample is supposedly taken.
#'
#' @import ggplot2
#' @import magrittr
#' @import dplyr
#'
#' @seealso \href{https://sean4andrew.github.io/safuncs/reference/Surv_Simul.html}{Link} for executed \bold{Examples} which includes any figure outputs.
#'
#' @export
#'
#' @examples
#' # Starting from an example mortality database, we first generate the complete survivor
#' # data using Surv_Gen()
#' data(mort_db_ex)
#' surv_dat = Surv_Gen(mort_db = mort_db_ex,
#'                     starting_fish_count = 100,
#'                     last_tte = 54)
#'
#' # Filter for the control group ("A") to get a reference hazard curve for simulations
#' surv_dat_A = surv_dat[surv_dat$Trt.ID == "A", ]
#'
#' # Estimate the hazard curve of the control group and get the associated hazard
#' # dataframe using bshazard::bshazard() or safuncs::Surv_Plots()$Hazard_DB
#' ref_haz_route_bshazard = bshazard::bshazard(data = surv_dat_A,
#'                                             survival::Surv(TTE, Status) ~ Tank.ID,
#'                                             nbin = max(surv_dat_A$TTE),
#'                                             verbose = FALSE)
#' ref_haz_route_bshazard = data.frame(summary(ref_haz_route_bshazard)$HazardEstimates)
#'
#' ref_haz_route_safuncs = safuncs::Surv_Plots(surv_db = surv_dat_A,
#'                                             data_out = TRUE)$Hazard_DB
#'
#' # Simulate! Sampled 10 fish per tank at 45 DPC, but otherwise default conditions.
#' Surv_Simul(haz_db = ref_haz_route_safuncs,
#'            treatments_hr = c(1, 0.8, 0.5),
#'            sampling_specs = data.frame(Amount = 10,
#'                                        TTE = 45))$surv_plots
#'
#' # Further, results of simulating multiple times are shown to better understand the
#' # chance that future samples accurately capture the truth/population. Specify n_sim!
#' Surv_Simul(haz_db = ref_haz_route_safuncs,
#'            treatments_hr = c(1, 0.8, 0.5),
#'            sampling_specs = data.frame(Amount = 10,
#'                                        TTE = 45),
#'            prog_notes = FALSE, #hide simulation progress notes for cleaner output
#'            n_sim = 4)$surv_plots
#'
#' # Surv_Simul() can handle even more complicated experimental designs. Below, I use
#' # different (across treatments) fish numbers per tank, tank numbers, and sampling
#' # designs.
#' Surv_Simul(haz_db = ref_haz_route_safuncs,
#'            fish_num_per_tank = c(50, 100, 100), #for Ctrl., Trt.A, B, respectively
#'            tank_num_per_trt = c(1, 1, 2),       #Ctrl., A, B
#'            treatments_hr = c(1, 0.8, 0.5),      #Ctrl., A, B
#'            sampling_specs = data.frame(TTE = c(20, 40, 50),
#'                                        Amount = c(0, 20, 5), #0 sample for Ctrl.
#'                                        Trt.ID = c("Control", "A", "B")),
#'            prog_notes = FALSE,
#'            n_sim = 4)$surv_plots
#'
#' # What if we want to compare power of the global log-rank test (shown in the plot)
#' # across different experimental setups with different fish numbers per treatment?
#' # Below, I setup a Surv_Simul() to answer this question.
#' Surv_Simul(haz_db = haz_db_ex,
#'            fish_num_per_tank = list(30, 100),
#'            tank_num_per_trt = 3,
#'            treatments_hr = c(1, 0.6),
#'            prog_notes = FALSE,
#'            n_sim = 30)$surv_plots
#'
#' # Plot[[1]] and [[2]] shows the results from fish_num_per_tank = 30 and 100,
#' # respectively. Additionally, the simulated data output (...$surv_simul_db) can be
#' # supplied to safuncs::Surv_Power() (under development) to calculate power for
#' # various tests (e.g. log-rank global, pairwise with(out) correction) or tests based
#' # on statistical models with various forms (e.g. with(out) tank-variation)).
Surv_Simul = function(haz_db,
                      fish_num_per_tank = 100,
                      tank_num_per_trt = 4,
                      treatments_hr = c(1, 1, 1, 1),
                      logHR_sd_intertank = 0,
                      sampling_specs = NULL,
                      exp_design = "between-tank",
                      n_sim = 1,
                      prog_notes = TRUE,
                      plot_out = TRUE,
                      pop_out = TRUE,
                      theme = "ggplot2",
                      plot_save = TRUE) {
  #Track time elapsed
  time_start = Sys.time()

  #Making sure input data has correct (lower case) column names
  colnames(haz_db) = tolower(colnames(haz_db))

  #Initialize objects to store second output type (across list elements and loops)
  output2 = list(surv_plots = list(), simul_surv_db = data.frame(), population_surv_db = data.frame())
  list_var_check = c()

  #Validation Check(s)
  if(plot_out == TRUE) {
    if(logHR_sd_intertank > 0) {
      print("NOTE: You specified a tank effect/contribution to variation, but the power shown in the plot is based of the logrank test. This test assumes no such tank effects. Adding tank-variation tends to decrease power of the logrank when the treatment effect is strong-modest (see Examples in Surv_Simul()'s documentation). Despite the decrease, power of the logrank will still be greater than that of other statistical tests which considers tank variation (assuming the experimental design is 'between-tank' because this matters). At the cost of having the greater power, logrank suffers from a greater false positive rate (> 5%). To calculate power of other tests that account for tank variation (hence keeping FPR at ~5%), use the coxph_glmm model option in the function Surv_Power() from package safuncs.")
    }
  }

  #Finding the input variable (var_name) that is a list and store info in var_list
  #First we stop the function if we find more than 1 list
  var_names = c("fish_num_per_tank", "tank_num_per_trt", "treatments_hr", "logHR_sd_intertank", "sampling_specs")
  for (var_name_check in var_names) {
    if(is.list(get(var_name_check, envir = environment())) & !is.data.frame(get(var_name_check, envir = environment()))){
      list_var_check = c(list_var_check, var_name_check)
    }
  }
  if(length(list_var_check) > 1) {stop("You specified more than 1 argument/parameter as a list. Currently, this is not allowed.")}

  for (var_name in var_names) {
    ifelse(is.list(get(var_name, envir = environment()))
           & !is.data.frame(get(var_name, envir = environment())),
           list_var <- get(var_name, envir = environment()),
           list_var <- "empty")
    if(length(list_var) > 1) {
      break
    }
  }

  #Track progress
  prog = 0

  #Change var_name based on list_var elements
  for (ele_num in 1:length(list_var)) {

    #if you have list elements, assign and print. If not just jump straight to old code
    if(length(list_var) > 1) { #if you have list elements, assign and print, otherwise just go to old code.
      assign(var_name, list_var[[ele_num]]) #assign
    }

    #Old code below. Will only run once if there is no list (i.e. length(list_var) = 1)

    #Initialize objects to store loop results
    surv_samps = data.frame() #for plotting purposes
    cens_db = data.frame() #for plotting purposes
    pvalues = c() #for plotting purposes
    Surv_simul_outDB = data.frame() #for dataoutput

    #Simulate survival dataframe
    for(loopnum in 1:n_sim) {

      CDF_Yval = c()
      Trt.ID = c()
      Tank.ID = c()
      Tank_num2 = 0
      iTT = 0

      if(exp_design == "between-tank") { #simulation procedure for between-tanks experimental design

        for(Treatment_Term in treatments_hr) {
          iTT = iTT + 1

          for(Tank_num in 1:ifelse(length(tank_num_per_trt) > 1, tank_num_per_trt[iTT], tank_num_per_trt)) {
            Tank_num2 = Tank_num2 + 1

            #Random sampling
            Tank_eff = rnorm(n = 1, mean = 0, sd = logHR_sd_intertank)

            U = runif(n = ifelse(length(fish_num_per_tank) > 1, fish_num_per_tank[iTT], fish_num_per_tank), min = 0, max = 1)

            CDF_Yval_temp = -log(U) * exp(-(log(Treatment_Term) + Tank_eff))
            CDF_Yval = append(CDF_Yval, CDF_Yval_temp)

            Trt.ID = c(Trt.ID, rep(c("Control", LETTERS[1:(length(treatments_hr) - 1)])[iTT], length(CDF_Yval_temp)))
            Tank.ID = c(Tank.ID, rep(Tank_num2, length(CDF_Yval_temp)))
          }
        }

      }
      if(exp_design == "within-tank") { #simulation procedure for within-tank experimental design. Similar but with flipped Tank-Trt loops.

        for(Tank_num in 1:tank_num_per_trt) { #only 1 tank num can be specified for the within-tank design.
          Tank_num2 = Tank_num2 + 1

          #Tank effect
          Tank_eff = rnorm(n = 1, mean = 0, sd = logHR_sd_intertank)

          iTT = 0
          for(Treatment_Term in treatments_hr) {
            iTT = iTT + 1

            #Simulate fish numbers per treatment for each tank. Write down in description that it must be per treatment per tank.
            #Can be as vector for treatment specific numbers or can be a single value if the same across treatments.
            U = runif(n = ifelse(length(fish_num_per_tank) > 1, fish_num_per_tank[iTT], fish_num_per_tank), min = 0, max = 1)

            CDF_Yval_temp = -log(U) * exp(-(log(Treatment_Term) + Tank_eff))
            CDF_Yval = append(CDF_Yval, CDF_Yval_temp)

            Trt.ID = c(Trt.ID, rep(c("Control", LETTERS[1:(length(treatments_hr) - 1)])[iTT], length(CDF_Yval_temp)))
            Tank.ID = c(Tank.ID, rep(Tank_num2, length(CDF_Yval_temp)))
          }
        }
      }

      #Get Time to Event
      TTE = approx(x = cumsum(haz_db$hazard), y = haz_db$time, xout = CDF_Yval, method = "linear")$y
      TTE = round(TTE, digits = 0)

      #Turn NA (from out of bound CDF_Yval) to the last follow up time
      TTE = ifelse(is.na(TTE), max(haz_db$time), TTE)

      #Label Status (1 - dead, or 0 - survived) given TTE, and create survival dataframe
      Surv_simul_DB = data.frame(TTE = TTE,
                                 Status = ifelse(TTE == max(haz_db$time), 0, 1),
                                 Trt.ID = Trt.ID,
                                 Tank.ID = Tank.ID,
                                 n_sim = loopnum)

      #Transform TTE and Status (to 0) in certain rows due to sampling
      if(!is.null(sampling_specs)) {

        if(!"Trt.ID" %in% colnames(sampling_specs)) {
          Trt_levels = unique(Surv_simul_DB$Trt.ID)
          sampling_specs = data.frame(TTE = rep(sampling_specs$TTE, each = length(Trt_levels)),
                                      Amount = rep(sampling_specs$Amount, each = length(Trt_levels)),
                                      Trt.ID = rep(unique(Trt_levels, times = nrow(sampling_specs))))
        }

        #Put in Tank.ID and replicate accordingly
        sampling_specs2 = merge(sampling_specs, unique(Surv_simul_DB[, 3:4]))

        #Run through every row of sampling_specs and sample accordingly
        for(samp_row in 1:nrow(sampling_specs2)) {

          rows_samp_space = which(Surv_simul_DB$Trt.ID == sampling_specs2$Trt.ID[samp_row] &
                                    Surv_simul_DB$Tank.ID == sampling_specs2$Tank.ID[samp_row] &
                                    Surv_simul_DB$TTE > sampling_specs2$TTE[samp_row])

          #Catch over sampling situation and print message
          if(length(rows_samp_space) < sampling_specs2$Amount[samp_row]) {
            print(paste(sep = "", "In simulation set-", loopnum, " Trt.ID-", sampling_specs2$Trt.ID[samp_row], ", Tank.ID-",
                        sampling_specs2$Tank.ID[samp_row],
                        ", you requested more samples than the fish alive! All remaining (living) fish sampled."))

            #Modify sampling amount
            sampling_specs2$Amount[samp_row] = length(rows_samp_space)
          }

          #Select rows that were sampled
          rows_sel = sample(x = rows_samp_space,
                            size = sampling_specs2$Amount[samp_row],
                            replace = FALSE)

          #Change Status and Time for sampled individuals
          Surv_simul_DB$TTE[rows_sel] = sampling_specs2$TTE[samp_row]
          Surv_simul_DB$Status[rows_sel] = 0
        }
      }

      #Get p-value for plots
      pvalues = append(pvalues, survival::survdiff(survival::Surv(TTE, Status) ~ Trt.ID, Surv_simul_DB)$pvalue)

      #Simulated survival data to be provided as output
      if(length(list_var) > 1){Surv_simul_DB$list_element_num <- ele_num}
      Surv_simul_outDB = rbind(Surv_simul_outDB, Surv_simul_DB)

      #Transform simulated survival data for plotting purposes
      surv_obj = survival::survfit(survival::Surv(TTE, Status) ~ Trt.ID, data = Surv_simul_DB)
      if(length(levels(as.factor(Surv_simul_DB$Trt.ID))) > 1) {
        attributes(surv_obj$strata)$names <- levels(as.factor(Surv_simul_DB$Trt.ID))
      } else {
        surv_obj$strata = length(surv_obj$surv)
        attributes(surv_obj$strata)$names <- levels(as.factor(Surv_simul_DB$Trt.ID))
      }

      surv_samps_temp = data.frame(Trt.ID = summary(surv_obj)$strata,
                                   surv_prob = summary(surv_obj)$surv,
                                   time = summary(surv_obj)$time,
                                   type = paste("Sample set (n = ", n_sim, ")", sep = ""),
                                   n_sim = loopnum,
                                   alpha = 1 - (0.0001 ^ (1/n_sim)))
      if(length(list_var) > 1){surv_samps_temp$list_element_num <- ele_num}

      surv_samps_ends = data.frame(surv_samps_temp %>%
                                     dplyr::group_by(Trt.ID) %>%
                                     dplyr::reframe(surv_prob = c(1, min(surv_prob)),
                                                    time = c(floor(min(haz_db$time)), ceiling(max(haz_db$time))),
                                                    n_sim = loopnum,
                                                    alpha = 1 - (0.0001 ^ (1/n_sim))))
      surv_samps_ends$type = paste("Sample set (n = ", n_sim, ")", sep = "")
      if(length(list_var) > 1){surv_samps_ends$list_element_num <- ele_num}

      surv_samps = rbind(surv_samps, surv_samps_temp, surv_samps_ends)

      if(!is.null(sampling_specs)){
        #Get survival probability at mid censoring
        cens_db_temp  = data.frame(Trt.ID = summary(surv_obj, time = sampling_specs$TTE)$strata,
                                   surv_prob = summary(surv_obj, time = sampling_specs$TTE)$surv,
                                   time = summary(surv_obj, time = sampling_specs$TTE)$time,
                                   n_sim = loopnum,
                                   type = as.factor(paste("Sample set (n = ", n_sim, ")", sep = "")))
        cens_db = rbind(cens_db, cens_db_temp)
      }

      #Print progress
      if(prog_notes != FALSE) {cat("\rSimulated", prog <- prog + 1, "of", n_sim * length(list_var), "sample sets")}
    } #close loopnum

    #Get "population" survival dataset by exponentiating the negative cumulative hazard
    pop_haz_db = data.frame(approx(x = haz_db$time, y = haz_db$hazard, xout = seq(min(haz_db$time), max(haz_db$time), 0.1), method = "linear"))
    colnames(pop_haz_db) = c("time", "hazard")

    #For use with old surv_prob method (revived)
    surv_pop = data.frame(Trt.ID = as.factor(rep(c("Control", LETTERS[1:(length(treatments_hr) - 1)]), each = length(haz_db$hazard))),
                          #cumhaz_prob = as.vector(apply((haz_db$hazard) %*% t(treatments_hr), 2, cumsum)),
                          surv_prob = exp(-as.vector(apply(haz_db$hazard %*% t(treatments_hr), 2, cumsum))),
                          time = rep(haz_db$time, times = length(treatments_hr)),
                          type = "Population / truth",
                          n_sim = 1,
                          alpha = 1)
    if(length(list_var) > 1){surv_pop$list_element_num <- ele_num}

    #To the end of creating survival plots
    surv_comb = rbind(surv_samps, surv_pop)
    surv_comb$type = factor(surv_comb$type, levels = c(paste("Sample set (n = ", n_sim, ")", sep = ""), "Population / truth"))
    surv_comb$Trt.ID = factor(surv_comb$Trt.ID, levels = rep(c("Control", LETTERS[1:(length(treatments_hr) - 1)])))

    #Get end_sr for population plots and sample plots
    end_db = data.frame(surv_comb %>%
                          dplyr::group_by(type, Trt.ID, n_sim) %>%
                          dplyr::summarise(surv_prob = min(surv_prob), time = max(TTE), .groups = "drop") %>%
                          dplyr::group_by(type, Trt.ID) %>%
                          dplyr::summarise(surv_prob = mean(surv_prob), time = max(time), .groups = "drop"))

    #Get % significance (i.e. power) for plotting
    perc_sf = paste(round(100 * sum(pvalues < 0.05) / length(pvalues), digits = 0), "%", sep = "")

    #Ggplot
    surv_plots = ggplot(data = surv_comb, aes(x = time, y = surv_prob, colour = Trt.ID, group = interaction(n_sim, Trt.ID))) +
      facet_wrap(~ type) +
      geom_step(aes(alpha = alpha)) +
      scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), labels = scales::percent) +
      scale_x_continuous(breaks = seq(0, max(surv_pop$time), max(round(max(surv_pop$time) / 12), 1))) +
      ylab("Survival Probability (%)") +
      xlab("Time to Event") +
      scale_alpha(range = c(min(surv_comb$alpha), 1)) +
      guides(alpha = "none") +
      coord_cartesian(clip = "off") +
      theme(plot.margin = margin(5.5, 20, 5.5, 5.5))

    if(n_sim == 1) {
      surv_plots = surv_plots +
        geom_text(data = end_db, aes(x = time, y = surv_prob, label = round(surv_prob * 100, digits = 0)),
                  vjust = -0.3, hjust = 0.8, show.legend = FALSE, size = 3.3) +
        annotation_custom(grob = grid::textGrob(paste(c(paste("The sample has a", sep = ""),
                                                        paste("p-value = ", signif(pvalues, digits = 2), sep = ""),
                                                        "(global test of Trt.)"), collapse = "\n"),
                                                x = grid::unit(1.05, "npc"),
                                                y = grid::unit(0.08, "npc"),
                                                hjust = 0,
                                                gp = grid::gpar(fontsize = 9)))
    } else {
      surv_plots = surv_plots +
        geom_text(data = end_db[end_db$type == "Population / truth",],
                  aes(x = time, y = surv_prob, label = round(surv_prob * 100, digits = 0)),
                  vjust = -0.3, hjust = 0.8, show.legend = FALSE, size = 3.3) +
        annotation_custom(grob = grid::textGrob(paste(c(paste(perc_sf, " of the sample", sep = ""),
                                                        paste("sets (n) has p<0.05", sep = ""),
                                                        "(global test of Trt.)"), collapse = "\n"),
                                                x = grid::unit(1.03, "npc"),
                                                y = grid::unit(0.08, "npc"),
                                                hjust = 0,
                                                gp = grid::gpar(fontsize = 9)))
    }

    #Add censoring points
    if(!is.null(sampling_specs)) {
      merged_db = merge(sampling_specs2, Surv_simul_DB)
      merged_db = merged_db[merged_db$Status == 0, ]
      cens_db = cens_db[interaction(cens_db$Trt.ID, cens_db$time) %in% interaction(merged_db$Trt.ID, merged_db$TTE),]
      surv_plots = surv_plots +
        geom_point(data = cens_db, aes(x = time, y = surv_prob, colour = Trt.ID), shape = 3, size = 0.7, stroke = 1)
    }

    #Plot theme
    if(theme == "prism") {surv_plots = surv_plots + ggprism::theme_prism()}

    #Plot title
    if(length(list_var > 1)) {

      surv_plots = surv_plots + labs(title = paste("List Element", ele_num))
    }

    #Save plot
    if(plot_save == TRUE){
      ggsave(paste("Simul_Surv_Plot",
                   ifelse(length(list_var) == 1, "_", paste("_Element", ele_num, "_", sep ="")),
                   Sys.Date(), ".tiff", sep = ""), dpi = 900, width = 7, height = 4, plot = surv_plots)
    }

    #remove "alpha" column from data output
    surv_pop = surv_pop[, -6]

    #Return R output if list_var length = 1 (i.e. no list)
    if(length(list_var) == 1) {
      if(plot_out == FALSE & pop_out == FALSE) {
        #Print time elapsed
        print(paste("Time elapsed:", substr(hms::as_hms(Sys.time() - time_start), 1, 8), "(hh:mm:ss)"))
        return(simul_surv_db = Surv_simul_outDB)

      } else {

        output = list(simul_surv_db = Surv_simul_outDB)

        if(pop_out == TRUE) {output$population_surv_db <- surv_pop}
        if(plot_out == TRUE) {output$surv_plots <- surv_plots}

        #Print time elapsed
        print(paste("Time elapsed:", substr(hms::as_hms(Sys.time() - time_start), 1, 8), "(hh:mm:ss)"))
        return(output)
      }
    }

    if(length(list_var) > 1) {

      #Store 2nd output if list_var length >1
      output2$surv_plots[[ele_num]] = surv_plots
      output2$simul_surv_db = rbind(output2$simul_surv_db, Surv_simul_outDB)
      output2$population_surv_db = rbind(output2$population_surv_db, surv_pop)
    }

    #Old code (non-lists stuff) ends here

  } #This closes the loop that deals with lists

  if(length(list_var) > 1){

    if(plot_out == FALSE) {
      output2$surv_plots = NULL
    }

    if(pop_out == FALSE) {
      output2$population_surv_db = NULL
    }

    #Print time elapsed
    print(paste("Time elapsed:", substr(hms::as_hms(Sys.time() - time_start), 1, 8), "(hh:mm:ss)"))
    return(output2)
  }
}

#################################################### Function 4 - theme_Publication() ##################################################

#' @title Publication theme for ggplot2
#'
#' @description A theme function to add to a ggplot2 object for publication style plots. Function adapted from \href{https://rdrr.io/github/HanjoStudy/quotidieR/src/R/theme_publication.R}{HanjoStudy/quotidieR}.
#'
#' @param base_size size of text in graph
#'
#' @return theme function to add to a ggplot2 object
#'
#' @import grid
#' @import dplyr
#' @import ggthemes
#'
#' @seealso \href{file:///C:/Users/sean4/Documents/GitHub/safuncs/docs/reference/theme_Publication.html}{Link} for executed \bold{Examples} which includes any figure outputs.
#'
#' @export
#'
#' @examples
#' # Load an example dataset
#' data(iris)
#'
#' # Create a ggplot modified with theme_Publication()
#' library(ggplot2)
#' ggplot(data = iris, aes(x = Species, colour = Species, y = Petal.Length)) +
#'    geom_boxplot() +
#'    theme_Publication()
theme_Publication = function(base_size = 14) {

  (theme_foundation(base_size = base_size)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle= 90,vjust = 2),
            axis.title.x = element_text(vjust = -0.6),
            axis.text = element_text(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_blank(),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.4, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(size = 12),
            plot.margin = unit(c(10, 15, 5, 5),"mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
    ))
}

################################################## Function 5 - Surv_Pred() #######################################################

#' @title Predict End Survival Rate
#'
#' @description Predict survival rate for a given survival dataset provided a reference survival database used to estimate a reference hazard curve. Prediction done separately by treatment group.
#'
#' @details P
#'
#' @param pred_db Placeholder
#' @param ref_db Placeholder
#' @param predsr_tte Placeholder
#' @param method Placeholder
#' @param coxph_mod Placeholder
#' @param lambda_pred Placeholder
#' @param phi_pred Placeholder
#'
#' @return Placeholder
#'
#' @import dplyr
#' @import ggplot2
#'
#' @seealso \href{https://sean4andrew.github.io/safuncs/reference/Surv_Pred.html}{Link} for executed \bold{Examples} which includes any figure outputs.
#'
#' @export
#'
#' @examples
#' #Placeholder
Surv_Pred = function(pred_db, #Data from ongoing study, with SR to be predicted. See \bold{Details} for specifics.
                     ref_db, #Reference survival data from to create the reference hazard function.
                     predsr_tte, #The day at which SR is to be predicted. Minimum is Day 5 post challenge.
                     method = 2, #SR prediction method. See \bold{Details} for more info.
                     coxph_mod = "GLMM", #Model used to estimate HR. Can be either "GLMM" or "GEE". See \bold{Details} for more info.
                     lambda_pred = NULL, #Lambda parameter for the bshazard curve of the predicted dataset.
                     phi_pred = NULL)
{
  #Ensure positive TTE
  pred_db = pred_db[pred_db$TTE > 0, ]
  ref_db = ref_db[ref_db$TTE > 0, ]

  #Get reference level hazard curve
  ref_id = levels(as.factor(ref_db$Trt.ID))
  ref_db$Trt.ID = paste("ref", ref_id)
  ref_bshaz = bshazard::bshazard(data = ref_db, survival::Surv(TTE, Status) ~ Tank.ID,
                                 verbose = FALSE)

  #Initialize dataframes
  pred_SR_DB = data.frame()
  haz_db = data.frame()

  #Loop for every treatment
  for(pred_trt in levels(as.factor(pred_db$Trt.ID))) {

    #Combine reference and sample hazard dataframes
    comb_db = rbind(ref_db, pred_db[pred_db$Trt.ID == pred_trt,])
    comb_db$Trt.ID = relevel(as.factor(comb_db$Trt.ID), ref = paste("ref", ref_id))

    #Get HR and SR for every predictable day
    for(SR_Day in 5:(predsr_tte-1)) {
      comb_db2 = survival::survSplit(comb_db, cut = SR_Day, end = "TTE", event = "Status", episode = "Obs")
      comb_db2 = comb_db2[comb_db2$Obs == 1, ]

      #HR Calculation
      if(coxph_mod == "GLMM"){cox_comp = coxme::coxme(survival::Surv(TTE, Status) ~ Trt.ID + (1|Tank.ID),
                                                        data = comb_db2)}
      if(coxph_mod == "GEE"){cox_comp = survival::coxph(survival::Surv(TTE, Status) ~ Trt.ID, cluster = Tank.ID,
                                                          data = comb_db2)}
      pred_HR = exp(coef(cox_comp))

      #SR Calculation
      ref_bshaz_t = data.frame(hazard = ref_bshaz$hazard[ref_bshaz$time < predsr_tte],
                                 time = ref_bshaz$time[ref_bshaz$time < predsr_tte])

        #Uses the shape of the reference hazard curve throughout
        if(method == 1) {
          cumhaz = DescTools::AUC(x = c(ref_bshaz_t$time),
                                  y = c(ref_bshaz_t$hazard) * pred_HR)
          pred_SR = 100 * exp(-cumhaz)
        }

        #Uses the shape of the reference and supplied/observable hazard curve
        if(method == 2) {

          #Split survival data for calculating cumulative hazard using the reference and the observable hazard curve
          pred_db2 = survival::survSplit(pred_db, cut = SR_Day, end = "TTE", event = "Status", episode = "Obs")
          pred_db2 = pred_db2[pred_db2$Obs == 1, ]

          #observable hazard curve
          pred_bshaz = bshazard::bshazard(data = droplevels(pred_db2[pred_db2$Trt.ID == pred_trt,]),
                                          survival::Surv(TTE, Status) ~ Tank.ID, verbose = FALSE, lambda = lambda_pred)
          cumhaz_precut = DescTools::AUC(x = c(pred_bshaz$time),
                                         y = c(pred_bshaz$hazard))

          #reference hazard curve
          ref_bshaz_t2 = ref_bshaz_t[ref_bshaz_t$time > SR_Day,]
          cumhaz_postcut = DescTools::AUC(x = c(ref_bshaz_t2$time),
                                          y = c(ref_bshaz_t2$hazard) * pred_HR)

          if(is.na(cumhaz_postcut)) {cumhaz_postcut <- 0}
          pred_SR = 100 * exp(-(cumhaz_precut + cumhaz_postcut))
        }

      #Get prediction database
      pred_SR_DB = rbind(pred_SR_DB, data.frame(Trt.ID = pred_trt,
                                                Observable_SR_Day = SR_Day,
                                                pred_SR,
                                                pred_HR))
    }
  }

  row.names(pred_SR_DB) = NULL

  #Create plots representing survival and hazard predictions over time
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

  return(list(survival_prediction = pred_SR_DB,
              surv_pred_time_plot = Pred_SR_Plot,
              hr_pred_time_plot = Pred_HR_Plot))
}

################################################## Function 6 - Surv_Gen() ########################################################

#' @title Generate Survivor Data
#'
#' @description Produces survival data that includes rows for every surviving fish based on the starting number of fish and mortality data. To generate survivor data for tanks absent in the input mortality dataframe, specify the arguments \code{tank_without_mort} and \code{trt_without_mort}. To generate survivor data with tank specific starting numbers of fish, input a dataframe into the argument \code{starting_fish_count} instead of a single value; details in \bold{Arguments}.
#'
#' @details The mort dataframe supplied as input should consist of the following 4 columns at minimum:
#' * "Trt.ID" = Labels for treatment groups in the study.
#' * "Tank.ID" = Labels for tanks in the study (each tank must have a unique label).
#' * "TTE" = Time to Event. Event could be fish death or being sampled and removed depending on "Status".
#' * "Status" = Value indicating what happened at TTE. 1 for dead fish, 0 for those sampled and removed.
#'
#' Each row should represent one fish.
#'
#' For an example dataframe, execute \code{data(mort_db_ex)} and view.
#' @md
#'
#' @param mort_db A mort dataframe as described in \bold{Details}.
#' @param starting_fish_count Value representing the starting number of fish for every tank. Alternatively, a dataframe containing the columns "Trt.ID", "Tank.ID", and "starting_fish_count" to allow for different fish starting numbers per tank.
#' @param last_tte Value representing the time-to-event the fish survived to, assigned to every row of survivor data generated.
#' @param tank_without_mort A vector of strings specifying the tanks absent from \code{mort_db}; used to generate survivor data for those tanks. Argument ignored if \code{starting_fish_count} is a dataframe.
#' @param trt_without_mort A vector of strings corresponding to \code{tank_without_mort}. Keep their order the same. Argument ignored if \code{starting_fish_count} is a dataframe.
#' @param output_prism Whether to generate and save a prism ready survival csv. Defaults to FALSE.
#' @param output_prism_date The starting date to be used in the prism file. Please specify date in "dd-Mmm-yyyy" syntax (e.g. "08-Aug-2024").
#'
#' @return A dataframe produced by combining the input mort data and generated rows of survivor data.
#'
#' @import magrittr
#' @import devtools
#'
#' @seealso \href{https://sean4andrew.github.io/safuncs/reference/Surv_Gen.html}{Link} for executed \bold{Examples} which includes any figure outputs.
#'
#' @export
#'
#' @examples
#' # First, we load an example mortality database available from the safuncs package
#' data(mort_db_ex)
#'
#' # Next, we input this data into Surv_Gen() as well as the study details to generate
#' # entries (rows) for survivors in the output - a "complete" dataframe for further
#' # survival analysis and data visualization.
#' Surv_Data_Output = Surv_Gen(mort_db = mort_db_ex,
#'                    starting_fish_count = 100,
#'                    last_tte = 54,
#'                    tank_without_mort = c("C99", "C100"),
#'                    trt_without_mort = c("A", "B"))
#'
#' # Below, the bottom 5 rows of the output is displayed to show the rows of survivor
#' # data generated.
#' tail(Surv_Data_Output, n = 5)
#'
#' # Below is another example, this time showing how to specify tank-specific fish
#' # numbers. First, create the database containing information on the starting fish
#' # counts for the different tanks. You must include all tanks that are in your mort
#' # database to get a proper output. For the example below, we will later trim down
#' # the mort database to only 4 tanks for simplicity. Use these 3 column names:
#' # starting_fish_count, Tank.ID and Trt.ID.
#' count_db = data.frame(starting_fish_count = c(100, 100, 120, 120),
#'                       Tank.ID = c("C1", "C6", "C5", "C8"),
#'                       Trt.ID = c("B", "B", "D", "D"))
#'
#' filtered_mort_db = mort_db_ex[mort_db_ex$Tank.ID %in% c("C1", "C6", "C5", "C8"),]
#'
#' # We then use 'count_db' as input to the argument 'starting_fish_count' in Surv_Gen():
#' Surv_Data_Output = Surv_Gen(mort_db = filtered_mort_db,
#'                             starting_fish_count = count_db,
#'                             last_tte = 54)
Surv_Gen = function(mort_db,
                    starting_fish_count,
                    last_tte,
                    tank_without_mort = NULL,
                    trt_without_mort = NULL,
                    output_prism = FALSE,
                    output_prism_date = NULL) {

  #Remove NA rows
  mort_db[mort_db == "#N/A"] = NA
  mort_db = na.omit(mort_db)

  #Count the number of rows in mort_db, for each combination of treatment and tank ID
  DB_Mort_Gensum = data.frame(mort_db %>%
                                dplyr::group_by(Trt.ID, Tank.ID) %>%
                                dplyr::summarise(Num_dead = dplyr::n()))

  #Include tanks without morts in the count database
  if(!is.null(tank_without_mort) && !is.null(trt_without_mort)) {
    WM_DB = data.frame(Trt.ID = trt_without_mort,
                       Tank.ID = tank_without_mort,
                       Num_dead = 0)
    DB_Mort_Gensum = rbind(DB_Mort_Gensum, WM_DB)
  }

  #Use tank-specific starting fish count if information provided
  if(is.data.frame(starting_fish_count)) {
    DB_Mort_Gensum = base::merge(DB_Mort_Gensum, starting_fish_count, all.y = TRUE)
    DB_Mort_Gensum$Num_dead[is.na(DB_Mort_Gensum$Num_dead)] = 0
    DB_Mort_Gensum$Num_alive = DB_Mort_Gensum$starting_fish_count - DB_Mort_Gensum$Num_dead
    DB_Mort_Gensum = DB_Mort_Gensum[, -3]
  } else {DB_Mort_Gensum$Num_alive = starting_fish_count - DB_Mort_Gensum$Num_dead}

  #Generate rows of data representing survivors
  DB_Mort_Genalive = data.frame(lapply(DB_Mort_Gensum, rep, DB_Mort_Gensum$Num_alive))
  DB_Mort_Genalive$Status = 0
  DB_Mort_Genalive$TTE = last_tte
  DB_Mort_Gencomb = plyr::rbind.fill(mort_db, DB_Mort_Genalive[, -c(3:4)])

  #Create prism output
  if(output_prism == TRUE){

    prism_db = data.frame(blank = rep("", nrow(DB_Mort_Gencomb)))

    #Get treatment specific column entries for Status
    for(col_nm in levels(factor(DB_Mort_Gencomb$Trt.ID))) {
      temp_db = data.frame(ifelse(DB_Mort_Gencomb$Trt.ID == col_nm, DB_Mort_Gencomb$Status, ""))
      colnames(temp_db) = col_nm
      prism_db = cbind(prism_db, temp_db)
    }

    #Organize prism data and save
    if(is.null(output_prism_date)) {
      prism_db = cbind(data.frame(DB_Mort_Gencomb[, -which(colnames(DB_Mort_Gencomb) == "Status")]),
                       prism_db[, -1])
    } else {
      prism_db = cbind(data.frame(starting_date = format(as.Date(output_prism_date, format = "%d-%b-%Y"), "%d-%b-%Y"),
                                  ending_date = format(as.Date(output_prism_date, format = "%d-%b-%Y") + DB_Mort_Gencomb$TTE, "%d-%b-%Y")),
                       data.frame(DB_Mort_Gencomb[, -which(colnames(DB_Mort_Gencomb) == "Status")]),
                       prism_db[, -1])
    }

    prism_db = data.frame(prism_db %>% dplyr::arrange(Trt.ID, Tank.ID))[, -which(colnames(prism_db) == "Trt.ID")]
    write.csv(prism_db, paste("Surv_Gen Prism - last TTE ", last_tte, ".csv", sep = ""))
  }

  print(paste("Your total number of tanks is:", length(levels(factor(DB_Mort_Gencomb$Tank.ID)))))
  print(paste("Your total number of treatment groups is:", length(levels(factor(DB_Mort_Gencomb$Trt.ID)))))
  print(paste("Your total number of fish in the output data is:", nrow(DB_Mort_Gencomb)))
  return(DB_Mort_Gencomb)
}

################################################# Function 7 - Surv_Plots() #######################################################

#' Generate Survival Plots
#'
#' @description Produces a Kaplan-Meier Survival Plot and/or Hazard Time Plot from survival data. Each plot contains multiple curves for the different treatment groups. Plots saved automatically to working directory.
#'
#' @details The survival dataset should be a dataframe containing at least 4 different columns:
#' * "Trt.ID" = Labels for treatment groups in the study.
#' * "Tank.ID" = Labels for tanks in the study (each tank must have a unique label).
#' * "TTE" = Time to Event. Event depends on "Status".
#' * "Status" = Value indicating what happened at TTE. 1 for dead fish, 0 for survivors or those sampled and removed.
#'
#' Each row should represent one fish. For an example dataframe, execute \code{data(surv_db_ex)} and view.
#'
#' For details on the statistical methodology used by \code{bshazard::bshazard()}, refer to: \href{https://www.researchgate.net/publication/287338889_bshazard_A_Flexible_Tool_for_Nonparametric_Smoothing_of_the_Hazard_Function}{here}.
#'
#' General concept: h(t) the hazard function is considered in an count model with the number of deaths as the response variable. I.e, death_count(t) = h(t) * P(t) where P(t) is the number alive as a function of time and h(t) is modeled over time using basis splines. The basis spline curvature\bold{s} is assumed to have a normal distribution with mean 0 (a random effect). Based on this assumption, the author found that the variance of curvatures (i.e. smoothness) is equal to the over-dispersion (phi) of the death counts related (divided) by some smoothness parameter (lambda). Phi and lambda can be estimated from the data or specified by the user. Specification can be helpful in low sample size situations where overdispersion (phi) estimates have been found to be unreliable and clearly wrong (based on my understanding of realistic estimates and what was estimated in past data with adequate, large sample sizes).
#' @md
#'
#' @param surv_db A survival dataframe as described in \bold{Details}.
#' @param plot_prefix A string specifying the prefix for the filename of the saved plots.
#' @param xlim A vector specifying the plots x-axis lower and upper limits, respectively.
#' @param ylim A vector specifying the Survival Plot y-axis lower and upper limits, respectively.
#' @param xlab A string specifying the plot x-axis label.
#' @param lambda Smoothing value for the hazard curve. Higher lambda produces greater smoothing. Defaults to NULL where \code{bshazard::bshazard()} uses the provided survival data to estimate lambda; NULL specification is recommended for large sample size situations which usually occurs on our full-scale studies with many mortalities and tank-replication. At low sample sizes, the lambda estimate can be unreliable. Choosing a lambda of 10 (or anywhere between 1-100) probably produces the most accurate hazard curve for these situations. In place of choosing lambda, choosing \code{phi} is recommended; see below.
#' @param phi Dispersion parameter for the count model used in hazard curve estimation. Defaults to NULL where \code{bshazard()} uses the provided survival data to estimate phi; NULL specification is recommended for large sample size situations. At low sample sizes, the phi estimate can be unreliable. Choosing a phi value of 1 for low sample sizes is recommended. This value of 1 (or close) seems to be that estimated in past Tenaci data (QCATC997; phi ~ 0.8-1.4) where there are large sample sizes with tank-replication. The phi value of 1 indicates the set of counts (deaths) over time have a Poisson distribution, following the different hazard rates along the curve and are not overdispersed (phi > 1).
#' @param dailybin Whether to set time bins at daily (1 TTE) intervals. Refer to the \code{bshazard()} documentation for an understanding on the role of bins to hazard curve estimation. Please set to TRUE at low sample sizes and set to FALSE for large sample sizes (often with tank replication), although at large sample sizes either TRUE or FALSE produces similar results usually. Defaults to TRUE.
#' @param plot Which plot to output. Use "surv" for the Kaplan-Meier Survival Curve, "haz" for the Hazard Curve, or "both" for both. Defaults to "both".
#' @param colours Vector of color codes for the different treatment groups in the plot. Defaults to ggplot2 default palette.
#' @param theme A string specifying the graphics theme for the plots. Theme "ggplot2" and "prism" currently available. Defaults to "ggplot2".
#' @param trt_order Vector representing the order of treatment groups in the plots. Defaults to NULL where alphabetical order is used.
#' @param data_out Whether to print out the survival and/or hazard databases illustrated by the plots. Defaults to FALSE.
#' @param plot_dim Vector representing the dimensions (width, height) with which to save the plot in .tiff and .pptx.
#'
#' @seealso \href{https://sean4andrew.github.io/safuncs/reference/Surv_Plots.html}{Link} for executed \bold{Examples} which includes any figure outputs.
#'
#' @return Rturns a list containing the Kaplan-Meier Survival Curve and the Hazard Curve if {\code{plot = "both"}}. If only one plot is to be calculated and shown, set either \code{plot = "haz"} or \code{plot = "surv"}.
#'
#' If \code{data_out = TRUE}, returns dataframes associated with the survival plots.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' # Starting from an example mortality database, we first generate the complete survivor
#' # data using Surv_Gen()
#' data(mort_db_ex)
#' surv_dat = Surv_Gen(mort_db = mort_db_ex,
#'                     starting_fish_count = 100,
#'                     last_tte = 54)
#'
#' # Create plot by inputting surv_dat into Surv_Plots()!
#' Surv_Plots(surv_db = surv_dat,
#'            plot_prefix = "QCATC777",
#'            xlim = c(0, 54),
#'            ylim = c(0, 1),
#'            xlab = "TTE",
#'            plot = "both",
#'            dailybin = FALSE)
Surv_Plots = function(surv_db,
                      plot_prefix = "plot_prefix",
                      xlim = NULL,
                      ylim = c(0, 1),
                      xlab = "Days Post Challenge",
                      lambda = NULL,
                      phi = NULL,
                      dailybin = TRUE,
                      plot = "both",
                      colours = NULL,
                      theme = "ggplot",
                      trt_order = NULL,
                      data_out = FALSE,
                      plot_dim = c(6, 4)) {

  if(is.null(xlim)) {xlim <- c(0, max(surv_db$TTE))}
  if(!is.null(trt_order)){surv_db$Trt.ID = factor(surv_db$Trt.ID, levels = trt_order)}

  if(plot == "surv" | plot == "both") {
  surv_obj = survminer::surv_fit(survival::Surv(TTE, Status) ~ Trt.ID, data = surv_db)

    if(length(levels(as.factor(surv_db$Trt.ID))) > 1) {
      attributes(surv_obj$strata)$names <- levels(as.factor(surv_db$Trt.ID))
    } else {
      surv_obj$strata = length(surv_obj$surv)
      attributes(surv_obj$strata)$names <- levels(as.factor(surv_db$Trt.ID))
    }

  surv_dat = data.frame(Trt.ID = summary(surv_obj)$strata,
                        Survprob = summary(surv_obj)$surv,
                        Time = summary(surv_obj)$time)

  surv_plot = survminer::ggsurvplot(surv_obj,
                                    conf.int = FALSE,
                                    ggtheme = theme(plot.background = element_rect(fill = "white")),
                                    break.y.by = 0.1,
                                    break.x.by = max(round(max(xlim) / 13), 1),
                                    xlim = xlim,
                                    ylim = ylim,
                                    xlab = xlab,
                                    surv.scale = "percent")
  Survival_Plot = surv_plot$plot + theme(legend.position = "right") + guides(color = guide_legend("Trt.ID"))
  if(theme == "prism") {Survival_Plot = Survival_Plot + ggprism::theme_prism()}
  eoffice::topptx(figure = Survival_Plot, filename = paste(plot_prefix, "Survival-Curve.pptx", sep = "-"), width = plot_dim[1], height = plot_dim[2])

  if(!is.null(colours)) {Survival_Plot = Survival_Plot + scale_color_manual(values = colours)}
  ggsave(paste(plot_prefix, "Survival-Curve.tiff", sep = "-"), dpi = 300, width = plot_dim[1], height = plot_dim[2], plot = Survival_Plot)

  }

  if(dailybin == TRUE) {dbin <- max(surv_db$TTE)}
  if(dailybin == FALSE) {dbin <- NULL}

  #create Haz_list
  if(plot == "haz" | plot == "both") {
  Haz_list = list()
  for(Haz_Trt in levels(as.factor(surv_db$Trt.ID))) {
    surv_db_trt = surv_db[surv_db$Trt.ID == Haz_Trt,]
    if(sum(surv_db_trt$Status) == 0){
      Haz_list[[Haz_Trt]] = data.frame(Hazard = 0,
                                       Time = rep(0, max(surv_db$TTE), 1))
    } else {
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
  }

  haz_db = dplyr::bind_rows(Haz_list, .id = "Trt.ID")
  if(!is.null(trt_order)){haz_db$Trt.ID = factor(haz_db$Trt.ID, levels = trt_order)}
  Hazard_Plot = ggplot(data = haz_db, aes(x = Time, y = Hazard, color = Trt.ID)) +
    geom_line(linewidth = 1) +
    geom_point() +
    xlab(xlab) +
    scale_x_continuous(breaks = seq(from = min(xlim),
                                    to = max(xlim),
                                    by = max(round(max(xlim) / 13), 1)),
                       limits = xlim)

  if(!is.null(colours)) {Hazard_Plot = Hazard_Plot + scale_color_manual(values = colours)}
  if(theme == "prism") {Hazard_Plot = Hazard_Plot + ggprism::theme_prism()}
  ggsave(paste(plot_prefix, "Hazard-Curve.tiff", sep = "-"), dpi = 300, width = plot_dim[1], height = plot_dim[2], plot = Hazard_Plot)
  eoffice::topptx(figure = Hazard_Plot, filename = paste(plot_prefix, "Hazard-Curve.pptx", sep = "-"), width = plot_dim[1], height = plot_dim[2])
  }

  if(data_out == TRUE) {
    if(plot == "surv") {return(list(Survival_Plot = Survival_Plot, Survival_DB = surv_dat))}
    if(plot == "haz") {return(list(Hazard_Plot = Hazard_Plot, Hazard_DB = haz_db))}
    if(plot == "both") {return(list(Survival_Plot = Survival_Plot, Survival_DB = surv_dat, Hazard_Plot = Hazard_Plot, Hazard_DB = haz_db))}
  } else {
    if(plot == "surv") {return(Survival_Plot = Survival_Plot)}
    if(plot == "haz") {return(Hazard_Plot = Hazard_Plot)}
    if(plot == "both") {return(list(Survival_Plot = Survival_Plot, Hazard_Plot = Hazard_Plot))}
  }
}

################################################# Function 8 - GG_Colour_Hue() #######################################################

#' Get Default Colours by ggplot
#'
#' @description Not my function but it is useful so here it is! \href{https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette}{Origin}
#'
#' @param n Number of colour groups
#'
#' @return Returns a vector representing the default colour codes assigned to each group by ggplot.
#' @export
#'
#' @seealso \href{https://sean4andrew.github.io/safuncs/reference/GG_Colour_Hue.html}{Link} for executed \bold{Examples} which includes any figure outputs.
#'
#' @examples
#' # Get colour codes used for 6 categorical groups
#' GG_Colour_Hue(6)
#'
GG_Colour_Hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#################################################### Function 9 - Label_Gen() ####################################################

#' @title Generate Texts for Labels
#'
#' @description Combines texts specified in a list. List should contain multiple variables that holds the texts in a value or vector format. All combinations of texts across variables are computed. Output text (combinations) is sorted in order of variables given in the list (default behavior) or as specified using \code{sort_by} argument. The output combinations are tabulated and saved in a .csv in your working directory.
#'
#' @param input_list A list of named variables, each containing one or more text/number(s). See \bold{Examples} for examples.
#' @param sort_by A value or vector representing the variable(s) to sort the output by. For each variable, sorts according to the order of text in the variable. When multiple variables is given, prioritizes sorting based on the order of variables; leftmost = highest priority. Defaults to NULL where sorting is based on \code{input_list} orders.
#' @param n_col The number of columns in the output table. It should match the number of columns in the label paper. Defaults to 6.
#' @param fill_by_row Whether combinations should fill the output table by row (otherwise column). Defaults to TRUE.
#' @param save_name Name of the saved .csv. Defaults to NULL where the file name is "Label_Gen" and today's date (YYYY-MM-DD).
#'
#' @return A .csv containing all possible combinations. Additionally, a printout describing the .csv file name and location. Another printout describing the total number of labels / combinations created.
#' @export
#'
#' @seealso \href{https://sean4andrew.github.io/safuncs/reference/Label_Gen.html}{Link} for web documentation.
#'
#' @examples
#' # Summarize the input variables in a list
#' input_variables = list(Time = c("Baseline", "1wpv"),
#'                        Animal = c("Oysters", "Lobsters"),
#'                        Tissue = c("Meat", "Shell", "Water", "Head"),
#'                        Replic_num = 1:3)
#'
#' # Run Label_Gen() using the input variables.
#' Label_Gen(input_list = input_variables,
#'           sort_by = c("Time", "Animal", "Tissue"),
#'           n_col = 6,
#'           fill_by_row = TRUE,
#'           save_name = NULL)
Label_Gen = function(input_list,
                     sort_by = NULL,
                     n_col = 6,
                     fill_by_row = TRUE,
                     save_name = NULL) {

  # Create combination data frame
  poss_grid = expand.grid(input_list)

  # Sort output
  if(is.null(sort_by)) {
    sort_by = names(input_list)
  }
  rev_sb = rev(sort_by)
  for(i in rev_sb){
    poss_grid = poss_grid[order(poss_grid[, which(colnames(poss_grid) == i)]),]
  }

  # Store ordered combinations
  ordered_combos = interaction(poss_grid, sep = ", ")
  extended_combos = c(paste(ordered_combos),
                      rep("", times = ceiling(length(ordered_combos)/n_col) * n_col - length(ordered_combos)))
  mat_combos = matrix(extended_combos, ncol = n_col, byrow = fill_by_row)
  colnames(mat_combos) = 1:n_col

  # Save and print outputs
  print(paste("You have", length(ordered_combos), "total labels"))

  if(is.null(save_name)){
    write.csv(x = mat_combos, file = paste("Label_Gen ", Sys.Date(), ".csv", sep = ""))
    print(paste("File saved as", paste("Label_Gen ", Sys.Date(), ".csv", sep = ""), "in", getwd()))
  } else {
    write.csv(x = mat_combos, file = paste(save_name, ".csv", sep = ""))
    print(paste("File saved as", paste(save_name, ".csv", sep = ""), "in", getwd()))
  }
}

##################################################### Function 10 - Surv_Power() ####################################################

#' Calculate Power for Survival Studies
#'
#' @description Calculates the power of global and/or pairwise hypothesis tests for survival studies with support over a range of experimental designs. This versatility is enabled by the simulation-based approach of the power calculation, using the modular \code{Surv_Simul()} to simulate survival data of various experimental designs. Power calculations can be made to account for inter-tank variation using a mixed cox proportional hazards model (set argument \code{model = "coxph_glmm"}). Additionally, power calculations can account for the multiplicity of pairwise comparisons using \code{pairwise_corr}. Users can compare power across different experimental designs by specifying each as a list element in \code{Surv_Simul()}. The results are returned as dataframes and plots.
#'
#' @details Power calculation follows the standard procedure for simulation-based approaches. First, the user simulates hypothetical future sample sets using \code{Surv_Simul()}. For each sample set, a p-value is calculated by \code{Surv_Power()}. The percentage of p-values below 0.05 (positives) were then calculated, representing power. The percent positives can also represent false positive rate if the population/truth from which different treatments are simulated are identical.
#'
#' @param simul_db The simulated survival dataframe from \code{Surv_Simul()} with the desired experimental design parameters.
#' @param global_test A character vector representing the method(s) to use for global hypothesis testing of significance of treatment. Methods available are: "logrank", "wald", "score", "LRT". "logrank" represents the global logrank test of significance. The latter three methods are standard global hypothesis testing methods for models. They are only available when the argument \code{model} is specified (i.e. not NULL)."wald" represents the Wald Chisquare Test (also known as joint test) which assesses whether model parameters (log(hazard ratios)) jointly are significantly different from 0 (i.e. HRs ≠ 1). Wald test can be done for various cox-proportional hazard models that could be relevant to our studies (glm, glmm, and gee). Due to its broad applicability, while also producing practically the same p-value most of the time compared to the other model tests, "wald" is the recommended option of the three. "score" represents the Lagrange multiplier or Score test. 'LRT' represents the likelihood ratio test. Defaults to "logrank" for now due to its ubiquity of use.
#' @param model A character vector representing the model(s) to fit for hypothesis testing. Models available are: "coxph_glm" and "coxph_glmm". "coxph_glm" represents the standard cox proportional hazard model fitted using \code{survival::coxph()} with Trt.ID as a fixed factor. "coxph_glmm" represents the mixed cox proportional hazard model fitted using \code{coxme::coxme()} with Trt.ID as a fixed factor and Tank.ID as a random factor to account for inter-tank variation. Defaults to NULL where no model is fitted for hypothesis testing.
#' @param pairwise_test A character vector representing the method(s) used for pairwise hypothesis tests. Use "logrank" to calculate power for logrank tests comparing different treatments. Use "EMM" to calculate power using Estimated Marginal Means based on model estimates (from 'coxph_glm' and/or 'coxph_glmm'). Defaults to "logrank".
#' @param pairwise_corr A character vector representing the method(s) used to adjust p-values for multiplicity of pairwise comparisons. For clarification, this affects the power of the pairwise comparisons. Methods available are: "tukey", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", and "none". In \bold{Details}, I discuss categories of the adjustment methods and provided a recommendation for "BH". Defaults to "none" for now.
#' @param prog_notes Whether to print the number of sample sets with p-values calculated. Defaults to TRUE.
#' @param plot_out Whether to display plot(s) illustrating power from global and/or pairwise hypothesis tests. Defaults to TRUE
#' @param plot_lines Whether to plot lines connecting points of the same "group" in the plot output. Defaults to TRUE.
#' @param xlab A string representing the x-axis title. Defaults to "List Element #".
#' @param xnames Vector of names for x-axis labels. Defaults to NULL where names are the list element numbers from \code{Surv_Gen()}.
#' @param plot_save Whether to save plots as a .tiff. Defaults to TRUE.
#'
#' @return Output. An SE for proportions (calculated using the binomial formula)
#'
#' @import magrittr
#' @import ggplot2
#'
#' @export
#'
#' @seealso \href{https://sean4andrew.github.io/safuncs/reference/Surv_Power.html}{Link} for web documentation.
#'
#' @examples To be made..
Surv_Power = function(simul_db = simul_db_ex,
                      global_test = "logrank",
                      model = NULL,
                      pairwise_test = "logrank",
                      pairwise_corr = "none",
                      prog_notes = TRUE,
                      plot_out = TRUE,
                      plot_lines = FALSE,
                      xlab = "List Element #",
                      xnames = NULL,
                      plot_save = TRUE){

  #Track time elapsed
  time_start = Sys.time()

  #Convert NULL pairwise_corr into "none"
  if(is.null(pairwise_corr)){pairwise_corr <- "none"}

  #Standardize simul_db as dataframe
  if(!is.data.frame(simul_db)){simul_db = data.frame(simul_db$simul_surv_db)}

  #Add a value of 1 for column list_element_num in case no value is present in simul_db
  if(!"list_element_num" %in% colnames(simul_db)){simul_db$list_element_num <- 1}

  #Validation checks (stops)
  if(sum(!global_test %in% c("logrank", "wald", "score", "LRT")) > 0) {
    stop(paste("The", global_test[!global_test %in% c("logrank", "wald", "score", "LRT")][1], "global test method is not in the list supported by Surv_Power(). Select any amount from 'logrank', 'wald', 'score', and/or 'LRT'. For no global test to be done, select NULL. "))
  }

  if(sum(!model %in% c("coxph_glm", "coxph_glmm")) > 0) {
    stop(paste("The", model[!model %in% c("coxph_glm", "coxph_glmm")][1], "model is currently not in the list supported by Surv_Power(). Select any amount from 'coxph_glm' and/or 'coxph_glmm'. For no model to be fitted, select NULL."))
  }

  if(sum(!pairwise_test %in% c("logrank", "EMM")) > 0) {
    stop(paste("The", pairwise_test[!pairwise_test %in% c("logrank", "EMM")][1], "pairwise test method is not in the list supported by Surv_Power(). Select any amount from 'logrank' and/or 'EMM'. For no pairwise test to be done, select NULL."))
  }

  pairwise_corr_options = c("tukey", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none")
  if(sum(!pairwise_corr %in% pairwise_corr_options) > 0) {
    stop(paste("The", pairwise_corr[!pairwise_corr %in% pairwise_corr_options][1], "pairwise correction method is not in the list supported by Surv_Power(). Select any amount from 'tukey', 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', and/or 'none'."))
  }

  #Initialize objects to store values
  power_glob = data.frame()
  power_pair = data.frame()
  prog = 0

  #Subset data by list_element_num. Loop through each.
  for(ele_num in unique(simul_db$list_element_num)) {
    simul_db_temp0 = simul_db[simul_db$list_element_num == ele_num,] #filter for ele_num

    #Clear stored p_values for every ele_num
    p_pair = data.frame()
    p_glob = list()

    #Calculate a p-value for every loopnum
    for(simnum in unique(simul_db_temp0$n_sim)) {
      simul_db_temp = simul_db_temp0[simul_db_temp0$n_sim == simnum,] #filter for loopnum

      #Logrank tests
      #Global
      if("logrank" %in% global_test){
        p_glob[["N/Ap"]][["logrank"]][simnum] = survival::survdiff(survival::Surv(TTE, Status) ~ Trt.ID, simul_db_temp)$pvalue
      }

      #Pairwise
      for(pairwise_corr_id0 in pairwise_corr[pairwise_corr %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")]){
        if("logrank" %in% pairwise_test) {
          pair_lr_res = survminer::pairwise_survdiff(survival::Surv(TTE, Status) ~ Trt.ID,
                                                     simul_db_temp, p.adjust.method = pairwise_corr_id0)
          temp_pair2 = data.frame(as.table(pair_lr_res$p.value))
          temp_pair2 = temp_pair2[-which(is.na(temp_pair2$Freq)),]

          p_pair = rbind(p_pair, data.frame(pair = interaction(temp_pair2$Var2, temp_pair2$Var1, sep = " - "),
                                            pvalues = temp_pair2$Freq,
                                            model = "N/Ap",
                                            pairwise_test = "logrank",
                                            corr = pairwise_corr_id0))
        }
      }

      #Model fits
      if("coxph_glm" %in% model){
        coxph_glm = survival::coxph(survival::Surv(TTE, Status) ~ Trt.ID, simul_db_temp)
        coxph_glm_sum = summary(coxph_glm)
      }
      if("coxph_glmm" %in% model){coxph_glmm <- coxme::coxme(survival::Surv(TTE, Status) ~ Trt.ID + (1|Tank.ID), simul_db_temp)}

      #Repeat for every model
      for(mod_id in model){ #for every model...

        #Repeat for every pairwise comparison correction setting
        for(pairwise_corr_id in pairwise_corr) { #for every pairwise comparison setting..
          if("EMM" %in% pairwise_test) {
            temp_pair = data.frame(emmeans::emmeans(mget(mod_id, envir = environment())[[1]],
                                                    pairwise ~ Trt.ID, adjust = pairwise_corr_id)$contrasts)

            p_pair = rbind(p_pair, data.frame(pair = temp_pair$contrast,
                                              pvalues = temp_pair$p.value,
                                              model = mod_id,
                                              pairwise_test = "EMM",
                                              corr = pairwise_corr_id))
          }
        }

        #Repeat for every global_test setting
        for(glob_id in global_test){
          if(glob_id == "wald"){p_glob[[mod_id]][[glob_id]][simnum] <- emmeans::joint_tests(mget(mod_id,
                                                                                                 envir = environment())[[1]])$p.value}
          if(glob_id == "score"){
            if(mod_id == "coxph_glm"){p_glob[[mod_id]][[glob_id]][simnum] <- coxph_glm_sum$waldtest["pvalue"]}
            if(mod_id == "coxph_glmm"){p_glob[[mod_id]][[glob_id]][simnum] <- NA} #Method not available/allowed
          }

          if(glob_id == "LRT"){
            if(mod_id == "coxph_glm"){p_glob[[mod_id]][[glob_id]][simnum] <- coxph_glm_sum$logtest["pvalue"]}
            if(mod_id == "coxph_glmm"){
              p_glob[[mod_id]][[glob_id]][simnum] =
                anova(coxph_glmm, coxme::coxme(survival::Surv(TTE, Status) ~ 1 + (1|Tank.ID), simul_db_temp))$`P(>|Chi|)`[2]
            }
          }
        }
      }

      #Print progress
      if(prog_notes == TRUE) {cat("\rCalculated p-values for", prog <- prog + 1, "of",
                                  max(simul_db$list_element_num) * max(simul_db_temp0$n_sim), "sample sets")}
    } #Close loop for simnum

    #Create power tables from p-values for each ele_num
    #For global test pvalues
    if(length(p_glob) > 0){
      p_glob_unlist = stack(unlist(p_glob))
      p_glob_unlist$ind = gsub(pattern = "[0-9]", x = p_glob_unlist$ind, replacement = "")
      p_glob_db = data.frame(tidyr::separate(data = p_glob_unlist, col = "ind",
                                             into = c("model", "global_test"), sep = "\\."))
      power_glob_temp = data.frame(p_glob_db %>%
                                     dplyr::group_by(model, global_test) %>%
                                     dplyr::summarise(percent_signif = 100 * sum(values < 0.05)/length(values),
                                                      datasets_n = length(values), .groups = "drop"))
      power_glob_temp$percent_signif_se = sqrt(power_glob_temp$percent_signif *
                                                 (100-power_glob_temp$percent_signif) /
                                                 (power_glob_temp$datasets_n))
      power_glob_temp = power_glob_temp[, c("model", "global_test", "percent_signif", "percent_signif_se", "datasets_n")]
      power_glob_temp$element_num = ele_num
      power_glob = rbind(power_glob, power_glob_temp)
    }

    #For pairwise test pvalues
    if(length(p_pair) > 0){
      power_pair_temp = data.frame(p_pair %>%
                                     dplyr::group_by(pair, model, pairwise_test, corr) %>%
                                     dplyr::summarise(percent_signif = 100 * sum(pvalues < 0.05)/length(pvalues),
                                                      datasets_n = length(pvalues), .groups = "drop"))
      power_pair_temp$percent_signif_se = sqrt(power_pair_temp$percent_signif *
                                                 (100-power_pair_temp$percent_signif) /
                                                 (power_pair_temp$datasets_n))
      power_pair_temp = power_pair_temp[, c("pair", "model", "pairwise_test", "corr", "percent_signif",
                                            "percent_signif_se", "datasets_n")]
      power_pair_temp$element_num = ele_num
      power_pair = rbind(power_pair, power_pair_temp)
    }
  } #Close loop for ele_num

  #Outermost steps
  #Plot #1 (global test)
  if(length(p_glob) > 0){
    mod_col = c("N/Ap" = "#F8766D", "coxph_glm" = "#00BA38", "coxph_glmm" = "#00BFC4")
    power_glob$model = factor(power_glob$model, levels = c("N/Ap", "coxph_glm", "coxph_glmm"))
    power_glob$global_test = factor(power_glob$global_test, levels = c("logrank", "wald", "score", "LRT"))
    test_shapes = c(15:18)
    names(test_shapes) = c("logrank", "wald", "score", "LRT")

    glob_plot = ggplot(data = na.omit(power_glob), aes(x = as.numeric(element_num), y = percent_signif/100,
                                                       colour = model, group = interaction(model, global_test))) +
      geom_errorbar(aes(ymin = (percent_signif - percent_signif_se)/100,
                        ymax = (percent_signif + percent_signif_se)/100),
                    position = position_dodge(width = 0.12), width = 0.1) +
      geom_point(aes(shape = global_test), position = position_dodge(width = 0.12)) +
      scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.1), limits = c(0, 1),
                         name = "% of significant results (p<0.05)") +
      scale_x_continuous(breaks = seq(1, max(power_glob$element_num), 1),
                         name = ifelse(is.null(xlab), "List Element #", xlab),
                         labels = if(is.null(xnames)){waiver()}
                         else{stringr::str_wrap(xnames, width = round(24/max(power_pair$element_num)))}) +
      labs(color = "Model", shape = "Test", title = "Global Test of Significance") +
      theme(plot.title = element_text(hjust = 0)) +
      guides(linetype = "none",
             shape = guide_legend(order = 2),
             colour = guide_legend(order = 1)) +
      scale_color_manual(values = mod_col) +
      scale_shape_manual(values = test_shapes)

    if(plot_lines == TRUE) {glob_plot <- glob_plot + geom_line(aes(linetype = global_test), position = position_dodge(width = 0.12))}

  } else {glob_plot <- NULL}

  #Plot #2 (pairwise test)
  if(length(p_pair) > 0){
    power_pair$pair = gsub(" - ", " vs. ", power_pair$pair)
    u_pairs = length(unique(power_pair$pair))
    n_col = ceiling(u_pairs/2)
    test_and_corr_combos = rev(levels(interaction(c("logrank", "EMM"), pairwise_corr_options, sep = " & ")))[-1]
    power_pair$test_and_corr = interaction(power_pair$pairwise_test, power_pair$corr, sep = " & ")
    power_pair$test_and_corr = factor(power_pair$test_and_corr, levels = test_and_corr_combos)
    power_pair$model = factor(power_pair$model, levels = c("N/Ap", "coxph_glm", "coxph_glmm", "coxph_gee"))
    test_and_corr_shapes = c(3, rep(c(15:18, 4, 7, 25), each = 2))
    names(test_and_corr_shapes) = test_and_corr_combos

    pair_plot = ggplot(data = power_pair, aes(x = as.numeric(element_num), y = percent_signif/100, colour = model,
                                              group = interaction(model, test_and_corr))) +
      facet_wrap(~pair, ncol = n_col) +
      geom_errorbar(aes(ymin = (percent_signif - percent_signif_se)/100,
                        ymax = (percent_signif + percent_signif_se)/100),
                    position = position_dodge(width = 0.30), width = 0.1) +
      geom_point(aes(shape = test_and_corr), position = position_dodge(width = 0.30)) +
      scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.1), limits = c(0, 1),
                         name = "% of significant results (p<0.05)") +
      scale_x_continuous(breaks = seq(1, max(power_pair$element_num), 1),
                         name = ifelse(is.null(xlab), "List Element #", xlab),
                         labels = if(is.null(xnames)){waiver()}
                         else{stringr::str_wrap(xnames, width = round(48/(max(power_pair$element_num) * n_col)))}) +
      labs(color = "Model", shape = "Test & Correction", title = "Pairwise Test of Significance") +
      theme(plot.title = element_text(hjust = 0)) +
      guides(linetype = "none",
             shape = guide_legend(order = 2),
             colour = guide_legend(order = 1)) +
      scale_color_manual(values = mod_col) +
      scale_shape_manual(values = test_and_corr_shapes)


    if(plot_lines == TRUE) {pair_plot <- pair_plot + geom_line(aes(linetype = test_and_corr), position = position_dodge(width = 0.30))}
  } else {pair_plot <- NULL}

  #Return output
  output = list()
  if(!is.null(global_test)){
    output[["power_global_db"]] = power_glob
    if(plot_out == TRUE & length(power_glob) > 0) {
      output[["power_global_plot"]] = glob_plot
      if(plot_save == TRUE) {
        ggsave(paste("Power_Global_Test_", Sys.Date(), ".tiff", sep =""),
               dpi = 900, width = 7, height = 5.5, plot = glob_plot)
      }
    }
  }
  if(!is.null(pairwise_test)){
    output[["power_pairwise_db"]] = power_pair[, -which((colnames(power_pair) == "test_and_corr"))]
    if(plot_out == TRUE & length(pair_plot) > 0) {
      output[["power_pairwise_plot"]] = pair_plot
      if(plot_save == TRUE) {
        ggsave(paste("Power_Pairwise_Test_", Sys.Date(), ".tiff", sep =""),
               dpi = 900, width = 1.5 + (2.75 * ceiling(u_pairs / n_col)),
               height = 2.75 * n_col, plot = pair_plot)
      }
    }
  }

  #Validation checks (print notes)
  ##For global tests
  if("coxph_glmm" %in% model & "score" %in% global_test){
    print("NOTE: No off-the-shelf function for conducting a Likelihood Ratio Test on 'coxph_glmm' models in R. No power value is returned for the global test using LRT on 'coxph_glmm' model.")
  }

  ##For pairwise tests
  if("logrank" %in% pairwise_test & "tukey" %in% pairwise_corr) {
    print("NOTE: Tukey pairwise correction is not available for log-rank tests. No power value is returned for such a combination of test and correction.")
  }

  if(is.null(model)){
    if("EMM" %in% pairwise_test){
      print("NOTE: No pairwise comparison of type 'EMM' was done since no model was specified.")
    }
    if(sum(global_test %in% c("wald", "score", "LRT")) > 0){
      print(paste("NOTE: No global_test of type", global_test, "was done since no model was specified."))
    }
  }

  if(is.null(pairwise_test)){
    print("NOTE: No pairwise test was done since pairwise_test was set to NULL.")
  }
  if(is.null(global_test)){
    print("NOTE: No global test was done since global_test was set to NULL.")
  }

  #Print time elapsed
  print(paste("Time elapsed:", substr(hms::as_hms(Sys.time() - time_start), 1, 8), "(hh:mm:ss)"))
  return(output)
}

##################################################### Data 1 - mort_db_ex #######################################################

#' Example Mort Data
#'
#' @description A subset of columns taken from the online excel mortality file in OneDrive.
#' @usage
#' data(mort_db_ex)
#' view(mort_db_ex)
#'
#' @format A data frame containing 399 rows and 4 columns:\tabular{lll}{
#'  \code{Tank.ID} \tab \tab Unique labels for the different tanks in the study \cr
#'  \code{Trt.ID} \tab \tab Unique labels for the different treatments in the study \cr
#'  \code{TTE} \tab \tab Time to Event. In this dataset, TTE = days post challenge \cr
#'  \code{Status} \tab \tab Value indicating what happened at TTE. In this dataset, Status = 1 indicating all events are death \cr
#' }
#'
"mort_db_ex"

##################################################### Data 2 - surv_db_ex #######################################################

#' Example Survival Data
#'
#' @description A complete survival dataset with survivors, based of \code{mort_db_ex}. Ready for proper survival analysis.
#' @usage
#' data(surv_db_ex)
#' view(surv_db_ex)
#'
#' @format A data frame containing 1200 rows and 4 columns:\tabular{lll}{
#'  \code{Tank.ID} \tab \tab Unique labels for the different tanks in the study \cr
#'  \code{Trt.ID} \tab \tab Unique labels for the different treatments in the study \cr
#'  \code{TTE} \tab \tab Time to Event. In this dataset, TTE = days post challenge \cr
#'  \code{Status} \tab \tab Value indicating what happened at TTE. In this dataset, Status = 1 or 0 indicating death or survival, respectively \cr
#' }
#'
"surv_db_ex"

##################################################### Data 3 - haz_db_ex #######################################################

#' Example Hazard Data
#'
#' @description A reference hazard dataframe created using \code{Surv_Plots(data_out = TRUE)$Hazard_DB} which uses \code{bshazard::bshazard()}. Contains hazard rates over time.
#' @usage
#' data(haz_db_ex)
#' view(haz_db_ex)
#'
#' @format A data frame containing 54 rows and 3 columns:\tabular{lll}{
#'  \code{Trt.ID} \tab \tab A label for the treatment group used in creating this reference hazard dataframe \cr
#'  \code{Hazard} \tab \tab Hazard values (rates) \cr
#'  \code{Time} \tab \tab Time / TTE in days post challenge. \cr
#' }
#'
"haz_db_ex"
