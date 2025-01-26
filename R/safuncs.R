# This is an R package containing useful functions for my work.

# Available functions with documentation:
# 1. Con_Simul() -- simulates contingency tables based on the multinomial distribution.
# 1b. Con_Simul_PR() -- calculates positive rates for statistical tests on contingency tables.
# 3. Surv_Simul() -- simulate survival data based on a reference hazard function, the specified hazard ratio(s), and inter-tank variation.
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
#' @param prog_show Whether to display the progress of \code{Surv_Simul()} by printing the number of simulations completed. Defaults to TRUE.
#' @param n_sim Number of survival dataset to simulate. Defaults to 1. For serious power calculations, an n_sim < 2000 is likely the most you'll ever need for precise results.
#' @param plot_out Whether to output the information plot (further details in \bold{Value}). Defaults to TRUE.
#' @param pop_out Whether to output a dataframe containing the survival probability values for the population. Defaults to TRUE.
#' @param theme A string specifying the graphics theme for the plots. Theme "ggplot2" and "prism" currently available. Defaults to "ggplot2".
#' @param plot_save Whether to save plot as .tiff in the working directory. Defaults to TRUE.
#'
#' @return Returns a list that, at minimum, contains the simulated survival dataframe (\code{surv_simul_db}) which has at least 5 columns: TTE (Time to Event), Status (0 / 1), Trt.ID, Tank.ID, and n_sim which represents the simulation number for the data subsets. Additionally, can contain a column named "list_element_num" which represents the list element number when an input argument to \code{Surv_Simul()} is specified as a list.
#'
#' If \code{plot_out = TRUE}, the list includes Kaplan-Meier survival plots. The left faceted plot represents the survival curves for the simulated sample set, while the right represents that for the population/truth. Numbers at the end of survival curves represent the end survival rate for that treatment. If the number of simulated sample sets exceed 1, multiple survival curves are drawn with each representing a sample set. In such cases, a statement is also provided informing of the probability to detect the effect of Treatment using a global logrank test from \code{survival::survdiff()}, i.e. the power or false positive rate.
#'
#' If \code{pop_out = TRUE}, the list includes a dataframe (\code{surv_pop_db}) representing the survival rates (probabilities) for the population / truth from which the sample is simulated. The dataframe contains three columns: Trt.ID, surv_prob which represents the survival probabilities, and TTE.
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
#'            prog_show = FALSE, #hide simulation progress notes for cleaner output
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
#'            prog_show = FALSE,
#'            n_sim = 4)$surv_plots
#'
#' # What if we want to compare power of the global log-rank test (shown in the plot)
#' # across different experimental setups with different fish numbers per treatment?
#' # Below, I setup a Surv_Simul() to answer this question.
#' Surv_Simul(haz_db = haz_db_ex,
#'            fish_num_per_tank = list(30, 100),
#'            tank_num_per_trt = 3,
#'            treatments_hr = c(1, 0.6),
#'            prog_show = FALSE,
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
                      prog_show = TRUE,
                      plot_out = TRUE,
                      pop_out = TRUE,
                      theme = "ggplot2",
                      plot_save = TRUE) {

  #Track time elapsed
  time_start = Sys.time()

  #Making sure input data has correct (lower case) column names
  colnames(haz_db) = tolower(colnames(haz_db))

  #Validation checks
  if(length(levels(factor(haz_db$tr.id))) > 1) {stop("Please use only one Trt.ID as reference in the supplied hazard dataframe.")}
  if(n_sim >= 10000) {
    print("NOTE: It is not recommended to use n_sim > 10000 as this may create impractically large simulated data files (>1 GB).")
  }
  if(plot_out == TRUE & n_sim > 100) {
    print("Plot too slow to render with n_sim > 50. No plot output is provided")
    plot_out = FALSE
  }

  #Add blank for last time point in haz_db, for convenience on later calculations
  haz_db = rbind(haz_db, data.frame(trt.id = haz_db$trt.id[1], hazard = last(cumsum(haz_db$hazard)) * .Machine$double.eps * 10,
                                    time = round(last(haz_db$time))))
  #Initialize objects to store second output type (across list elements and loops)
  output2 = list(surv_plots = list(), surv_simul_db = data.frame(), surv_pop_db = data.frame())
  list_var_check = c()

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

    #Establish treatment names
    Trt_names = c("A (Control)", LETTERS[2:(length(treatments_hr))])

    #if you have list elements, assign and print. If not just jump straight to old code
    if(length(list_var) > 1) { #if you have list elements, assign and print, otherwise just go to old code.
      assign(var_name, list_var[[ele_num]]) #assign
    }

    #Initialize objects to store loop results
    surv_samps = data.frame() #for plotting purposes
    cens_db = data.frame() #for plotting purposes
    pvalues = c() #for plotting purposes
    Surv_simul_outDB = as.list(seq_len(n_sim))

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
            CDF_Yval = c(CDF_Yval, CDF_Yval_temp)

            Trt.ID = c(Trt.ID, rep(Trt_names[iTT], length(CDF_Yval_temp)))
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
            CDF_Yval = c(CDF_Yval, CDF_Yval_temp)

            Trt.ID = c(Trt.ID, rep(Trt_names[iTT], length(CDF_Yval_temp)))
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
                                      Trt.ID = rep(unique(Trt_levels), times = nrow(sampling_specs)))
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
      pvalues = c(pvalues, survival::survdiff(survival::Surv(TTE, Status) ~ Trt.ID, Surv_simul_DB)$pvalue)

      #Simulated survival data to be provided as output
      if(length(list_var) > 1){Surv_simul_DB$list_element_num <- ele_num}
      Surv_simul_outDB[[loopnum]] = Surv_simul_DB

      #Transform simulated survival data for plotting purposes
      if(plot_out == TRUE) {
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
                                                      alpha = max(c(0.1, 1 - (0.0001 ^ (1/n_sim))))))
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
      }

      #Print progress
      if(prog_show != FALSE) {cat("\rSimulated", prog <- prog + 1, "of", n_sim * length(list_var), "sample sets")}
    } #close loopnum
    Surv_simul_outDB = do.call(rbind, Surv_simul_outDB)

    #Get "population" survival dataset by exponentiating the negative cumulative hazard
    pop_haz_db = data.frame(approx(x = haz_db$time, y = haz_db$hazard, xout = seq(min(haz_db$time), max(haz_db$time), 0.1), method = "linear"))
    colnames(pop_haz_db) = c("time", "hazard")

    Trt_names = c("A (Control)", LETTERS[2:(length(treatments_hr))])
    surv_pop = data.frame(Trt.ID = as.factor(rep(Trt_names, each = length(haz_db$hazard))),
                          surv_prob = exp(-as.vector(apply(haz_db$hazard %*% t(treatments_hr), 2, cumsum))),
                          time = rep(haz_db$time, times = length(treatments_hr)),
                          type = "Population / truth",
                          n_sim = 1,
                          alpha = 1)
    if(length(list_var) > 1){surv_pop$list_element_num <- ele_num}

    #To the end of creating survival plots
    if(plot_out == TRUE) {
      surv_comb = rbind(surv_samps, surv_pop)
      surv_comb$type = factor(surv_comb$type, levels = c(paste("Sample set (n = ", n_sim, ")", sep = ""), "Population / truth"))
      surv_comb$Trt.ID = factor(surv_comb$Trt.ID, levels = rep(Trt_names))

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
        coord_cartesian(clip = "off", xlim = c(min(haz_db$time), max(haz_db$time) - 0.5)) +
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
      if(length(list_var) > 1) {

        surv_plots = surv_plots + labs(title = paste("List Element", ele_num))
      }

      #Save plot
      if(plot_save == TRUE){
        ggsave(paste("Surv_Simul_Plot",
                     ifelse(length(list_var) == 1, "_", paste("_Element", ele_num, "_", sep ="")),
                     Sys.Date(), ".tiff", sep = ""), dpi = 900, width = 7, height = 4, plot = surv_plots)
      }
    }
    #remove columns "alpha", "type", and "n_sim" from data output
    surv_pop = surv_pop[, -c(4:6)]
    colnames(surv_pop) = c("Trt.ID", "surv_prob", "TTE")

    #Return R output if list_var length = 1 (i.e. no list)
    Surv_simul_outDB = dplyr::arrange(Surv_simul_outDB, n_sim, Tank.ID, Trt.ID, TTE)
    if(length(list_var) == 1) {
      if(plot_out == FALSE & pop_out == FALSE) {
        #Print time elapsed
        print(paste("Time elapsed:", substr(hms::as_hms(Sys.time() - time_start), 1, 8), "(hh:mm:ss)"))
        return(surv_simul_db = Surv_simul_outDB)

      } else {

        output = list(surv_simul_db = Surv_simul_outDB)

        if(pop_out == TRUE) {output$surv_pop_db <- surv_pop}
        if(plot_out == TRUE) {output$surv_plots <- surv_plots}

        #Print time elapsed
        print(paste("Time elapsed:", substr(hms::as_hms(Sys.time() - time_start), 1, 8), "(hh:mm:ss)"))
        return(output)
      }
    }

    if(length(list_var) > 1) {

      #Store 2nd output if list_var length >1
      if(plot_out == TRUE) {output2$surv_plots[[ele_num]] <- surv_plots}
      output2$surv_simul_db = rbind(output2$surv_simul_db, Surv_simul_outDB)
      output2$surv_pop_db = rbind(output2$surv_pop_db, surv_pop)
    }
  } #This closes the loop that deals with lists

  if(length(list_var) > 1){

    if(plot_out == FALSE) {
      output2$surv_plots = NULL
    }

    if(pop_out == FALSE) {
      output2$surv_pop_db = NULL
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
#' @return Returns a list containing the Kaplan-Meier Survival Curve and the Hazard Curve if {\code{plot = "both"}}. If only one plot is to be calculated and shown, set either \code{plot = "haz"} or \code{plot = "surv"}.
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

#' Calculates Power for Survival Experiments
#'
#' @description Calculates the power of global and/or pairwise hypothesis tests for survival studies with support over a range of experimental designs. This versatility is enabled by the simulation-based approach of the power calculation, using the modular \code{Surv_Simul()} to simulate survival data of various experimental designs. Power calculations can be made to account for inter-tank variation using a mixed cox proportional hazards model (set argument \code{model = "coxph_glmm"}). Additionally, power calculations can account for the multiplicity of pairwise comparisons using \code{pairwise_corr}. Users can compare power across different experimental designs by specifying each as a list element in \code{Surv_Simul()}. A brief tutorial is written in \bold{Examples} to guide the user on how to use \code{Surv_Power()} to calculate power under various scenarios.
#'
#'
#' @param simul_db An output from \code{Surv_Simul()} which includes the survival dataframe simulated with the desired experimental design parameters.
#' @param global_test A character vector representing the method(s) to use for global hypothesis testing of significance of treatment. Methods available are: "logrank", "wald", "score", "LRT". "logrank" represents the global logrank test of significance. The latter three methods are standard global hypothesis testing methods for models. They are only available when the argument \code{model} is specified (i.e. not NULL). "wald" represents the Wald Chisquare Test which assesses whether model parameters (log(hazard ratios)) jointly are significantly different from 0 (i.e. HRs ≠ 1). Wald test can be done for various cox-proportional hazard models that could be relevant to our studies (glm, glmm, and gee). Due to its broad applicability, while also producing practically the same p-value most of the time compared to the other model tests, "wald" is the recommended option of the three. "score" represents the Lagrange multiplier or Score test. 'LRT' represents the likelihood ratio test. Defaults to "logrank" for now due to its ubiquity of use.
#' @param model A character vector representing the model(s) to fit for hypothesis testing. Models available are: "coxph_glm" and "coxph_glmm" ("cox_gee" may be supported upon request, omitted for reasons not discussed here). "coxph_glm" represents the standard cox proportional hazard model fitted using \code{survival::coxph()} with Trt.ID as a fixed factor. "coxph_glmm" represents the mixed cox proportional hazard model fitted using \code{coxme::coxme()} with Trt.ID as a fixed factor and Tank.ID as a random factor to account for inter-tank variation. Defaults to NULL where no model is fitted for hypothesis testing.
#' @param pairwise_test A character vector representing the method(s) used for pairwise hypothesis tests. Use "logrank" to calculate power for logrank tests comparing different treatments. Use "EMM" to calculate power using Estimated Marginal Means based on model estimates (from 'coxph_glm' and/or 'coxph_glmm'). Defaults to "logrank".
#' @param pairwise_corr A character vector representing the method(s) used to adjust p-values for multiplicity of pairwise comparisons. For clarification, this affects the power of the pairwise comparisons. Methods available are: "tukey", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", and "none". Under \bold{Details} (yet to be finished), I discuss the common categories of adjustment methods and provided a recommendation for "BH". Defaults to c("none", "BH").
#' @param prog_show Whether to display the progress of \code{Surv_Power()} by printing the number of sample sets with p-values calculated. Defaults to TRUE.
#' @param data_out Whether to output dataframes containing the power of global and/or pairwise hypothesis tests. Defaults to TRUE.
#' @param plot_out Whether to output plots illustrating the power of global and/or pairwise hypothesis tests. Defaults to TRUE.
#' @param plot_lines Whether to plot lines connecting points of the same hypothesis test in the plot output. Defaults to TRUE.
#' @param xlab A string representing the x-axis title. Defaults to "List Element #".
#' @param xnames A character vector of names for x-axis labels. Defaults to NULL where names are the list element numbers from \code{Surv_Simul()}.
#' @param plot_save Whether to save plots as a .tiff in the working directory. Defaults to TRUE.
#'
#' @return Outputs a list containing any of the following four items depending on input arguments:
#'
#' \itemize{
#'  \item When \code{global_test ≠ NULL} and \code{data_out = TRUE}, outputs a dataframe named \code{power_glob_db} containing power values calculated for global hypothesis tests. The dataframe consists of six columns: \tabular{lll}{
#'  \code{model} \tab \tab The type of model being evaluated in power calculations \cr
#'  \code{global_test} \tab \tab The global hypothesis test being evaluated \cr
#'  \code{power} \tab \tab The percentage of p-values below 0.05, i.e. power \cr
#'  \code{power_se} \tab \tab Standard error for percentages \cr
#'  \code{sample_sets_n} \tab \tab Number of sample sets used in calculating power \cr
#'  \code{list_element_num} \tab \tab The list element number associated with the power value calculated \cr
#'  }
#'  \item When \code{global_test ≠ NULL} and \code{plot_out = TRUE}, a plot showing power values for global hypothesis test. Plot corresponds to \code{power_glob_db}.
#'  \item When \code{pairwise_test ≠ NULL} and \code{data_out = TRUE}, outputs a dataframe named \code{power_pair_db} containing power values for pairwise hypothesis tests. The dataframe consists of eight columns: \tabular{lll}{
#'  \code{pair} \tab \tab The treatment groups to be compared \cr
#'  \code{model} \tab \tab The type of model being evaluated in power calculations \cr
#'  \code{pairwise_test} \tab \tab The type of pairwise_test being evaluated \cr
#'  \code{pairwise_corr} \tab \tab The method of p-value correction/adjustments for paiwrise comparisons \cr
#'  \code{power} \tab \tab The percentage of p-values below 0.05, i.e. power \cr
#'  \code{power_se} \tab \tab Standard error for percentages \cr
#'  \code{sample_sets_n} \tab \tab Number of sample sets used in calculating power \cr
#'  \code{list_element_num} \tab \tab The list element number associated with the power value calculated \cr
#'  }
#'  \item When \code{pairwise_test ≠ NULL} and \code{plot_out = TRUE}, a plot showing power values across list elements. Plot corresponds to \code{power_pair_db}.
#' }
#'
#' @import magrittr
#' @import ggplot2
#'
#' @export
#'
#' @seealso \href{https://sean4andrew.github.io/safuncs/reference/Surv_Power.html}{Link} for web documentation.
#'
#' @examples
#' # Below is a tutorial on how to calculate power using Surv_Power(). The function
#' # calculates power as a percentage of positive (p < 0.05) test conducted on simulated
#' # datasets. Hence, the first step is simulating the data
#'
#' # A past data is retrieved for simulating future datasets:
#' data(haz_db_ex)
#' haz_db_ex = haz_db_ex[haz_db_ex$Trt.ID == "A",] # filter for control fish
#' head(haz_db_ex, n = 5)
#'
#' # As may be clear from the above, we are using the past data's hazard curve properties.
#' # NOTE: Whether the shape of the hazard curve repeats in the future study only matters
#' # to the accuracy of the power calculation when the future study involves fish being
#' # dropped (e.g. due to sampling).
#'
#' # To begin simulating, input the past data to Surv_Simul() with the supposed future
#' # experiment parameters:
#' surv_sim_db_ex1 = Surv_Simul(haz_db = haz_db_ex,
#'                             fish_num_per_tank = 100,
#'                             tank_num_per_trt = 4,
#'                             treatments_hr = c(1, 0.7, 0.5),
#'                             logHR_sd_intertank = 0,
#'                             n_sim = 500,
#'                             prog_show = FALSE, #omit progress bar for cleaner output
#'                             plot_out = FALSE) #omit plotting for efficiency/speed
#'
#' # Above, we simulated 500 datasets with each having 400 fish per treatment and three
#' # treatments total with hazard ratios of 1, 0.7, and 0.5 relative to 'haz_db_ex.'
#'
#' # Next, the simulated data can be supplied to Surv_Power() to calculate power. Below,
#' # power is calculated for the logrank test:
#' Surv_Power(simul_db = surv_sim_db_ex1,
#'            global_test = "logrank",
#'            pairwise_test = "logrank",
#'            pairwise_corr = "none",
#'            prog_show = FALSE,
#'            data_out = FALSE) # remove data output for brevity
#'
#' # From the above, we can see that there is a high probability of detecting significance
#' # of treatment using the global test. On the other hand, power to detect differences
#' # between treatment B (HR = 0.7 relative to control) and C (HR = 0.5) is ~ 60%.
#' # Notably, the presented power for pairwise comparisons uses unadjusted p-values.
#'
#' # If there is a need to adjust p-values to provide a guarantee for, for example, the
#' # false discovery rate, then supply the chosen FDR method(s) (e.g. "BH") to the
#' # 'pairwise_corr' argument:
#' Surv_Power(simul_db = surv_sim_db_ex1,
#'            global_test = "logrank",
#'            pairwise_test = "logrank",
#'            pairwise_corr = c("none", "BH"),
#'            prog_show = FALSE,
#'            data_out = FALSE)
#'
#' # By default, FDR adjustment methods "control" or limit the false discovery rate to 5%;
#' # true FDR would be lower.
#'
#' # If for some reason there is a need to control for family-wise error rate due to a
#' # desire to conclude at the family-level (already achievable via global test), then p
#' # values can be adjusted using FWER methods (e.g. "bonferroni", "tukey"). Notably, due
#' # to the correlation between pairwise tests, the tukey method is possible and is more
#' # powerful than bonferroni. The tukey option is available if a model is specified:
#' Surv_Power(simul_db = surv_sim_db_ex1,
#'            model = "coxph_glm",
#'            global_test = c("logrank", "wald"),
#'            pairwise_test = c("logrank", "EMM"),
#'            pairwise_corr = c("none", "tukey", "bonferroni"),
#'            prog_show = FALSE,
#'            data_out = FALSE)
#'
#' # Based on the pairwise plot above, it seems that the tukey method is only marginally
#' # more powerful in this case compared to bonferroni. The model used above was the cox
#' # proportional hazards model. The global test option "wald" corresponds to a wald
#' # chi-square test of the significance of treatment as a model factor. The pairwise
#' # test option "EMM" corresponds to treatment comparisons using model estimated
#' # marginal means. Notably, the model methods appear to produced similar result to the
#' # logrank test, at least in this example.
#'
#' # More sophisticated models may be fitted to account for tank variation if any was
#' # introduced in the simulation process as below:
#' surv_sim_db_ex2 = Surv_Simul(haz_db = haz_db_ex,
#'                              fish_num_per_tank = 100,
#'                              tank_num_per_trt = 6,
#'                              treatments_hr = c(1, 1, 0.7),
#'                              logHR_sd_intertank = 0.2,
#'                              n_sim = 500,
#'                              prog_show = FALSE,
#'                              plot_out = FALSE)
#'
#' # To calculate power considering the tank variation, we fit a mixed model (coxph_glmm):
#' Surv_Power(simul_db = surv_sim_db_ex2,
#'            model = "coxph_glmm",
#'            global_test = c("logrank", "wald"),
#'            pairwise_test = c("logrank", "EMM"),
#'            pairwise_corr = c("none"),
#'            prog_show = FALSE,
#'            data_out = FALSE)
#'
#' # The above results show a property of the test that accounted for tank variation;
#' # it is weaker than the one that does not (logrank test) given the data has tank
#' # variability. However, this is for a reason. If you look into the positive rate of
#' # the test comparing Control vs Trt. B (where true HR = 1), it only slightly above 5%
#' # for the mixed model (a real bias for reasons not discussed here). This result is
#' # roughly consistent with the acclaimed alpha or false positive rate of 0.05. In
#' # contrast, the false positive rate for logrank test is much higher at ~20%.
#'
#' ## VARYING EXPERIMENTAL CONDITIONS
#' # As a precautionary measure, it may be important to study power under various scenarios
#' # (e.g. increasing or decreasing strength of pathogen challenge). To do this, specify
#' # each scenario as separate elements of a list in an argument of Surv_Simul():
#' HR_vec = c(1, 0.7, 0.5)
#' HR_list = list(1.5 * HR_vec, # strong challenge
#'                1.0 * HR_vec, # medium challenge
#'                0.5 * HR_vec) # weak challenge
#'
#' surv_sim_db_ex3 = Surv_Simul(haz_db = haz_db_ex,
#'                              fish_num_per_tank = 100,
#'                              tank_num_per_trt = 4,
#'                              treatments_hr = HR_list,
#'                              logHR_sd_intertank = 0,
#'                              n_sim = 500,
#'                              prog_show = FALSE,
#'                              plot_out = FALSE)
#'
#' # Next, calculate and compare power across scenarios:
#' Surv_Power(simul_db = surv_sim_db_ex3,
#'            global_test = "logrank",
#'            pairwise_test = "logrank",
#'            pairwise_corr = "none",
#'            prog_show = FALSE,
#'            data_out = FALSE,
#'            xlab = "Challenge Strength",
#'            xnames = c("Strong", "Medium", "Weak"))
#'
#' # The results show that power for the strong challenge model is generally high.
#' # However, for the weak challenge model, more sample size appears to be needed to
#' # achieve 80% power for some pairwise comparisons. We can investigate at what sample
#' # size would such power be achieved by specifying different sample sizes as separate
#' # list elements in Surv_Simul():
#' fish_num_vec = seq(from = 50, to = 200, by = 30)
#' fish_num_vec
#'
#' fish_num_list = as.list(fish_num_vec) # convert vector elements to list elements
#'
#' # Simulate the different experimental conditions:
#' surv_sim_db_ex4 = Surv_Simul(haz_db = haz_db_ex,
#'                              fish_num_per_tank = fish_num_list,
#'                              tank_num_per_trt = 4,
#'                              treatments_hr = c(1, 0.7, 0.5) * 0.5,
#'                              logHR_sd_intertank = 0,
#'                              n_sim = 500,
#'                              prog_show = FALSE,
#'                              plot_out = FALSE)
#'
#' # Compare power across sample sizes:
#' Surv_Power(simul_db = surv_sim_db_ex4,
#'            global_test = "logrank",
#'            pairwise_test = "logrank",
#'            pairwise_corr = "none",
#'            prog_show = FALSE,
#'            data_out = FALSE,
#'            xlab = "Fish number per tank",
#'            xnames = fish_num_vec,
#'            plot_lines = TRUE)
#'
#' # The results showed that even an increase in sample size (from 100 to 200) result in
#' # only modest gains in power for comparisons B to C. Perhaps a possible future action
#' # then is to ensure the challenge is medium to strong by having higher salinity for
#' # example. The above examples show some possible use case of Surv_Power(), but other
#' # cases can be simulated to understand power in various scenarios. I hope this tool
#' # can help you calculate power for your specific needs!
Surv_Power = function(simul_db = simul_db_ex,
                      global_test = "logrank",
                      model = NULL,
                      pairwise_test = "logrank",
                      pairwise_corr = c("none", "BH"),
                      prog_show = TRUE,
                      data_out = TRUE,
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
  if(!is.data.frame(simul_db)){simul_db = data.frame(simul_db$surv_simul_db)}

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
          temp_pair2 = na.omit(data.frame(as.table(pair_lr_res$p.value)))

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

        #Repeat for every global_test setting (maybe dont need this loop, just assess for inclusioon in the iefs)
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
      if(prog_show == TRUE) {cat("\rCalculated p-values for", prog <- prog + 1, "of",
                                 max(simul_db$list_element_num) * max(simul_db$n_sim), "sample sets")}
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
                         name = "Power (%)") +
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
                         name = "Power (%)") +
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

  if(length(power_pair) == 8) {
    colnames(power_pair) = c("pair", "model", "pairwise_test", "pairwise_corr", "power", "power_se", "sample_sets_n", "list_element_num")
  }

  #Return output
  output = list()
  if(!is.null(global_test)){
    #Rename columns
    colnames(power_glob) = c("model", "global_test", "power", "power_se", "sample_sets_n", "list_element_num")
    #Create output
    if(data_out == TRUE) {output[["power_global_db"]] <- power_glob}
    if(plot_out == TRUE & length(power_glob) > 0) {
      output[["power_global_plot"]] = glob_plot
      if(plot_save == TRUE) {
        ggsave(paste("Power_Global_Test_", Sys.Date(), ".tiff", sep =""),
               dpi = 900, width = 7, height = 5.5, plot = glob_plot)
      }
    }
  }
  if(!is.null(pairwise_test)){
    power_pair = power_pair[, -which((colnames(power_pair) == "test_and_corr"))]
    #Rename columns
    colnames(power_pair) = c("pair", "model", "pairwise_test", "pairwise_corr", "power", "power_se", "sample_sets_n", "list_element_num")
    #Create output
    if(data_out == TRUE) {output[["power_pairwise_db"]] <- power_pair}
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

##################################################### Function 11 - MultiVar() ####################################################

#' @title Analyze and Visualize Multivariate Data
#'
#' @description Expedite multivariate analysis and gain comprehensive insights to the data through PCA, LDA, and MANOVA. Outputs include tables of statistical results, PCA plots, and Boxplots of each dependent variable. Plots can be customized using \code{pca_} and \code{boxplots_} prefixed arguments. Supports one- or two-factor analyses. In two-factor analyses, additional plots may be created with facets and/or pooling of values across levels of a selected factor, chosen using the arguments \code{factors_pool} and \code{factors_facet}. By default, saves results in a Word document, but allows exports as .pptx, .png, and/or R objects for further edits. A tutorial on how to generate various outputs from \code{MultiVar()} is available under \bold{Examples}.
#'
#' @details Several functions published on CRAN are used by \code{MultiVar()} for various types of statistical analyses:
#'
#' \describe{
#'  \item{PCA}{Uses \code{stats::prcomp()}. Pre-PCA, outcome variables may be scaled or centered using \code{scale} and \code{center} arguments, respectively. Missing values may have been imputed using \code{missMDA::imputePCA()} and ncp parameter \code{ncp = missMDA::estim_ncpPCA()$ncp}.}
#'  \item{LDA}{Uses \code{MASS::lda()}. Pre-LDA, values may be scaled, centered, and missing values imputed, with the same methods as described for PCA.}
#'  \item{MANOVA}{Uses \code{stats::manova()} and subsequently \code{car::Anova()} with the argument \code{type = 3} - relevant for analyses with two or more factors. Missing values may have been imputed with the same method described for PCA. By default, \code{stats::manova()} omits all rows with missing values. To achieve this in \code{MultiVar()}, set \code{missing_method} argument to "na_omit".}
#'  \item{PERMANOVA}{Output available upon request. Uses \code{RVAideMemoire::adonis.II()} with the argument \code{method = "euclidean"}. Missing values may be imputed using the same method described for PCA. By default, \code{RVAideMemoire::adonis.II()} omits all rows with missing values. To achieve this in \code{MultiVar()}, set \code{missing_method} argument to "na_omit".}
#'  \item{ANOVA}{Uses \code{stats::anova()} and subsequently \code{car::Anova()} with the argument \code{type = 3} - relevant for analyses with two or more factors.}
#'  \item{Pairwise Comparisons}{Uses \code{emmeans::emmeans()}. P-values are corrected using the Benjamini Hochberg procedure to control for a false Discovery Rate of 5 percent.}
#' }
#'
#' @param multivar_db A dataframe with two types of columns. The first holds numeric values of the multivariate outcome variables, with each column containing one variable. The second type holds categories (factor levels), with each column containing one factor. An example dataframe can be viewed by running \code{View(data(multivar_db_ex))} in the R console.
#' @param values_cols A numeric vector specifying the order number of columns containing the outcome variables.
#' @param factors_cols A numeric vector specifying the order number of columns containing the factors. Maximum of two numbers (i.e. factors).
#' @param factors_pool A character vector indicating the factors which levels are to be pooled across in additional plots. Choose any combination of "col1" and "col2" which refers to the first and second column in \code{factors_cols}. Defaults to c("col1", "col2").
#' @param factors_facet A character vector indicating the factors which levels are to be faceted across in additional plots. Choose any combination of "col1" and "col2" which refers to the first and second column in \code{factors_cols}. Defaults to "none" which creates no additional plots.
#' @param pca_ellipse A character vector representing the type of ellipses to draw in PCA plots. Generates a plot for every specified type. Choose any combination of "confidence", "distribution", "convexhull", and/or "none". "confidence" draws ellipses representing the 95 percent confidence interval about the center of multivariate normal data (principal component scores); drawn using \code{ggpubr::stat_conf_ellipse()}. "distribution" represents ellipses expected to cover 95 percent of all multivariate normal data; drawn using \code{ggplot2::stat_ellipse()} with argument \code{type = "norm"}. "convexhull" represents the smallest convex polygon enclosing all points; drawn using \code{ggpubr::stat_chull()}. For plots without ellipses, include "none". Defaults to c("confidence").
#' @param pca_labels A character vector representing the labels to draw in PCA plots. Choose any combination of "ind" and/or "var". "ind" represents individual point labels by their row number. "var" represents variable loadings drawn as arrows; the arrow length and direction are calculated as in \code{factoextra::fviz_pca_biplot()}. Defaults to NULL (no labels drawn).
#' @param pca_shapes TRUE/FALSE indicating whether to use different shapes for factor levels in PCA plots. Defaults to FALSE. Different shapes are not provided for plots with greater than 6 factor levels.
#' @param missing_method A string representing the method to address missing values in \code{values_cols}. Choose from "imputation" or "na_omit". "imputation" fills in missing values with values created (imputed) based on the correlation between variables essentially; accomplished using \code{missMDA::imputePCA()} with the ncp parameter \code{missMDA::estim_ncpPCA()}. "na_omit" removes entire rows of data when at least one NA value is present. This method may result in a significant loss of data. Defaults to "imputation". The choice of \code{missing_method} would affect PCA, LDA, and MANOVA results but likely only to a small degree with few missing values. Has no impact on boxplots and ANOVAs.
#' @param scale TRUE/FALSE indicating whether to scale variable values (such that SD = 1 for each variable) before PCA or LDA. A common procedure in the z-score normalization of values that commonly precede PCA. It is not recommended to set this to FALSE, unless justified. Defaults to TRUE.
#' @param center TRUE/FALSE indicating whether to center variable values (such that mean = 0 for each variable) before PCA or LDA. A common procedure in the z-score normalization of values that commonly precede PCA. It is not recommended to set this to FALSE, unless justified. Defaults to TRUE.
#' @param boxplot_filled TRUE/FALSE indicating whether to color the insides of boxplots and points (i.e. fill them). If FALSE, boxplots and points are hollow with colored borders using \code{colour} in \code{ggplot2::aes()}, instead of using \code{fill}. Defaults to TRUE.
#' @param boxplot_x_angle A number describing the degree of tilt in the x-axis labels of the boxplots. Defaults to NULL (horizontal labels).
#' @param boxplot_x_wrap The maximum number of characters on a single line that would be split if a space bar is available between them. Defaults to NULL (no text wrapping).
#' @param boxplot_x_lab TRUE/FALSE indicating whether to include a title for the x-axis of the boxplots. Defaults to FALSE.
#' @param boxplot_legend TRUE/FALSE indicating whether to include the legend for the different groups in the boxplots. Defaults to TRUE.
#' @param boxplot_var_sep TRUE/FALSE indicating whether to include boxplots made separately for every variable. Defaults to FALSE.
#' @param colours A named character vector specifying the colors to use for different factor levels. E.g. For a factor with levels "A", "B", and "C", the \code{colours} vector may look like \code{c('A' = "brown", B = 'blue', C = '#f8e723')}. Defaults to NULL (default ggplot2 colours).
#' @param plot_out_png TRUE/FALSE indicating whether to save plots as .png in the working directory. Defaults to FALSE.
#' @param plot_out_pptx TRUE/FALSE indicating whether to save plots as editable forms in .pptx in the working directory. Defaults to FALSE.
#' @param plot_out_R TRUE/FALSE indicating whether to output plots as ggplot2 objects in a list in R. Defaults to FALSE.
#'
#' @return Returns a Word document containing at least the following types of results (in the given order):
#'
#' \strong{Multivariate Analyses:}
#'
#' \enumerate{
#' \item PC plot(s) illustrating the separation (if any) of PC scores between groups
#' \item Correlation plot between individual variables and principal components
#' \item Correlation table summarizing pearson correlation coefficients between all pairs of variables and principal components
#' \item Contribution table summarizing contribution values of variables to principal components
#' \item Contribution table summarizing contribution values of variables to major linear discriminants (i.e. to the separation between groups)
#' \item MANOVA table summarizing statistical evidence for any effect of factor(s)
#' }
#'
#' \strong{(Bonus) Univariate Analyses:}
#'
#' \enumerate{
#' \item Boxplot(s) comparing values of individual variables across groups
#' \item ANOVA results table (Two-Factor or One-Factor)
#' \item Pairwise comparisons table
#' }
#'
#' @import ggplot2
#' @import magrittr
#'
#' @export
#'
#' @seealso \href{https://sean4andrew.github.io/safuncs/reference/MultiVar.html}{Link} for web documentation.
#'
#' @examples
#' # Below is a brief tutorial to help you get various outputs from MultiVar()!
#'
#' # Lets start simple, with the aim of creating a report for a one-factor dataset
#' # called iris:
#' data(iris)
#'
#' # Data structure:
#' head(iris, n = 5)
#'
#' # We can see that columns #1 to #4 contain dependent variables, while column #5 holds
#' # the factor. Simply feed this information into Multivar()!
#' MultiVar(multivar_db = iris,
#'          values_cols = 1:4,
#'          factors_cols = 5)
#'
#' # A report is now saved in your working directory, check it out! You can locate your
#' # working directory using:
#' getwd()
#'
#' # Suppose you wanted a more equipped PCA plot. We can customize our plots using a
#' # set of arguments indicated with the proper prefix; 'pca_' for PCA plots and
#' # 'boxplots_' for Boxplots:
#' MultiVar(multivar_db = iris,
#'          values_cols = 1:4,
#'          factors_cols = 5,
#'          pca_ellipse = c("none", "confidence", "distribution", "convexhull"),
#'          pca_labels = c("var"),
#'          plot_out_R = TRUE)$pca$distribution
#'
#' # If you want to customize the plot beyond the capabilities of MultiVar(), save them
#' # as editable forms in .pptx or manipulate them further in R:
#' Sepal_L_box = MultiVar(multivar_db = iris,
#'                        values_cols = 1:4,
#'                        factors_cols = 5,
#'                        boxplot_x_angle = 45,
#'                        plot_out_pptx = TRUE,
#'                        plot_out_R = TRUE)$box$none$'Sepal Length'
#'
#' Sepal_L_box + ggthemes::theme_tufte() # Use the thufte theme
#'
#' ggplot2::ggsave(filename = "box_SepalLength.tiff",
#'                 width = 5,
#'                 height = 5,
#'                 dpi = 600) # Save with desired plot dimensions and resolution
#'
#' # MultiVar() can also handle missing values which are common in real life data. By
#' # default it imputes (creates) missing values based on the correlational properties
#' # between dependent variables (essentially). Alternatively, rows with missing values
#' # can be omitted using missing_method = "na_omit". Below is a demo showing the use
#' # of the imputation method.
#' iris[1:25, 1] = NA # replace half of Setosa's Sepal.Length with NA
#'
#' MultiVar(multivar_db = iris,
#'          values_cols = 1:4,
#'          factors_cols = 5,
#'          missing_method = "imputation",
#'          pca_ellipse = "distribution",
#'          plot_out_R = TRUE)$pca$distribution
#'
#' # The plot is not clearly different than before which indicates the imputed values
#' # substitute NAs effectively, without clear bias.
#'
#' ## TWO FACTOR CASES
#' # Real data often have two factors, and MultiVar() produces plots tailored to that.
#' # At minimum, MultiVar() produces plots with the interaction between two factors as
#' # x-axis labels or legend labels. This can lead to a lot of levels in the axis or
#' # legends which look unsightly. Run the code below to see how this looks like! Below
#' # I am using an example fish mucus dataset with one real factor (Treatment) and one
#' # fake (Fruits).
#' data(multivar_db_ex)
#'
#' head(multivar_db_ex, n = 5)
#'
#' MultiVar(multivar_db = multivar_db_ex,
#'          values_cols = 2:8,
#'          factors_cols = c(1, 9),
#'          factors_pool = "none",
#'          factors_facet = "none")
#'
#' # Check out the report in your working directory!
#'
#' # Additional plots can be created with facets or pooling of values across levels of a
#' # factor. Below we explore the default parameters which consists of 'factors_pool =
#' # c('col1', 'col2')' and its output:
#' MultiVar(multivar_db = multivar_db_ex,
#'          values_cols = 2:8,
#'          factors_cols = c(1, 9),
#'          factors_pool = c("col1", "col2"),
#'          plot_out_R = TRUE)$pca$confidence[c("pooled1", "pooled2")]
#'
#' # The "col1" specification instructs MultiVar() to create a plot with pooling across
#' # levels of the first factor specified in 'values_cols', while "col2" refers to the
#' # second factor.
#'
#' # Plots with facets can also be created. Below, we facet against levels of the second
#' # factor (Fruits) which is considered as a nuisance factor in this case:
#' MultiVar(multivar_db = multivar_db_ex,
#'          values_cols = 2:8,
#'          factors_cols = c(1, 9),
#'          factors_pool = c("col1", "col2"),
#'          factors_facet = c("col2"),
#'          plot_out_R = TRUE)$pca$confidence["facet2"]
#'
#' # Check your word document output and you will see that it contains facets and pooling
#' # for both PCA and boxplots. End of tutorial. Hope it helps!
MultiVar = function(multivar_db,
                    values_cols,
                    factors_cols,
                    factors_pool = c("col1", "col2"),
                    factors_facet = "none",
                    pca_ellipse = c("confidence"),
                    pca_labels = NULL,
                    pca_shapes = FALSE,
                    missing_method = "imputation",
                    scale = TRUE,
                    center = TRUE,
                    boxplot_filled = TRUE,
                    boxplot_x_angle = NULL,
                    boxplot_x_wrap = NULL,
                    boxplot_x_lab = FALSE,
                    boxplot_legend = TRUE,
                    boxplot_var_sep = FALSE,
                    colours = NULL,
                    plot_out_png = FALSE,
                    plot_out_pptx = FALSE,
                    plot_out_R = FALSE) {

  # Initialize and define objects
  results_doc = officer::read_docx() %>%
    officer::body_add_par(value = "Analyzing the Outcome Variables as a Multivariate Response", style = "heading 1") %>%
    officer::body_add_par("PCA Plots", style = "heading 2") %>%
    officer::body_add_par("") %>%
    officer::body_add_par("")
  results_doc_boxplot = officer::read_docx() %>%
    officer::body_add_par(value = "Analyzing the Outcome Variables Separately", style = "heading 1") %>%
    officer::body_add_par("All Variables in One Plot", style = "heading 2") %>%
    officer::body_add_par("") %>%
    officer::body_add_par("")
  tempf = tempfile(fileext = ".docx")
  pca_list = list(list())
  boxplot_list = list(list())
  dataframe_list = list()
  boxplot_name_vec = c()
  box_height2_vec = c()
  univariate_tests = TRUE # defined in case it is to be converted to an optional output
  knife <- function(x) { # for formatting p-values (mostly)
    if(x < 1e-4) {
      return(format(x, scientific = TRUE, digits = 4))
    } else {
      return(round(x, 4))
    }}

  # Set path and file names
  if(plot_out_png == TRUE) {
    img_path = file.path(getwd(), paste("MultiVar PNG Figures", format(Sys.Date(), "%d%b%Y")))
    suppressWarnings(dir.create(img_path))
  } else {
    img_path = file.path(tempdir(), "Temp MultiVar PNG Figures")
    suppressWarnings(dir.create(img_path))
  }

  if(plot_out_pptx == TRUE) {
    pptx_name = paste(sep = "", "MultiVar Editable Figures ",
                      format(Sys.Date(), "%d%b%Y"), ".pptx")
    if(file.exists(pptx_name)) {
      file.remove(pptx_name)
    }
  }

  # Get matrix of values
  matrix_values_ori = multivar_db[, values_cols]
  matrix_values = multivar_db[, values_cols]

  # Validation check(s)
  if(length(values_cols) < 2){stop("At least two columns are needed with variable values for proper multivariate analysis. Please add more columns to the 'values_cols' argument.")}
  if(length(factors_cols) > 2){stop("A maximum of two factors are allowed. Please reduce the number of columns specified in 'factors_cols' to ≤ 2.")}
  if(length(factors_cols) != 2) {factors_pool <- "none"}
  if(length(factors_cols) != 2) {factors_facet <- "none"}

  # Address NAs in matrix values
  if(missing_method == "imputation") {
    if(sum(is.na(matrix_values)) > 0) {
      matrix_values = data.frame(missMDA::imputePCA(matrix_values, ncp = missMDA::estim_ncpPCA(matrix_values)$ncp)$completeObs)
    }
  }
  if(missing_method == "na_omit") {
    multivar_db_ori = multivar_db
    multivar_db = multivar_db[stats::complete.cases(matrix_values),]
    matrix_values = multivar_db[, values_cols]
  }

  # Get PC values
  pc_values = stats::prcomp(matrix_values, scale = scale, center = center)
  matrix_values$PC1 = PC1 <- pc_values$x[, 1]
  matrix_values$PC2 = PC2 <- pc_values$x[, 2]

  # Get loadings
  load_db = data.frame(load1 = load1 <- pc_values$rotation[, 1] * pc_values$sdev[1],
                       load2 = load2 <- pc_values$rotation[, 2] * pc_values$sdev[2],
                       vars = stringr::str_wrap(width = 15, gsub("\\.", " ", rownames(pc_values$rotation))))

  # Create base database for plotting
  plot_db = cbind(multivar_db[, factors_cols, drop = FALSE], matrix_values)

  # Create base PCA plot
  pca_plot_base = ggplot(data = plot_db) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_minimal() +
    xlab(paste("Principal Component 1 (", round(summary(pc_values)$importance[2, 1] * 100, 1), "%)", sep = "")) +
    ylab(paste("Principal Component 2 (", round(summary(pc_values)$importance[2, 2] * 100, 1), "%)", sep = "")) +
    theme(axis.line = element_line(color = "black", linewidth = 0.4),
          strip.background = element_rect(fill = "grey", color = NA))

  # Create base boxplot
  boxplot = ggplot() +
    geom_boxplot(size = 0.8, outlier.size = -1, na.rm = TRUE) +
    geom_jitter(width = 0.2, height = 0, shape = 21, na.rm = TRUE) +
    theme_minimal() +
    theme(axis.line = element_line(color = "black", linewidth = 0.4),
          strip.background = element_rect(fill = "grey", colour = NA),
          axis.text.x = element_text(hjust = 0.5))
  if(!is.null(boxplot_x_angle)) {boxplot <- boxplot + theme(axis.text.x = element_text(angle = boxplot_x_angle, hjust = 1))}
  if(!is.null(boxplot_x_wrap)) {boxplot <- boxplot + scale_x_discrete(labels = scales::wrap_format(boxplot_x_wrap))}

  # Colour plots if applicable
  if(!is.null(colours)){
    pca_plot_base = pca_plot_base +
      scale_color_manual(values = colours) +
      scale_fill_manual(values = colours)

    boxplot = boxplot +
      scale_color_manual(values = colours) +
      scale_fill_manual(values = colours)
  }

  # Set default plot dimensions
  pca_width = 6.3
  pca_height = 4.45
  varplot_width = 6.4
  varplot_height = 4.5

  # Create factor conditions
  full_conds = c("none", "pooled1", "pooled2", "facet1", "facet2")
  treatment_conds = as.vector(factor(unique(c("none", gsub("col", "pooled", factors_pool),
                                              gsub("col", "facet", factors_facet))), levels = full_conds))

  # Create multiple plots (potentially) to deal with multi factor scenarios
  for(treatment_conds_i in treatment_conds) {

    # Create different databases for plotting in different factor scenarios
    if(treatment_conds_i == "none") {
      plot_db$base_factor = interaction(multivar_db[, factors_cols], sep = " ")
      base_factor_name = ifelse(length(factors_cols) == 1,
                                colnames(multivar_db)[factors_cols], "Combination")
    }
    if(treatment_conds_i %in% c("facet2", "pooled2")) {
      plot_db$base_factor = multivar_db[, factors_cols[1]]
      plot_db$second_factor = multivar_db[, factors_cols[2]]
      base_factor_name = gsub("\\.", " ", colnames(multivar_db)[factors_cols[1]])
      second_factor_name = gsub("\\.", " ", colnames(multivar_db)[factors_cols[2]])
    }
    if(treatment_conds_i %in% c("facet1", "pooled1")) {
      plot_db$base_factor = multivar_db[, factors_cols[2]]
      plot_db$second_factor = multivar_db[, factors_cols[1]]
      base_factor_name = gsub("\\.", " ", colnames(multivar_db)[factors_cols[2]])
      second_factor_name = gsub("\\.", " ", colnames(multivar_db)[factors_cols[1]])
    }

    # Create facet dimensions depending on treatment_conds_i
    facet_num = ifelse(treatment_conds_i %in% c("none", "pooled2", "pooled1"),
                       length(values_cols),
                       length(unique(plot_db$second_factor)))
    facet_col = ifelse(facet_num %in% 1:3, 2, 3)
    facet_row = ceiling(facet_num / facet_col)

    # Modify PCA Plot based on treatment_conds_i
    pca_plot = pca_plot_base
    pca_plot$data = plot_db
    pca_plot$mapping = aes(x = PC1, y = PC2, colour = base_factor, fill = base_factor)
    pca_plot$guides = guides(colour = guide_legend(base_factor_name), fill = guide_legend(base_factor_name),
                             shape = guide_legend(base_factor_name))
    if(pca_shapes == TRUE) {
      if(length(unique(plot_db$base_factor)) <= 6) {pca_plot$layers[[1]] <- geom_point(aes(shape = base_factor), size = 2)}
      if(length(unique(plot_db$base_factor)) > 6) {pca_plot$layers[[1]] <- geom_point(size = 2)}
    } else {pca_plot$layers[[1]] <- geom_point(size = 2)}

    # Add PCA labels
    if("ind" %in% pca_labels) {
      pca_plot = pca_plot + ggrepel::geom_text_repel(show.legend = FALSE, max.overlaps = 20, aes(label = rownames(plot_db)))
    }
    if("var" %in% pca_labels) {
      # Calculate arrow length multiplier as in factoextra::fviz_pca_biplot()
      r0.7 = min((max(PC1) - min(PC1)/(max(load1) - min(load1))), (max(PC2) - min(PC2)/(max(load2) - min(load2)))) * 0.7
      pca_plot$layers = append(pca_plot$layers, geom_segment(data = load_db,
                                                             aes(fill = NULL, color = NULL, x = 0, y = 0,
                                                                 xend = load1 * r0.7, yend = load2 * r0.7),
                                                             color = "black", alpha = 0.8,
                                                             arrow = arrow(length = unit(0.2, "cm"), type = "open"),
                                                             show.legend = FALSE), after = 0)
      text_size_algo = ifelse(treatment_conds_i %in% c("none", "pooled2", "pooled1"),
                              3, 3.4 - (0.6 * facet_col))
      pca_plot = pca_plot + ggrepel::geom_text_repel(show.legend = FALSE, max.overlaps = 20,
                                                     data = load_db, aes(fill = NULL, color = NULL,
                                                                         x = load1 * r0.7, y = load2 * r0.7,
                                                                         label = vars),
                                                     color = "black", size = text_size_algo)
    }

    # Add PCA facets
    if(treatment_conds_i %in% c("facet2", "facet1")) {
      pca_plot = pca_plot + facet_wrap(~ second_factor, ncol = facet_col, nrow = facet_row)

      # Set plot height based on facet
      pca_height = 1.8 + 1.2 * facet_row
    }

    # Add ellipses and store PCA plots in a list
    if("none" %in% pca_ellipse) {pca_list[["none"]][[treatment_conds_i]] <- pca_plot}
    if("confidence" %in% pca_ellipse) {pca_list[["confidence"]][[treatment_conds_i]] <- pca_plot + ggpubr::stat_conf_ellipse(geom = "polygon", alpha = 0.15)}
    if("distribution" %in% pca_ellipse) {pca_list[["distribution"]][[treatment_conds_i]] <- pca_plot + ggplot2::stat_ellipse(type = "norm", geom = "polygon", alpha = 0.15)}
    if("convexhull" %in% pca_ellipse) {pca_list[["convexhull"]][[treatment_conds_i]] <- pca_plot + ggpubr::stat_chull(geom = "polygon", alpha = 0.15)}

    # Save PCA plots
    for(ellipse_i in pca_ellipse) {
      pca_name = paste("PCA", ifelse(treatment_conds_i == "none", "",
                                     paste(" Facet-", second_factor_name, sep = "")),
                       " Ellipse-", tools::toTitleCase(ellipse_i), sep = "")
      if(treatment_conds_i %in% c("pooled2", "pooled1")) {
        pca_name = paste("PCA Pooled Across-", second_factor_name, " Ellipse-", tools::toTitleCase(ellipse_i), sep = "")
      }

      if(plot_out_pptx){
        eoffice::topptx(figure = pca_list[[ellipse_i]][[treatment_conds_i]],
                        filename = pptx_name,
                        width = pca_width, height = pca_height, append = TRUE)
      }
      ggsave(plot = pca_list[[ellipse_i]][[treatment_conds_i]], filename = paste(sep = "", pca_name, ".png"),
             path = img_path, dpi = 600, width = pca_width, height = pca_height)
      results_doc = officer::body_add_img(x = results_doc, sr = file.path(img_path, paste(pca_name, ".png", sep = "")),
                                          width = pca_width, height = pca_height) %>%
        officer::body_add_par(value = pca_name, style = "graphic title") %>%
        officer::body_add_par("")
      if((last(treatment_conds) == treatment_conds_i) & (last(pca_ellipse) == ellipse_i)){} else {
        results_doc <- results_doc %>% officer::body_add_par("")
      }
    }

    # Create variables plot, contributions table & variable correlations plot
    if(last(treatment_conds) == treatment_conds_i) {
      # Note for PCA
      results_doc = results_doc %>%
        officer::body_add_par(paste(sep = "", "The ", round(summary(pc_values)$importance[2, 1] * 100, 1), "% and ",
                                    round(summary(pc_values)$importance[2, 2] * 100, 1), "% values in the axes of the figures above describe the percent variation in the data explained by the first and second principal components, respectively."), style = "Normal") %>%
        officer::body_add_par("") %>%
        officer::body_add_par("")

      # Create variable plots
      varplot_name = "Correlation of Variables to PCs"
      varplot = factoextra::fviz_pca_var(pc_values, repel = TRUE) +
        xlab(expression(paste(sep = "", "Correlation (", italic(r), ") to PC1"))) +
        ylab(expression(paste(sep = "", "Correlation (", italic(r), ") to PC2"))) +
        theme(plot.title = element_blank(),
              axis.line = element_line(color = "black", linewidth = 0.4),
              plot.margin = unit(c(0, 0.18, 0, 0.18), "in"))
      varplot_layer3 = varplot$layers[[3]]
      varplot$layers[[3]] = ggrepel::geom_text_repel(data = varplot$data,
                                                     aes(x = x, y = y,
                                                         label = stringr::str_wrap(gsub("\\.", " ", name), width = 15)),
                                                     nudge_x = 0, nudge_y = 0, size = 3.75, color = "darkblue",
                                                     box.padding = 0, min.segment.length = 0.3, segment.colour = "darkblue")
      varplot$layers[[1]] = varplot_layer3

      # Save plot
      if(plot_out_pptx == TRUE) {
        eoffice::topptx(figure = varplot, filename = pptx_name,
                        width = varplot_width, height = varplot_height, append = TRUE)
      }
      ggsave(plot = varplot, filename = paste(sep = "", varplot_name, ".png"), path = img_path,
             dpi = 600, width = varplot_width, height = varplot_height)

      # Create correlations table
      correl_db = data.frame(cor(matrix_values)) %>% tibble::rownames_to_column(" ")
      colnames(correl_db) = gsub("\\.", " ", colnames(correl_db))
      correl_db[, 1] = gsub("\\.", " ", correl_db[, 1])

      correl_tab = flextable::flextable(correl_db) %>%
        flextable::set_caption(caption = "Table 1. Pearson correlation coefficient (r) between variables and principal components") %>%
        flextable::colformat_double(digits = ifelse(length(values_cols) < 5, 4, 3)) %>%
        flextable::line_spacing(space = ifelse(length(values_cols) < 5, 1, 0.8), part = "all") %>%
        flextable::fontsize(size = ifelse(length(values_cols) < 5, 11, 9), part = "all") %>%
        flextable::set_table_properties(layout = "autofit", width = 1)

      # Create contribution table
      contrib_db = data.frame(Variable = gsub("\\.", " ", names(pc_values$rotation[,1])),
                              ContributionPC1 = as.vector(100 * pc_values$rotation[,1]^2),
                              EigenvectorPC1 = as.vector(pc_values$rotation[,1]),
                              ContributionPC2 = as.vector(100 * pc_values$rotation[,2]^2),
                              EigenvectorPC2 = as.vector(pc_values$rotation[,2])) %>% dplyr::arrange(desc(ContributionPC1))

      contrib_tab = flextable::flextable(contrib_db) %>%
        flextable::set_caption(caption = "Table 2. Contribution of variables to principal components") %>%
        flextable::set_table_properties(layout = "autofit", width = 1)

      # Save the three results (1 plots, 2 table) in the word document
      results_doc = results_doc %>%
        officer::body_add_img(sr = file.path(img_path, paste(varplot_name, ".png", sep = "")),
                              width = varplot_width, height = varplot_height) %>%
        officer::body_add_par(value = varplot_name, style = "graphic title") %>%
        officer::body_add_break() %>%
        officer::body_add_par("Statistical Tables", style = "heading 2") %>%
        flextable::body_add_flextable(value = correl_tab, align = "left", topcaption = TRUE, split = FALSE) %>%
        officer::body_add_break() %>%
        flextable::body_add_flextable(value = contrib_tab, align = "left", topcaption = TRUE, split = FALSE) %>%
        officer::body_add_par("The contribution describes how much (in %) of the PC is composed of that variable; a more precise definition is that it is the variable's eigenvector^2 when the sum of eigenvectors^2 (over all variables) is equal to 1 (100%) for each PC (this condition is implicit in PCA). The eigenvector represents both the magnitude and direction of change of the principal component for each unit increase of the variable. This variable is in the standardized scale (z-score normalized) if scale and center arguments are set to TRUE.", style = "Normal") %>%
        officer::body_add_break()
    }

    # Create long database for boxplots
    if(missing_method == "na_omit") {
      if(treatment_conds_i %in% c("none")){base_factor <- interaction(multivar_db_ori[, factors_cols], sep = " ")}
      if(treatment_conds_i %in% c("pooled2", "facet2")) {
        base_factor = multivar_db_ori[, factors_cols[1]]
        second_factor = multivar_db_ori[, factors_cols[2]]
      }
      if(treatment_conds_i %in% c("pooled1", "facet1")) {
        base_factor = multivar_db_ori[, factors_cols[2]]
        second_factor = multivar_db_ori[, factors_cols[1]]
      }
    } else {
      base_factor = plot_db$base_factor
      second_factor = plot_db$second_factor}

    # Create long database from plot_db
    plot_db_long = data.frame(cbind(base_factor = base_factor,
                                    second_factor = if(treatment_conds_i == "none"){"none"} else {second_factor},
                                    matrix_values_ori) %>%
                                tidyr::pivot_longer(cols = 3:(2 + length(values_cols)),
                                                    names_to = "variable", values_to = "values"))
    plot_db_long$variable = gsub("\\.", " ", plot_db_long$variable)

    # Create boxplots with variables facet
    point_size_algo = 2.4 - (facet_col * length(unique(plot_db_long$base_factor)) - 4)/ (33 / 0.9)
    box_height = 1.7 + 1.1 * facet_row + ifelse(boxplot_legend == FALSE, 0.8, 0)

    boxplot$data = plot_db_long
    boxplot$mapping = aes(y = values, x = base_factor, colour = base_factor)
    boxplot$guides = guides(colour = guide_legend(base_factor_name), fill = guide_legend(base_factor_name),
                            shape = guide_legend(base_factor_name))
    boxplot$layers[[2]] = geom_jitter(width = 0.2, height = 0, shape = 21, size = point_size_algo, na.rm = TRUE)
    boxplot = boxplot +
      facet_wrap(~ variable, ncol = facet_col, nrow = facet_row, scales = "free_y") +
      ylab("Value") +
      theme(plot.margin = unit(c(0.07, 0.09, 0.07, 0.09), "in"),
            strip.text = element_text(size = 8.8, margin = margin(0.155, 0.1, 0.155, 0.1, unit = "cm")),
            strip.background = element_rect(fill = "grey", color = NA))
    if(boxplot_x_lab == TRUE) {boxplot <- boxplot + xlab(base_factor_name)}
    if(boxplot_filled == TRUE) {boxplot$mapping <- aes(y = values, x = base_factor, fill = base_factor)}
    if(boxplot_legend == FALSE) {boxplot = boxplot + theme(legend.position = "none")}

    if(treatment_conds_i %in% c("pooled1", "pooled2")){
      boxplot$mapping = aes(y = values, x = base_factor, fill = second_factor)
      boxplot$guides = guides(colour = guide_legend(second_factor_name), fill = guide_legend(second_factor_name),
                              shape = guide_legend(second_factor_name))
      boxplot$layers[[2]] = geom_point(position = position_jitterdodge(jitter.width = 0.25, jitter.height = 0),
                                       shape = 21, size = point_size_algo, na.rm = TRUE)

      save_name = ifelse(treatment_conds_i == "pooled1", "facetvariables1", "facetvariables2")
      boxplot_name = paste(sep = "", "Boxplot of Values Facet-Variables Group-", second_factor_name)
      boxplot_list[[treatment_conds_i]][[save_name]] = boxplot

      if(plot_out_pptx == TRUE) {
        eoffice::topptx(figure = boxplot, filename = pptx_name,
                        width = 6.4, height = box_height, append = TRUE)
      }

      ggsave(plot = suppressWarnings(boxplot), filename = paste(sep = "", boxplot_name, ".png"), path = img_path,
             dpi = 600, width = 6.4, height = box_height)
      results_doc_boxplot = results_doc_boxplot %>%
        officer::body_add_img(sr = file.path(img_path, paste(sep = "", boxplot_name, ".png")),
                              width = 6.4, height = box_height) %>%
        officer::body_add_par(value = boxplot_name, style = "graphic title") %>%
        officer::body_add_par("") %>%
        officer::body_add_par("")

      #Return to original params.
      boxplot$mapping = aes(y = values, x = base_factor, fill = base_factor)
      boxplot$guides = guides(colour = guide_legend(base_factor_name), fill = guide_legend(base_factor_name),
                              shape = guide_legend(base_factor_name))
      boxplot$layers[[2]] = geom_jitter(width = 0.2, height = 0, shape = 21, size = point_size_algo, na.rm = TRUE)
    }

    if(treatment_conds_i %in% c("none", "pooled1", "pooled2")) {
      # Save the boxplots with variables facet
      boxplot_name = ifelse(treatment_conds_i == "none", "Boxplot of Values Facet-Variables",
                            paste("Boxplot of Values Facet-Variables Pooled Across-", second_factor_name, sep = ""))
      boxplot_list[[treatment_conds_i]][["facetvariables"]] = boxplot

      if(plot_out_pptx == TRUE) {
        eoffice::topptx(figure = boxplot, filename = pptx_name,
                        width = 6.4, height = box_height, append = TRUE)
      }

      ggsave(plot = suppressWarnings(boxplot), filename = paste(sep = "", boxplot_name, ".png"), path = img_path,
             dpi = 600, width = 6.4, height = box_height)
      results_doc_boxplot = results_doc_boxplot %>%
        officer::body_add_img(sr = file.path(img_path, paste(sep = "", boxplot_name, ".png")),
                              width = 6.4, height = box_height) %>%
        officer::body_add_par(value = boxplot_name, style = "graphic title") %>%
        officer::body_add_par("") %>%
        officer::body_add_par("")
    }

    if(treatment_conds_i %in% c("facet1", "facet2")) {
      results_doc_boxplot = results_doc_boxplot %>%
        officer::body_add_par(paste("NOTE: Plots with facets by",
                                    second_factor_name, "can only be created by removing facets by variables. To achieve this, use 'boxplot_var_sep' == TRUE.")) %>%
        officer::body_add_par("")
    }

    # Create boxplots without variables facet, for each separate variable.
    if(boxplot_var_sep == TRUE) {
      box_num = length(unique(plot_db_long$base_factor)) * ifelse(treatment_conds_i %in% c("facet1", "facet2"), facet_col, 1)
      box_row = ifelse(treatment_conds_i %in% c("facet1", "facet2"), facet_row, 1)
      point_size_algo = 3 - (box_num - 2)/ 24

      box_height2 = 2.05 + (1.1 * box_row) + ifelse(box_num <= 6, 0.4, 0) + ifelse(boxplot_legend == FALSE, 0.8, 0)
      for(values_cols_i in values_cols) {
        varname = gsub("\\.", " ", colnames(multivar_db)[values_cols_i])

        # Filter for a variable
        plot_db_long_filtered = plot_db_long[plot_db_long$variable == varname, ]

        # Change boxplot data, point size, facet and names
        boxplot = boxplot +
          ylab(varname) +
          theme(strip.text = element_blank(),
                strip.background = element_blank(),
                plot.margin = unit(c(0.0 + (0.04 * box_num), 0.25, 0.05 + (0.01 * box_num), 0.25), "in"))
        boxplot$data = plot_db_long_filtered
        boxplot$layers[[2]] = geom_jitter(width = 0.2, height = 0, shape = 21, size = point_size_algo, na.rm = TRUE)

        if(treatment_conds_i %in% c("facet2", "facet1")) {
          boxplot$facet = facet_wrap(~ second_factor, ncol = facet_col, nrow = facet_row)
          boxplot = boxplot +
            theme(strip.text = element_text(size = 8.8, margin = margin(0.155, 0.1, 0.155, 0.1, unit = "cm")),
                  strip.background = element_rect(fill = "grey", color = NA),
                  plot.margin = unit(c(0.0 + (0.04 * box_num), 0.25, 0.05 + (0.01 * box_num), 0.25), "in"))
          if(boxplot_legend == FALSE) {boxplot = boxplot + theme(legend.position = "none")}
          boxplot_name = paste("Boxplot of ", varname, " Facet-", second_factor_name, sep = "")
        }
        if(treatment_conds_i == "none") {boxplot_name <- paste("Boxplot of ", varname, sep = "")}
        if(treatment_conds_i %in% c("pooled2", "pooled1")) {boxplot_name <- paste("Boxplot of ", varname, " Pooled Across-", second_factor_name, sep = "")}
        boxplot_name_vec = c(boxplot_name_vec, boxplot_name)
        box_height2_vec = c(box_height2_vec, box_height2)

        # Save these boxplots
        boxplot_list[[treatment_conds_i]][[varname]] = boxplot
        if(plot_out_pptx == TRUE) {
          eoffice::topptx(figure = boxplot, filename = pptx_name,
                          width = 6.4, height = box_height2, append = TRUE)
        }
        ggsave(plot = suppressWarnings(boxplot), filename = paste(boxplot_name, ".png", sep = ""),
               dpi = 600, width = 6.4, height = box_height2, path = img_path)
      }

      # Add them to results_doc
      if(treatment_conds_i == last(treatment_conds)){
        i = 0
        for(boxplot_name_i in boxplot_name_vec) {
          i = i + 1
          if(i == 1) {
            results_doc_boxplot = results_doc_boxplot %>%
              officer::body_add_par(value = "Variables in Separate Plots", style = "heading 2")
          }

          results_doc_boxplot = results_doc_boxplot %>%
            officer::body_add_img(sr = file.path(img_path, paste(boxplot_name_i, ".png", sep = "")),
                                  width = 6.4, height = box_height2_vec[i]) %>%
            officer::body_add_par(value = boxplot_name_i, style = "graphic title") %>%
            officer::body_add_par("") %>%
            officer::body_add_par("")
        }
      }
    }

    # Create LDA, MANOVA and potentially ANOVA tables
    if(last(treatment_conds) == treatment_conds_i) {
      results_doc_boxplot = results_doc_boxplot %>% officer::body_add_break()

      # Create LDA table results
      lda_db = cbind(multivar_db[, factors_cols, drop = FALSE],
                     scale(matrix_values[,!colnames(matrix_values) %in% c("PC1", "PC2")],
                           scale = scale, center = center))
      factors_vec = paste(collapse = " * ", colnames(multivar_db[, factors_cols, drop = FALSE]))
      values_vec = paste(collapse = " + ", colnames(lda_db[, (1 + length(factors_cols)):ncol(lda_db)]))

      tab_counter = 2
      for(lda_factor in colnames(multivar_db[, factors_cols, drop = FALSE])){
        tab_counter = tab_counter + 1

        # Fit LDA mod
        lda_mod = MASS::lda(formula = as.formula(paste(lda_factor, "~", values_vec)),
                            data = lda_db)

        # Store LDA results in a table and list
        lda_contrib = (scale(lda_mod$scaling, F, sqrt(colSums(lda_mod$scaling^2))))^2
        lda_contrib_db = data.frame(Variable = gsub("\\.", " ", rownames(lda_contrib)),
                                    ContributionLD1 = 100 * lda_contrib[, 1],
                                    CoefficientLD1 = lda_mod$scaling[, 1])
        if(ncol(lda_contrib) > 1) {
          lda_contrib_db$ContributionLD2 = 100 * lda_contrib[, 2]
          lda_contrib_db$CoefficientLD2 = lda_mod$scaling[, 2]
        }

        lda_contrib_db = lda_contrib_db %>% dplyr::arrange(dplyr::desc(ContributionLD1))
        lda_contrib_tab = flextable::flextable(lda_contrib_db) %>%
          flextable::set_caption(caption = paste("Table ", tab_counter, ". Contribution of variables to linear discriminants of ",
                                                 lda_factor, sep = "")) %>%
          flextable::set_table_properties(layout = "autofit", width = 1)
        lda_exp = round(100 * (lda_mod$svd^2 / sum(lda_mod$svd^2)), 1)

        # Add results to officer
        results_doc = results_doc %>%
          flextable::body_add_flextable(value = lda_contrib_tab, align = "left", topcaption = TRUE, split = FALSE) %>%
          officer::body_add_par(value = paste(sep = "", "LD1 and LD2 explain ", lda_exp[1], "% and ", lda_exp[2],
                                              "% of the separation between ", lda_factor, ", respectively."), style = "Normal") %>%
          officer::body_add_par("")
      }

      results_doc = results_doc %>%
        officer::body_add_par("Tabulated values are based on Linear Discriminant Analysis (LDA). The contribution describes the percentage of the linear discriminant composed of that variable; a more exact definition is that it is the variable's coefficient of linear discriminant^2 when the sum of coefficients^2 (across variables) is normalized to 100%. The coefficient of linear discriminant describes the magnitude and direction of change in the linear discriminant score (towards a particular factor level) per unit increase of the variable. Apart from the coefficients relative signs and magnitude, they may be difficult to interpret without an associated LDA plot.", style = "Normal")

      # Get MANOVA results
      manova_db = cbind(multivar_db[, factors_cols, drop = FALSE],
                        matrix_values[,!colnames(matrix_values) %in% c("PC1", "PC2")])
      values_vec2 = paste("cbind(",
                          paste(collapse = ",", colnames(manova_db[, (1 + length(factors_cols)):ncol(manova_db)])),
                          ")", sep = "")

      options(contrasts = c("contr.sum", "contr.poly"))
      permanova = suppressWarnings(RVAideMemoire::adonis.II(formula = as.formula(paste(values_vec, "~", factors_vec)),
                                                            data = manova_db, method = "euclidean"))
      manova_mod = manova(formula = as.formula(paste(values_vec2, "~", factors_vec)),
                          data = manova_db)
      manova = car::Anova(type = 3, manova_mod)
      options(contrasts = c("contr.treatment", "contr.poly"))
      manova_p = suppressWarnings(as.numeric(na.omit(as.numeric(sub(".*\\s*(\\d+\\.[0-9e-]+)\\s*[*.]*",
                                                                    "\\1", capture.output(manova))))))
      # Create MANOVA table
      manova_db = data.frame(Factor = manova$terms[!manova$terms == "(Intercept)"],
                             PERMANOVA = as.numeric(na.omit(permanova$`Pr(>F)`)),
                             P_value = sapply(manova_p[!manova$terms == "(Intercept)"], knife))
      manova_db$PERMANOVA = ifelse(manova_db$PERMANOVA < 0.001, formatC(manova_db$PERMANOVA, format = "e", digits = 4), manova_db$PERMANOVA)
      manova_db = manova_db[, -2] #exclude PERMANOVA results for now to prevent confusion. May add some time in the future.

      manova_tab = flextable::flextable(manova_db) %>%
        flextable::set_caption(caption = paste("Table ", tab_counter + 1, ". MANOVA Results", sep = "")) %>%
        flextable::set_table_properties(layout = "autofit", width = 1)

      results_doc = results_doc %>%
        officer::body_add_break() %>%
        flextable::body_add_flextable(value = manova_tab, align = "left", topcaption = TRUE, split = FALSE) %>%
        #officer::body_add_par(ifelse(length(factors_cols) == 2, "Type 3 SS MANOVA was conducted.", "")) %>%
        officer::body_add_par("") %>%
        officer::body_add_break()

      # Run Univariate Tests
      if(univariate_tests == TRUE) {
        p_no_adj = c()
        contrast_db = data.frame()
        contrast_db2 = data.frame()
        cld_db = data.frame()
        cld_db2 = data.frame()

        # Get p-values for each variable
        varnames_vec = colnames(multivar_db[, values_cols, drop = FALSE])
        options(contrasts = c("contr.sum", "contr.poly"))
        for(values_vec_i in varnames_vec) {
          # ANOVA and model fitting
          lm_mod = lm(data = multivar_db, formula = as.formula(paste(values_vec_i, "~", factors_vec)))
          Anova_res = car::Anova(lm_mod, type = 3, test.statistic = "F")
          p_no_adj = c(p_no_adj, as.numeric(na.omit(Anova_res$`Pr(>F)`[-1])))

          # Pairwise comparisons
          factor_conds = colnames(multivar_db[, factors_cols, drop = FALSE])
          for(factor_conds_i in factor_conds){
            ## One Factor situations
            if(length(factors_cols) == 1){
              emmeans_obj = emmeans::emmeans(lm_mod, as.formula(paste("pairwise ~", factor_conds_i)), adjust = "BH")
              cld_obj = data.frame(multcomp::cld(emmeans_obj, Letter = letters, adjust = "BH"))
            }
            ## Two Factor situations
            if(length(factors_cols) == 2) {
              emmeans_obj = emmeans::emmeans(lm_mod, as.formula(paste("pairwise ~", factor_conds_i)), adjust = "BH")
              emmeans_obj2 = emmeans::emmeans(lm_mod, as.formula(paste("pairwise ~", factor_conds_i,
                                                                       "|", factor_conds[factor_conds != factor_conds_i])))
              cld_obj = data.frame(multcomp::cld(emmeans_obj, Letter = letters, adjust = "BH"))
              cld_obj2 = data.frame(multcomp::cld(emmeans_obj2, Letter = letters, adjust = "BH"))

              cld_db2 = rbind(cld_db2, data.frame(Factor = factor_conds_i,
                                                  Variable = values_vec_i,
                                                  Condition = cld_obj2[, 2],
                                                  Compared_Levels = cld_obj2[, 1],
                                                  Statistical_Classes = cld_obj2[, 8]))

              contrast_db2 = rbind(contrast_db2, data.frame(Factor = factor_conds_i,
                                                            Variable = values_vec_i,
                                                            Condition = data.frame(emmeans_obj2$contrasts)[, 2],
                                                            Contrast = data.frame(emmeans_obj2$contrasts)[, 1],
                                                            Difference = data.frame(emmeans_obj2$contrasts)[, 3],
                                                            P_value = sapply(data.frame(emmeans_obj2$contrasts)[, 7], knife)))
            }

            # Create cld db
            cld_db = rbind(cld_db, data.frame(Factor = factor_conds_i,
                                              Variable = values_vec_i,
                                              Levels = cld_obj[, 1],
                                              Statistical_Classes = cld_obj[, 7]))

            # Create contrast db
            contrast_db = rbind(contrast_db, data.frame(Factor = factor_conds_i,
                                                        Variable = values_vec_i,
                                                        Contrast = data.frame(emmeans_obj$contrasts)[, 1],
                                                        Difference = data.frame(emmeans_obj$contrasts)[, 2],
                                                        P_value = sapply(data.frame(emmeans_obj$contrasts)[, 6], knife)))
          }
        }

        options(contrasts = c("contr.treatment", "contr.poly"))
        p_bh = p.adjust(p_no_adj, method = "BH")

        # Create ANOVA p-values table
        anova_fac = rownames(data.frame(Anova_res))[-c(1, length(rownames(data.frame(Anova_res))))]
        p_db = data.frame(Variable = rep(varnames_vec, each = length(anova_fac)),
                          Factor = rep(anova_fac, times = length(varnames_vec)),
                          P_value = sapply(p_no_adj, knife))
        p_db$Variable = gsub("\\.", " ", p_db$Variable)

        colnames(p_db) = gsub("_", " ", colnames(p_db))
        p_tab = flextable::flextable(p_db) %>%
          flextable::set_caption(caption = paste("Table ", tab_counter + 2, ". ANOVA Results", sep = "")) %>%
          flextable::set_table_properties(layout = "autofit", width = 1) %>%
          flextable::merge_v(j = "Variable")

        # Create cld table (simple)
        cld_db = cld_db[order(cld_db$Factor),]
        cld_db$Variable = gsub("\\.", " ", cld_db$Variable)

        colnames(cld_db) = gsub("_", " ", colnames(cld_db))
        cld_tab = flextable::flextable(cld_db) %>%
          flextable::set_caption(caption = paste("Table ", tab_counter + 3,
                                                 ". Statistical Classes from Simple Comparisons", sep = "")) %>%
          flextable::set_table_properties(layout = "autofit", width = 1) %>%
          flextable::merge_v(j = c("Factor", "Variable"))

        # Create contrast table (simple)
        contrast_db = contrast_db[order(contrast_db$Factor),]
        contrast_db$Variable = gsub("\\.", " ", contrast_db$Variable)

        colnames(contrast_db) = gsub("_", " ", colnames(contrast_db))
        contrast_tab = flextable::flextable(contrast_db) %>%
          flextable::set_caption(caption = paste("Table ", tab_counter + 4,
                                                 ". Simple Pairwise Comparison Results", sep = "")) %>%
          flextable::colformat_double(digits = 5) %>%
          flextable::set_table_properties(layout = "autofit", width = 1) %>%
          flextable::merge_v(j = c("Factor", "Variable"))

        if(length(factors_cols) == 2){

          # Create cld table (conditional)
          cld_db2 = cld_db2[order(cld_db2$Factor),]
          cld_db2$Variable = gsub("\\.", " ", cld_db2$Variable)

          colnames(cld_db2) = gsub("_", " ", colnames(cld_db2))
          cld_db2 = cld_db2[, -1]
          cld_tab2 = flextable::flextable(cld_db2) %>%
            flextable::set_caption(caption = paste("Table ", tab_counter + 5,
                                                   ". Statistical Classes from Conditional Comparisons", sep = "")) %>%
            flextable::set_table_properties(layout = "autofit", width = 1) %>%
            flextable::merge_v(j = c("Variable", "Condition"))

          # Create contrast table (conditional)
          contrast_db2 = contrast_db2[order(contrast_db2$Factor),]
          contrast_db2$Variable = gsub("\\.", " ", contrast_db2$Variable)

          colnames(contrast_db2) = gsub("_", " ", colnames(contrast_db2))
          contrast_db2 = contrast_db2[, -1]
          contrast_tab2 = flextable::flextable(contrast_db2) %>%
            flextable::set_caption(caption = paste("Table ", tab_counter + 6,
                                                   ". Conditional Pairwise Comparison Results", sep = "")) %>%
            flextable::colformat_double(digits = 5) %>%
            flextable::set_table_properties(layout = "autofit", width = 1) %>%
            flextable::merge_v(j = c("Variable"))
        }

        # Add results to word output
        results_doc_boxplot = results_doc_boxplot %>%
          officer::body_add_par("Statistical Tables", style = "heading 2") %>%
          flextable::body_add_flextable(value = p_tab, align = "left", topcaption = TRUE, split = FALSE) %>%
          officer::body_add_par(ifelse(length(factors_cols) == 2, "Type 3 SS ANOVAs were conducted.", "")) %>%
          officer::body_add_par("") %>%
          flextable::body_add_flextable(value = cld_tab, align = "left", topcaption = TRUE, split = FALSE) %>%
          officer::body_add_par(paste("Different letters denote significantly different levels of a factor. All levels within a factor were compared pairwise using the function emmeans() with the Benjamini Hocherbg correction to control for a False Discovery Rate of 5%.", ifelse(length(factors_cols) == 2, "emmeans() were based on two-factor models (with interactions) fitted separately by variable.", ""))) %>%
          officer::body_add_par("") %>%
          flextable::body_add_flextable(value = contrast_tab, align = "left", topcaption = TRUE, split = FALSE) %>%
          officer::body_add_par(paste("P-values were generated from pairwise comparisons of levels of a factor. Comparisons conducted using the function emmeans() with the Benjamini Hochberg correction to control for a False Discovery Rate of 5%.", ifelse(length(factors_cols) == 2, "emmeans() were based on two-factor models (with interactions) fitted separately by variable.", ""))) %>%
          officer::body_add_par("")

        if(length(factors_cols) == 2) {
          results_doc_boxplot = results_doc_boxplot %>%
            flextable::body_add_flextable(value = cld_tab2, align = "left", topcaption = TRUE, split = FALSE) %>%
            officer::body_add_par("Different letters denote significantly different levels of a factor. Levels were compared pairwise within a factor while holding a level from another factor constant. Comparisons conducted using the function emmenas() with the Benjamini Hochberg correction to control for a False Discovery Rate of 5%. emmeans() were based on two-factor models (with interactions) fitted separately by variable.") %>%
            officer::body_add_par("") %>%
            flextable::body_add_flextable(value = contrast_tab2, align = "left", topcaption = TRUE, split = FALSE) %>%
            officer::body_add_par("P-values were generated from pairwise comparisons of levels of a factor, while holding a level from another factor constant. Comparisons conducted using the function emmenas() with the Benjamini Hochberg correction to control for a False Discovery Rate of 5%. emmeans() were based on two-factor models (with interactions) fitted separately by variable.")
        }
      }
    }
  }

  # Save word document
  print(results_doc_boxplot, target = tempf)
  print(officer::body_add_docx(results_doc, src = tempf),
        target = paste(sep = "", "MultiVar Report ", format(Sys.Date(), "%d%b%Y"), ".docx"))

  # Delete temporary files
  file.remove(tempf)
  if(plot_out_png == FALSE) {unlink(img_path)}

  # Print plot output in R
  if(plot_out_R == TRUE) {
    return(list(pca = pca_list,
                correl = varplot,
                box = boxplot_list))
  }
}

##################################################### Data 1 - mort_db_ex #######################################################

#' Example Mort Data
#'
#' @description A subset of columns taken from the online excel mortality file in OneDrive.
#' @usage
#' data(mort_db_ex)
#' View(mort_db_ex)
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
#' View(surv_db_ex)
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
#' @description A reference hazard dataframe created using \code{Surv_Plots(..., dailybin = TRUE, data_out = TRUE)$Hazard_DB} which uses \code{bshazard::bshazard()} inside. Contains hazard rates over time.
#' @usage
#' data(haz_db_ex)
#' View(haz_db_ex)
#'
#' @format A data frame containing 54 rows and 3 columns:\tabular{lll}{
#'  \code{Trt.ID} \tab \tab A label for the treatment group used in creating this reference hazard dataframe \cr
#'  \code{Hazard} \tab \tab Hazard values (rates) \cr
#'  \code{Time} \tab \tab Time / TTE in days post challenge. \cr
#' }
#'
"haz_db_ex"

##################################################### Data 4 - multivar_db_ex #######################################################

#' Example Multivariate Data
#'
#' @description An example dataset representing the mucus chemistry of fish exposed to different treatments. Mucus chemistry is described by several variables from columns 2 to 8. In column 9, a fake factor "Fruits" was added for use as an example second factor. In column 4, row 10, the available value was replaced by NA to use as an example missing value.
#'
#' @usage
#' data(multivar_db_ex)
#' View(multivar_db_ex)
#'
"multivar_db_ex"

##################################################### Data 5 - surv_sim_db_ex #######################################################

#' Example Simulated Survival Object
#'
#' @description A list generated from \code{Surv_Simul()} containing the following elements: a simulated survival dataset (\code{surv_simul_db}), the dataset describing the survival characteristics of the population (\code{surv_pop_db}) and a plot illustrating both (\code{surv_plots}). The specified simulation parameters (arguments) in \code{Surv_Simul()} are \code{fish_num_per_tank = 100}, \code{tank_num_per_trt = 4}, \code{treatments_hr = c(1, 0.7, 0.5)}, \code{logHR_sd_intertank = 0}, \code{n_sim = 10}, \code{plot_out = TRUE}.
#'
#' @usage
#' data(surv_sim_db_ex)
#' View(surv_sim_db_ex$surv_simul_db)
#' View(surv_sim_db_ex$surv_pop_db)
#' View(surv_sim_db_ex$surv_plots)
#'
"surv_sim_db_ex"
