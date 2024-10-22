# This is an R package containing useful functions for my work.

# Available functions with documentation:
# 1. Con_Simul() -- simulates contingency tables based on the multinomial distribution.
# 1b. Con_Simul_PR() -- calculates positive rates for statistical tests on contingency tables.
# 6. Surv_Gen() -- generate rows of survivors given a starting number of fish per tank and data containing morts and sampled fish.
# 7. Surv_Plots() -- generate Kaplan-Meier survival curve and hazard curve from survival data.
# 8. GG_Color_Hue() -- returns the default colour codes assigned by ggplot to a given number of categorical groups (n)

# Available functions without documentation:
# 2. Simul_Con_MULT.FISH.ORD() -- simulates ordinal-distributed data across treatments and lesions with inter-fish variation in the PO.
# 3. Simul_Surv() -- simulate survival data based on a reference hazard function, the specified hazard ratio(s), and inter-tank variation.
# 4. theme_Publication() -- ggplot theme for generating publication-ready plots.
# 5. Surv_Pred() -- predict future survival rate(s) for ongoing experiment based on a reference hazard function from older data.

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
#' pwr::pwr.chisq.test(w = ES.w2(probs_mat), df = 2, sig.level = 0.05, N = 100)
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
#' @description Simulates survival data based on a set of user-specified experimental parameters and a reference hazard curve (e.g. hazard curve from the control group). Able to simulate data with inter-cluster (e.g. tank) variation based on the framework of the mixed cox proportional hazards model (\code{coxme::coxme}). Able to simulate right censored data (e.g. sampled fish) using \code{sampling_specs} argument. Optionally produces a plot illustrating the characteristics of the simulated data and that of the population / truth from which the data (sample) is simulated.
#'
#' @details Simulations are based on uniform-probability draws (\emph{U} ~ (0, 1)) from a set of events which can be expressed through time using the cumulative density function of failures (\emph{F(t)}, i.e. cumulative mort. curve). \emph{F(t)} can be transformed to the cumulative hazard function \emph{H(t)}, hence the relationship between \emph{H(t)} and uniform draws (from \emph{U}) is also known (derivation and equation in \href{https://epub.ub.uni-muenchen.de/1716/1/paper_338.pdf}{Bender et al. (2003)}. Because \emph{H(t)} is related (as the integral) to the hazard function \emph{h(t)}, and since \emph{h(t)} is related to effects (e.g. treatment or tank) based on the cox proportional hazards model, such effects can now be incorporated into the simulation process as they can be interacted with \emph{U}. The simulation process is as follows:
#'
#' \enumerate{
#' \item \code{Surv_Simul()} takes a random sample from \emph{U} (e.g. 0.7).
#'
#' \item U is then transformed into \emph{H} as their relationship is known (see \href{https://epub.ub.uni-muenchen.de/1716/1/paper_338.pdf}{Bender et al. 2003}) as:
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
#' @param haz_db A dataframe representing the reference hazard curve; can be generated from \code{bshazard()} or \code{Surv_Plots()}.
#' @param fish_num_per_tank The number of fish to simulate per tank. Defaults to 100.
#' @param tank_num_per_trt The number of tanks to simulate per treatment group. Defaults to 4.
#' @param treatments_hr A vector representing the hazard ratios of the treatment groups starting with the reference/control (HR = 1). Length of the vector represents the number of treatment groups. Defaults to \code{c(1, 1, 1, 1)}.
#' @param logHR_sd_intertank The standard deviation of inter-tank variation in the log(HR) scale according to the \code{coxme} framework. Defaults to 0 (no inter-tank variation) which has been and quite oftenly, the estimate for injected Trojan fish data. For reference 0.1 reflects a low inter-tank variation situation, while 0.35 is fairly high but can and has occurred in some experiments.
#' @param sampling_specs A dataframe representing the number / amount of right censored data (e.g. sampled fish) per tank at different times represented by two columns "Amount" and "TTE", respectively. See \bold{Examples} for example of use. Defaults to NULL (no sampling).
#' @param n_sim Number of survival dataset to simulate. Defaults to 1.
#' @param plot_out Whether to output the information plot (further details in \bold{return}). Defaults to TRUE.
#' @param pop_out Whether to output a dataframe containing the survival probability values for the population. Defaults to TRUE.
#' @param plot_name Character string specifying the name of the saved plot. Defaults to "Surv_Simul-Plot-Output".
#' @param theme Character string specifying the graphics theme for the plots. Theme "ggplot2" and "prism" currently available. Defaults to "ggplot2".
#'
#' @return At minimum, returns a simulated survival dataframe consisting of 5 columns: TTE (Time to Event), Status (0 / 1), Trt.ID, Tank.ID, and n_sim which represents the simulation number for the data subsets.
#'
#' If \code{plot_out = TRUE}, outputs a list that additionally contains a Kaplan-Meier survival plot. The plot illustrates the survival curve with end survival rates for the simulated dataset as well as the population. If the number of simulated dataset is greater than 1, multiple curves are drawn representing each and a statement is provided regarding the power to detect the effect of Treatment -- specifically, the percent positive (p < 0.05) from a global log-rank test using \code{survival::survdiff()}.
#'
#' If \code{pop_out = TRUE}, outputs a list that additionally contains a dataframe representing the survival probability values for the population / truth from which the sample is supposedly taken.
#'
#' @import ggplot2
#' @import magrittr
#' @import dplyr
#'
#' @seealso \href{file:///C:/Users/sean4/Documents/GitHub/safuncs/docs/reference/Surv_Simul.html}{Link} for executed examples which includes any figure outputs.
#'
#' @export
#'
#' @examples
#' #Starting from an example mortality database, we first generate the complete survivor data using Surv_Gen()
#' data(mort_db_ex)
#' surv_dat = Surv_Gen(mort_db = mort_db_ex,
#'                     starting_fish_count = 100,
#'                     last_tte = 54)
#'
#' #Filter for the control group ("A") to use as a reference hazard curve for simulations
#' surv_dat_A = surv_dat[surv_dat$Trt.ID == "A", ]
#'
#' #Estimate the hazard curve of the control group and get the associated hazard dataframe using either bshazard::bshazard() or safuncs::Surv_Plots()$Hazard_DB
#' ref_haz_route_bshazard = bshazard::bshazard(data = surv_dat_A, survival::Surv(TTE, Status) ~ Tank.ID, verbose = FALSE)
#' ref_haz_route_bshazard = data.frame(summary(ref_haz_route_bshazard)$HazardEstimates)
#'
#' ref_haz_route_safuncs = safuncs::Surv_Plots(surv_db = surv_dat_A, data_out = TRUE)$Hazard_DB
#'
#' #Simulate!
#' Surv_Simul(haz_db = ref_haz_route_safuncs,
#'            treatments_hr = c(1, 0.8, 0.5),
#'            sampling_specs = data.frame(Amount = 10,
#'                                        TTE = 45))$surv_plots #sampling of 10 fish at 45 DPC, but otherwise default experimental conditions.
#'
#' #Simulate multiple times to better see if samples are reliable to answer the question: are my future samples likely to be good approximates of the truth / population
#' Surv_Simul(haz_db = ref_haz_route_safuncs,
#'            treatments_hr = c(1, 0.8, 0.5),
#'            sampling_specs = data.frame(Amount = 10,
#'                                        TTE = 45),
#'            n_sim = 4)$surv_plots
Surv_Simul = function(haz_db,
                      fish_num_per_tank = 100,
                      tank_num_per_trt = 4,
                      treatments_hr = c(1, 1, 1, 1),
                      logHR_sd_intertank = 0,
                      sampling_specs = NULL,
                      n_sim = 1,
                      plot_out = TRUE,
                      pop_out = TRUE,
                      plot_name = "Surv_Simul-Plot-Output",
                      theme = "ggplot2") {

  #Initialize objects to store loop results
  surv_samps = data.frame()
  Surv_simul_outDB = data.frame()
  pvalues = c()
  colnames(haz_db) = tolower(colnames(haz_db))

  #Simulate survival dataframe
  for(loopnum in 1:n_sim) {

    CDF_Yval = c()

    for(Treatment_Term in treatments_hr) {

      for(Tank_num in 1:tank_num_per_trt) {

        #Random sampling
        Tank_eff = rnorm(n = 1, mean = 0, sd = logHR_sd_intertank)

        U = runif(n = fish_num_per_tank, min = 0, max = 1)

        CDF_Yval_temp = -log(U) * exp(-(log(Treatment_Term) + Tank_eff))
        CDF_Yval = append(CDF_Yval, CDF_Yval_temp)
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
                               Trt.ID = as.factor(rep(c("Control", LETTERS[2:length(treatments_hr)]),
                                                      each = fish_num_per_tank * tank_num_per_trt, times = n_sim)),
                               Tank.ID = as.factor(rep(1:(length(treatments_hr) * tank_num_per_trt),
                                                       each = fish_num_per_tank, times = n_sim)),
                               n_sim = loopnum)

    #Transform TTE and Status in certain rows due to sampling
    if(!is.null(sampling_specs)) {

      #Repeat for each specified sampling time and each tank
      for(samp_time in 1:length(sampling_specs$TTE)) {
        for(tank_num in levels(Surv_simul_DB$Tank.ID)) {

          rows_sel = sample(x = which(Surv_simul_DB[Surv_simul_DB$Tank.ID == tank_num, "TTE"] > sampling_specs$TTE[samp_time]),
                            size = sampling_specs$Amount[samp_time],
                            replace = FALSE)

          Surv_simul_DB$TTE[rows_sel] = sampling_specs$TTE[samp_time]
          Surv_simul_DB$Status[rows_sel] = 0
        }
      }
    }

    #Get p-value for plots
    pvalues = append(pvalues, survival::survdiff(survival::Surv(TTE, Status) ~ Trt.ID, Surv_simul_DB)$pvalue)

    #Simulated survival data to be provided as output
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
                                 type = paste("Sample (n = ", n_sim, ")", sep = ""),
                                 n_sim = loopnum,
                                 alpha = 1 - (0.0001 ^ (1/n_sim)))

    surv_samps_ends = data.frame(surv_samps_temp %>%
                                   dplyr::group_by(Trt.ID) %>%
                                   dplyr::reframe(surv_prob = c(min(1, max(surv_prob)), min(surv_prob)),
                                                  time = c(min(haz_db$time), max(haz_db$time)),
                                                  n_sim = loopnum,
                                                  alpha = 1 - (0.0001 ^ (1/n_sim))))
    surv_samps_ends$type = paste("Sample (n = ", n_sim, ")", sep = "")

    surv_samps = rbind(surv_samps, surv_samps_temp, surv_samps_ends)
  }

  #Get "population" survival dataset by exponentiating the negative cumulative hazard
  surv_pop = data.frame(Trt.ID = as.factor(rep(c("Control", LETTERS[2:length(treatments_hr)]), each = length(haz_db$hazard))),
                        surv_prob = exp(-as.vector(apply(haz_db$hazard %*% t(treatments_hr), 2, cumsum))),
                        time = rep(haz_db$time, times = length(treatments_hr)),
                        type = "Population / truth",
                        n_sim = 1,
                        alpha = 1)

  #To the end of creating survival plots
  surv_comb = rbind(surv_samps, surv_pop)
  surv_comb$type = factor(surv_comb$type, levels = c(paste("Sample (n = ", n_sim, ")", sep = ""), "Population / truth"))
  surv_comb$Trt.ID = factor(surv_comb$Trt.ID, levels = rep(c("Control", LETTERS[2:length(treatments_hr)])))

    #Get end_sr for population plots and sample plots
    end_db = data.frame(surv_comb %>%
                          dplyr::group_by(type, Trt.ID, n_sim) %>%
                          dplyr::summarise(surv_prob = min(surv_prob), time = max(TTE), .groups = "drop") %>%
                          dplyr::group_by(type, Trt.ID) %>%
                          dplyr::summarise(surv_prob = mean(surv_prob), time = max(time), .groups = "drop"))

    #Get % significance (i.e. power) for plotting
    perc_sf = paste(round(100 * sum(pvalues < 0.05) / length(pvalues), digits = 0), "%", sep = "")

    #Ggplot
    if(n_sim == 1) {
      surv_plots = ggplot(data = surv_comb, aes(x = time, y = surv_prob, colour = Trt.ID, group = interaction(n_sim, Trt.ID))) +
        facet_wrap(~ type) +
        geom_step() +
        scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), labels = scales::percent) +
        scale_x_continuous(breaks = seq(0, max(surv_pop$time), max(round(max(surv_pop$time) / 15), 1))) +
        ylab("Survival Probability (%)") +
        xlab("Time to Event") +
        scale_alpha(range = c(min(surv_comb$alpha), 1))+
        guides(alpha = "none") +
        geom_text(data = end_db, aes(x = time, y = surv_prob, label = round(surv_prob * 100, digits = 0)),
                  vjust = -0.3, hjust = 0.8, show.legend = FALSE, size = 3.3) +
        annotation_custom(grob = grid::textGrob(paste(c(paste("The sample has a", sep = ""),
                                                        paste("p-value = ", signif(pvalues, digits = 2), sep = ""),
                                                        "(global test of Trt.)"), collapse = "\n"),
                                                x = grid::unit(1.05, "npc"),
                                                y = grid::unit(0.08, "npc"),
                                                hjust = 0,
                                                gp = grid::gpar(fontsize = 9))) +
        coord_cartesian(clip = "off") +
        theme(plot.margin = margin(5.5, 20, 5.5, 5.5))

    } else {
      surv_plots = ggplot(data = surv_comb, aes(x = time, y = surv_prob, colour = Trt.ID, group = interaction(n_sim, Trt.ID))) +
        facet_wrap(~ type) +
        geom_step(aes(alpha = alpha)) +
        scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), labels = scales::percent) +
        scale_x_continuous(breaks = seq(0, max(surv_pop$time), max(round(max(surv_pop$time) / 15), 1))) +
        ylab("Survival Probability (%)") +
        xlab("Time to Event") +
        scale_alpha(range = c(min(surv_comb$alpha), 1)) +
        guides(alpha = "none") +
        geom_text(data = end_db[end_db$type == "Population / truth",],
                  aes(x = time, y = surv_prob, label = round(surv_prob * 100, digits = 0)),
                  vjust = -0.3, hjust = 0.8, show.legend = FALSE, size = 3.3) +
        annotation_custom(grob = grid::textGrob(paste(c(paste(perc_sf, " of the sample", sep = ""),
                                                        paste("sets (n) has p < 0.05", sep = ""),
                                                        "(global test of Trt.)"), collapse = "\n"),
                                                x = grid::unit(1.03, "npc"),
                                                y = grid::unit(0.08, "npc"),
                                                hjust = 0,
                                                gp = grid::gpar(fontsize = 9))) +
        coord_cartesian(clip = "off") +
        theme(plot.margin = margin(5.5, 20, 5.5, 5.5))

    }

  #Plot theme
  if(theme == "prism") {surv_plots = surv_plots + ggprism::theme_prism()}

  #Save plots and return outputs
  if(!is.null(plot_name)) {
    eoffice::topptx(figure = surv_plots, filename = paste(plot_name, ".pptx", sep = ""), width = 6, height = 4)
    ggsave(paste(plot_name, ".tiff", sep = ""), dpi = 900, width = 6, height = 4, plot = surv_plots)
  }

  if(plot_out == FALSE && pop_out == FALSE) {

    return(simul_surv_db = Surv_simul_outDB)

    } else {

      output = list(simul_surv_db = Surv_simul_outDB)
      if(plot_out == TRUE) {output$surv_plots <- surv_plots}
      if(pop_out == TRUE) {output$population_surv_db <- surv_pop}
      return(output)
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
#' @seealso \href{file:///C:/Users/sean4/Documents/GitHub/safuncs/docs/reference/theme_Publication.html}{Link} for executed examples which includes any figure outputs.
#'
#' @export
#'
#' @examples
#' #Load an example dataset
#' data(iris)
#'
#' #Create a ggplot modified with theme_Publication()
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
#' @description Produces survival data that includes rows for every surviving fish based on the starting number of fish and mortality data. To generate survivor data for tanks absent in the input mortality dataframe, specify the arguments \code{tank_without_mort} and \code{trt_without_mort}. To generate survivor data with tank specific starting numbers of fish, input a dataframe into the argument \code{starting_fish_count} instead of a single value; details in \code{Arguments}.
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
#'
#' @return A dataframe produced by combining the input mort data and generated rows of survivor data.
#'
#' @import magrittr
#' @import devtools
#'
#' @export
#'
#' @examples
#' #First, we load an example mortality database available from the safuncs package
#' data(mort_db_ex)
#'
#' #Next, we input this data into Surv_Gen() as well as the study details to generate entries (rows) for survivors in the output - a "complete" dataframe for further survival analysis and data visualization.
#' Surv_Gen(mort_db = mort_db_ex,
#'          starting_fish_count = 100,
#'          last_tte = 54,
#'          tank_without_mort = c("C99", "C100"),
#'          trt_without_mort = c("A", "B"))
Surv_Gen = function(mort_db,
                    starting_fish_count,
                    last_tte,
                    tank_without_mort = NULL,
                    trt_without_mort = NULL) {

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
#' For details on the statistical methodology used by \code{bshazard()}, refer to: \href{https://www.researchgate.net/publication/287338889_bshazard_A_Flexible_Tool_for_Nonparametric_Smoothing_of_the_Hazard_Function}{here}.
#'
#' General concept: h(t) the hazard function is considered in an count model with the number of deaths as the response variable. I.e, death_count(t) = h(t) * P(t) where P(t) is the number alive as a function of time and h(t) is modeled over time using basis splines. The basis spline curvature\bold{s} is assumed to have a normal distribution with mean 0 (a random effect). Based on this assumption, the author found that the variance of curvatures (i.e. smoothness) is equal to the over-dispersion (phi) of the death counts related (divided) by some smoothness parameter (lambda). Phi and lambda can be estimated from the data or specified by the user. Specification can be helpful in low sample size situations where overdispersion (phi) estimates have been found to be unreliable and clearly wrong (based on my understanding of realistic estimates and what was estimated in past data with adequate, large sample sizes).
#' @md
#'
#' @param surv_db A survival dataframe as described in \bold{Details}.
#' @param plot_prefix A string specifying the prefix for the filename of the saved plots.
#' @param xlim A vector specifying the plots x-axis lower and upper limits, respectively.
#' @param ylim A vector specifying the Survival Plot y-axis lower and upper limits, respectively.
#' @param xlab A string specifying the plot x-axis label.
#' @param lambda Smoothing value for the hazard curve. Higher lambda produces greater smoothing. Defaults to NULL where \code{bshazard()} uses the provided survival data to estimate lambda; NULL specification is recommended for large sample size situations which usually occurs on our full-scale studies with many mortalities and tank-replication. At low sample sizes, the lambda estimate can be unreliable. Choosing a lambda of 10 (or anywhere between 1-100) probably produces the most accurate hazard curve for these situations. In place of choosing lambda, choosing \code{phi} is recommended; see below.
#' @param phi Dispersion parameter for the count model used in hazard curve estimation. Defaults to NULL where \code{bshazard()} uses the provided survival data to estimate phi; NULL specification is recommended for large sample size situations. At low sample sizes, the phi estimate can be unreliable. Choosing a phi value of 1 for low sample sizes is recommended. This value of 1 (or close) seems to be that estimated in past Tenaci data (QCATC997; phi ~ 0.8-1.4) where there are large sample sizes with tank-replication. The phi value of 1 indicates the set of counts (deaths) over time have a Poisson distribution, following the different hazard rates along the curve and are not overdispersed (phi > 1).
#' @param dailybin Whether to set time bins at daily (1 TTE) intervals. Refer to the \code{bshazard()} documentation for an understanding on the role of bins to hazard curve estimation. Please set to TRUE at low sample sizes and set to FALSE at large sample sizes with tank-replication. Defaults to TRUE.
#' @param plot Which plot to output. Use "surv" for the Kaplan-Meier Survival Curve, "haz" for the Hazard Curve, or "both" for both. Defaults to "both".
#' @param colours Vector of color codes for the different treatment groups in the plot. Defaults to ggplot2 default palette.
#' @param theme Character string specifying the graphics theme for the plots. Theme "ggplot2" and "prism" currently available. Defaults to "ggplot2".
#' @param trt_order Vector representing the order of treatment groups in the plots. Defaults to NULL where alphabetical order is used.
#' @param data_out Whether to print out the survival and/or hazard databases illustrated by the plots. Defaults to FALSE.
#' @param plot_dim Vector representing the dimensions (width, height) with which to save the plot in .tiff and .pptx.
#'
#' @return By default, with argument {\code{plot = "both"}}, returns the Kaplan-Meier Survival Curve and the Hazard Curve. Output can be trimmed by setting \code{plot = "haz"} or \code{plot = "surv"}.
#'
#' If \code{data_out = TRUE}, returns a list of the plot(s) and the associated dataframes to create them.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' #Starting from an example mortality database, we first generate the complete survivor data using Surv_Gen()
#' data(mort_db_ex)
#' surv_dat = Surv_Gen(mort_db = mort_db_ex,
#'                     starting_fish_count = 100,
#'                     last_tte = 54)
#'
#' #Create plot by inputting surv_dat into Surv_Plots()!
#' Surv_Plots(surv_db = surv_dat,
#'            plot_prefix = "QCATC777",
#'            xlim = c(0, 54),
#'            ylim = c(0, 1),
#'            xlab = "TTE",
#'            phi = 1,
#'            plot = "both")
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
#' @examples
#' # Get colour codes used for 6 categorical groups
#' GG_Colour_Hue(6)
#'
GG_Colour_Hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
