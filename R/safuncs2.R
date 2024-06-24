# This is an R package containing useful functions for my work.

# Available functions:
# 1. Simul_Con_MULT() -- simulates contingency tables based on the multinomial distribution.
# 2. Simul_Con_MULT.FISH.ORD() -- simulates ordinal-distributed data across treatments and lesions with inter-fish variation in the PO.
# 3. Simul_Surv() -- simulate survival data based on a reference hazard function, the specified hazard ratio(s), and inter-tank variation.
# 4. theme_Publication() -- ggplot theme for generating publication-ready plots.
# 5. Predict_SR() -- predict future survival rate(s) for ongoing experiment based on a reference hazard function from older data.
# 6. Surv_Gen0() -- generate rows of survivors given a starting number of fish per tank and data containing morts and sampled fish.

# Some functions require datasets as inputs. To show how these datasets are structured, sample datasets can be loaded by executing
# "data(function.name_argument.name)" in the R console. For example, for function #6, type "data(Surv_Gen0_DB_Mort)" where "DB_Mort"
# is the argument for the input dataset. Once the sample dataset is loaded, the function can run on "nothing" to produce an example output.
# For example, I can type "Surv_Gen0()" in the R console and obtain the example output.

###################################################################################################################################
################################################## Function 1 - Simul_Con_MULT() ##################################################

Simul_Con_MULT = function(total_count = 750,
                          n_lesion = 3,
                          n_Trt. = 5) { #default conditions

  #Generate vector of counts
  Count_Vec1 = rmultinom(n = 1, size = total_count, p = rep(1/(n_lesion * n_Trt.), n_lesion * n_Trt.))

  #Create contingency table
  Con_Tab1 = matrix(Count_Vec1, nrow = n_Trt., ncol = n_lesion)
  dimnames(Con_Tab1) = list(Trt. = LETTERS[1:n_Trt.],
                            LS = 1:n_lesion)

  return(Con_Tab1)
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

Predict_SR = function(New_DB = Predict_SR_New_DB, #Data from ongoing study, with SR to be predicted. Need specific column names, see data(Predict_SR_New_DB)
                      Ref_DB = Predict_SR_Ref_DB, #Reference survival data to create the reference hazard function. See data(Predict_SR_Ref_DB)
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
################################################## Function 6 - Surv_Gen0() #######################################################

Surv_Gen0 = function(DB_Mort = Surv_Gen0_DB_Mort,  #Mort data with specific column names. See data(Surv_Gen0_DB_Mort).
                     Starting_Number_of_Fish_per_Tank = 70) {

  require(dplyr)
  require(plyr)

  DB_Mort_Gensum = data.frame(DB_Mort %>%
                                dplyr::group_by(Trt.ID, Tank.ID) %>%
                                dplyr::summarise(Num_dead = n()))

  DB_Mort_Gensum$Num_alive = Starting_Number_of_Fish - DB_Mort_Gensum$Num_dead
  DB_Mort_Genalive = data.frame(lapply(DB_Mort_Gensum, rep, DB_Mort_Gensum$Num_alive))
  DB_Mort_Genalive$Status = 0
  DB_Mort_Genalive$TTE = max(DB_Mort_Gen$TTE)
  DB_Mort_Gencomb = rbind.fill(DB_Mort_Gen, DB_Mort_Genalive[,-c(3:4)])

  return(DB_Mort_Gencomb)
}
