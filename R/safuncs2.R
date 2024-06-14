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

Simul_SURV = function(Cumuhaz_DB = Cumuhaz_Base_DB, #object from bshazard()
                      fish_num_per_tank = 100,
                      tank_num_per_trt = 4,
                      Treatments_HR = c(1, 1, 1, 1), #HR for treatment groups, starting from control (1)
                      logHR_sd_intertank = 0.35, #inter-tank variation (sd) in HR in the log-scale
                      round_digits_days = 0) { #Time is rounded to Days (i.e. fish cannot die at Day 1.4 or any other decimal)
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
  TTE = approx(x = cumsum(Cumuhaz_DB$hazard), y = Cumuhaz_DB$time, xout = CDF_Yval, method = "linear")$y
  TTE = round(TTE, digits = round_digits_days)

  #Turn NA (from out of bound CDF_Yval) to the last follow up time
  TTE = ifelse(is.na(TTE), max(Cumuhaz_DB$time), TTE)

  #Label Status (1 - dead, or 0 - survived) given TTE
  Surv_simul_DB = data.frame(TimetoEvent = TTE,
                             CDF_Yval = CDF_Yval,
                             Status = ifelse(TTE == max(Cumuhaz_DB$time), 0, 1),
                             Treatment_G = as.factor(rep(LETTERS[1:length(Treatments_HR)], each = fish_num_per_tank * tank_num_per_trt)),
                             Tank_num = as.factor(rep(1:(length(Treatments_HR) * tank_num_per_trt), each = fish_num_per_tank)))

  return(Surv_simul_DB)
}



theme_Publication = function(base_size = 14, base_family = "helvetica") {
  library(grid)
  library(ggthemes)
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
