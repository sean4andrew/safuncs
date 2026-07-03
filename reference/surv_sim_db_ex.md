# Example Simulated Survival Object

A list generated from
[`Surv_Simul()`](https://sean4andrew.github.io/safuncs/reference/Surv_Simul.md)
containing the following elements: a simulated survival dataset
(`surv_simul_db`), the dataset describing the survival characteristics
of the population (`surv_pop_db`) and a plot illustrating both
(`surv_plots`). The specified simulation parameters (arguments) in
[`Surv_Simul()`](https://sean4andrew.github.io/safuncs/reference/Surv_Simul.md)
are `fish_num_per_tank = 100`, `tank_num_per_trt = 4`,
`treatments_hr = c(1, 0.7, 0.5)`, `logHR_sd_intertank = 0`,
`n_sim = 10`, `plot_out = TRUE`.

## Usage

``` r
data(surv_sim_db_ex)
View(surv_sim_db_ex$surv_simul_db)
View(surv_sim_db_ex$surv_pop_db)
View(surv_sim_db_ex$surv_plots)
```

## Format

An object of class `list` of length 3.
