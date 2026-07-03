# Calculate Positive Rates for Contingency Table

Computes statistical power and optionally false positive rates for tests
applied to contingency tables based on simulations. Specify the
simulation process using
[`Con_Simul()`](https://sean4andrew.github.io/safuncs/reference/Con_Simul.md),
which serves as input. Positive rates are computed for the Chi-square
test and optionally for Fisher's exact test and the Wald test applied to
an ordinal regression model.

## Usage

``` r
Con_Simul_PR(
  Con_Simul_Object = Con_Simul(),
  add_fisher_exact = FALSE,
  add_ord = FALSE,
  sample_sizes = NULL,
  n_sim = 1000,
  FPR = TRUE
)
```

## Arguments

- Con_Simul_Object:

  Output from
  [`Con_Simul()`](https://sean4andrew.github.io/safuncs/reference/Con_Simul.md)
  with the argument `verbose` set to TRUE.

- add_fisher_exact:

  Whether to compute positive rates for Fisher's Exact test. May add \>1
  min of calculation time. Defaults to FALSE.

- add_ord:

  Whether to compute positive rates for Wald test on a fitted ordinal
  regression model. May add \>1 min of calculation time. Defaults to
  FALSE.

- sample_sizes:

  A vector of sample sizes over which false positive rates are to be
  calculated. A sample size is defined as the total number of counts in
  a contingency table. Defaults to total count received by
  `Con_Simul_Object`.

- n_sim:

  Number of contingency tables simulated for each positive rate
  calculation. Defaults to 1000.

- FPR:

  Whether to calculate false positive rate in addition to power.
  Defaults to TRUE.

## Value

Returns a list containing:

|  |  |  |
|----|----|----|
| `pos_rate` |  | A dataframe containing positive rates of various tests at specified sample sizes |
| `eff_mat` |  | The probability matrix used to calculate power |
| `null_mat` |  | The probability matrix used to calculate false positive rate |
| `plot` |  | A plot of positive rates of various tests at specified sample sizes |

## Details

Power is defined as the percentage of tests yielding positive results
(p-value \< 0.05) on a set of contingency tables simulated based on the
data-generating process and probability matrix specified in
[`Con_Simul()`](https://sean4andrew.github.io/safuncs/reference/Con_Simul.md).
The specified probability matrix should represent the parameters of the
population where there is a desire to detect a significant effect in the
sample. The simulated contingency tables then reflect the sample
outcomes of the specified population parameter.

False Positive Rate is defined as the percentage of tests yielding
positive results (p-value \< 0.05) on contingency tables simulated based
on a probability matrix without treatment effect. The odd ratios between
lesion scores for all treatment groups are assumed to be the same as
that for control (first row in the probability matrix).

P-values for the Chi-square test are computed using
[`stats::chisq.test()`](https://rdrr.io/r/stats/chisq.test.html) with
default parameters. P-values for Fisher's exact test are computed using
[`stats::fisher.test()`](https://rdrr.io/r/stats/fisher.test.html) with
`simulate.p.value` set to TRUE, alongside default parameters. P-values
for ordinal regression are computed from
[`stats::anova()`](https://rdrr.io/r/stats/anova.html) applied to the
output of [`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html).

## See also

[Link](https://sean4andrew.github.io/safuncs/reference/Con_Simul_PR.html)
for executed **Examples** which includes any figure outputs.

## Examples

``` r
# Below I show how we can perform a simple power calculation using this tool.
# Suppose I want to calculate power for Treatment B which halves the lesions in
# category 2 and 3. I then specify the following probability matrix and feed it into
# Con_Simul():
probs_mat = matrix(nrow = 2, ncol = 3, data = c(1/6, 1/3, 1/6, 1/12, 1/6, 1/12))
probs_mat
#>           [,1]       [,2]       [,3]
#> [1,] 0.1666667 0.16666667 0.16666667
#> [2,] 0.3333333 0.08333333 0.08333333

sim_tab = Con_Simul(probs_mat)

# Next, I feed the output of Con_Simul() into Con_Simul_PR():
Con_Simul_PR(sim_tab,
             sample_sizes = c(50, 100, 150),
             add_ord = TRUE,
             add_fisher_exact = TRUE)
#> $pos_rate
#>    Total_Count Percent_of_Significant_Results  Statistical_Test Class
#> 1           50                           54.3         Chisquare Power
#> 2          100                           88.0         Chisquare Power
#> 3          150                           96.8         Chisquare Power
#> 4           50                           59.0 OrdinalRegression Power
#> 5          100                           90.1 OrdinalRegression Power
#> 6          150                           96.8 OrdinalRegression Power
#> 7           50                           52.4      FishersExact Power
#> 8          100                           87.1      FishersExact Power
#> 9          150                           96.5      FishersExact Power
#> 10          50                            6.6         Chisquare   FPR
#> 11         100                            5.1         Chisquare   FPR
#> 12         150                            5.4         Chisquare   FPR
#> 13          50                            6.0 OrdinalRegression   FPR
#> 14         100                            5.0 OrdinalRegression   FPR
#> 15         150                            5.2 OrdinalRegression   FPR
#> 16          50                            6.2      FishersExact   FPR
#> 17         100                            5.2      FishersExact   FPR
#> 18         150                            5.3      FishersExact   FPR
#> 
#> $eff_mat
#>           [,1]       [,2]       [,3]
#> [1,] 0.1666667 0.16666667 0.16666667
#> [2,] 0.3333333 0.08333333 0.08333333
#> 
#> $null_mat
#>           [,1]      [,2]      [,3]
#> [1,] 0.1666667 0.1666667 0.1666667
#> [2,] 0.1666667 0.1666667 0.1666667
#> 
#> $plot

#> 
# Results: Power is ~55, 86, and 97% for the Chi-square test using total counts of
# 50, 100, and 150, respectively.

# The same power for Chi-square test can be calculated using Cohen's omega (w) method
# which is faster but has its own limitations; e.g. assumes one data generating
# process for the contingency table (the no fixed marginals).
library(pwr)
pwr::pwr.chisq.test(w = ES.w2(probs_mat), df = 2, sig.level = 0.05, N = 100)
#> 
#>      Chi squared power calculation 
#> 
#>               w = 0.3333333
#>               N = 100
#>              df = 2
#>       sig.level = 0.05
#>           power = 0.8563144
#> 
#> NOTE: N is the number of observations
#> 
# Results: Power is 85.6% for the Chi-square test at the total count of 100.
```
