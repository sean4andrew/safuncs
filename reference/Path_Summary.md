# Create Summary Statistics from Pathology Data

Takes a cleaned pathology dataframe and summarizes, for each
pathological sign, its prevalence, standard error, and n for each unique
combination of levels of factor(s). Additionally, returns p-values and
letter groups based on logistic regression and likelihood ratio tests,
with p-value adjustments to account for multiple pairwise comparisons;
the logistic regression model and letters considers only the first
specified factor and ignores all others.

## Usage

``` r
Path_Summary(
  path_db,
  path_cols,
  factors = c("Trt.ID"),
  p_adj = "BH",
  contrast_out = TRUE
)
```

## Arguments

- path_db:

  A cleaned pathology dataframe such as that from
  [`Path_Gen()`](https://sean4andrew.github.io/safuncs/reference/Path_Gen.md).
  Must contain at least one column as the factor, and one column for a
  pathological sign. Multiple signs path

- path_cols:

  A numeric vector representing the column indices of pathological
  signs.

- factors:

  A character vector specifying the factor(s) relevant in the dataframe.
  Defaults to "Trt.ID", Only the first factor is considered in
  statistical tests, however, summary statistics such as prevalence,
  standard errors, and n, are calculated for every possible combination
  of factor level.

- p_adj:

  A string specifying the p-value adjustment method to account for
  multiple comparisons. Defaults to "BH". Method based on
  [`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

- contrast_out:

  A logical (TRUE/FALSE) specifying whether to return an additional
  dataframe containing the p-value based on pairwise comparison of
  groups/level of the first specified factor in `factors`. Defaults to
  TRUE.

## Value

A dataframe or list (when `contrast_out = TRUE`) containing summary
statistics and results from univariate hypothesis test.

## Examples

``` r
#Generate pathology dataset:
data(path_db_ex)
path_db = Path_Gen(path_db = path_db_ex,
                   rel_cols = 1,
                   path_cols = 3:12)

path_sum = Path_Summary(path_db = path_db,
                        factor = c("Trt.ID"),
                        path_cols = 2:11,
                        contrast_out = TRUE)
#> [1] "Ascites"
#> [1] "Dark body"
#> [1] "Exopthlamia"
#> [1] "Eye haemorrhage"
#> [1] "Grey kidney"
#> [1] "Mouth haemorrhage"
#> [1] "Skin haemorrhage"
#> [1] "Swollen kidney"
#> [1] "Swollen spleen"
#> [1] "Visceral haemorrhage"
#Dataframe containing summary statistics and hypothesis test across factor levels of the first factor (in this case "Trt.ID"):
head(path_sum$summary_db, 5)
#>      Path Trt.ID Presence  n      prop sum_n    prop_se  pos   logis_p Letters
#> 1 Ascites      A   Absent 24 0.8888889    27 0.06048123 0.95 0.3718439       a
#> 2 Ascites      A  Present  3 0.1111111    27 0.06048123 0.05 0.3718439       a
#> 3 Ascites      B   Absent 20 0.7692308    26 0.08262864 0.95 0.3718439       a
#> 4 Ascites      B  Present  6 0.2307692    26 0.08262864 0.05 0.3718439       a
#> 5 Ascites      C   Absent 19 0.7600000    25 0.08541663 0.95 0.3718439       a

#Dataframe containing P-values from pairwise comparison of levels of the first factor
head(path_sum$contrast_db, 5)
#>      Path Contrast      P
#> 1 Ascites      A-B 0.6068
#> 2 Ascites      A-C 0.6068
#> 3 Ascites      A-D 0.7174
#> 4 Ascites      A-E 0.7174
#> 5 Ascites      B-C 0.9381
```
