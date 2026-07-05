# Create Summary Statistics from Pathology Data

Takes a cleaned pathology dataframe and summarizes, for each
pathological sign, its prevalence, standard error, and n across each
unique combination of factors. P-value and letter groups from logistic
regression and a wald-type test with p-value adjustments are also
produced; the logistic regression model and letters considers only the
first specified factor and ignores all others.

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
#placeholder
```
