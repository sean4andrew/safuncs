# Simulate Contingency Table

Simulate a contingency table with fish counts distributed across *n*
lesion categories and *n* treatment groups. Probability values for
generating counts in each cell (i.e. each factor level combination) can
be assigned using the `probs` argument. This function is designed for
use in power and/or false positive rate calculations; for details, see
[`Con_Simul_PR()`](https://sean4andrew.github.io/safuncs/reference/Con_Simul_PR.md).

## Usage

``` r
Con_Simul(
  probs = "equal",
  total_count = 750,
  n_lesion = 3,
  n_Trt. = 5,
  margin_fixed_Trt. = FALSE,
  verbose = TRUE
)
```

## Arguments

- probs:

  Matrix of probability values created using
  [`matrix()`](https://rdrr.io/r/base/matrix.html). Each row in the
  matrix should represent a treatment group and each column a lesion
  category. All probability values in the matrix should sum to 1.
  Default = equal probability across all cells.

- total_count:

  Total number of counts in the contingency table. Defaults to 750.

- n_lesion:

  Number of lesion categories. Ignored if `probs` specified. Defaults to
  3.

- n_Trt.:

  Number of treatment groups. Ignored if `probs` specified. Defaults to
  5.

- margin_fixed_Trt.:

  Whether margins are fixed per treatment group (i.e. fixed number of
  fish per treatment). Default = FALSE. See **Details** for further
  information on marginals.

- verbose:

  Whether to print the parameters and probability matrix used. Default =
  TRUE.

## Value

Returns a list containing:

|  |  |  |
|----|----|----|
| `sim_tab` |  | The simulated contingency table containing counts across different treatment groups (rows) and lesion categories (columns) |
| `params` |  | The simulation parameters as a vector |
| `probs` |  | The probability matrix used for simulation |

## Details

Counts are simulated from a multinomial distribution using
[`rmultinom()`](https://rdrr.io/r/stats/Multinom.html). Counts may be
assumed to have a fixed total in the marginals (e.g. per treatment
group) or no fixed total in row or column marginals.

For further discussion into the types of marginals in contingency
tables, refer to:
[here](https://www.uvm.edu/~statdhtx/StatPages/More_Stuff/Chi-square/Contingency-Tables.pdf)
and the comments on **Arguments**.

## See also

[Link](https://sean4andrew.github.io/safuncs/reference/Con_Simul.html)
for executed **Examples** which includes any figure outputs.

## Examples

``` r
# Simulate table with uniform probabilities across cells
Con_Simul(total_count = 750, n_lesion = 3, n_Trt. = 5)
#> $sim_tab
#>     Lesion_Category
#> Trt.  1  2  3
#>    A 43 59 46
#>    B 48 60 49
#>    C 45 49 51
#>    D 44 50 45
#>    E 48 55 58
#> 
#> $params
#>       total_count          n_lesion            n_Trt. margin_fixed_Trt. 
#>               750                 3                 5                 0 
#> 
#> $probs
#>            [,1]       [,2]       [,3]
#> [1,] 0.06666667 0.06666667 0.06666667
#> [2,] 0.06666667 0.06666667 0.06666667
#> [3,] 0.06666667 0.06666667 0.06666667
#> [4,] 0.06666667 0.06666667 0.06666667
#> [5,] 0.06666667 0.06666667 0.06666667
#> 

# Simulate table with specified probabilities across cells
Con_Simul(probs = matrix(nrow = 2, ncol = 3, c(1/6, 3/12, 1/6, 1/6, 1/6, 1/12)))
#> $sim_tab
#>     Lesion_Category
#> Trt.   1   2   3
#>    A 127 117 112
#>    B 174 154  66
#> 
#> $params
#>       total_count          n_lesion            n_Trt. margin_fixed_Trt. 
#>               750                 3                 2                 0 
#> 
#> $probs
#>           [,1]      [,2]       [,3]
#> [1,] 0.1666667 0.1666667 0.16666667
#> [2,] 0.2500000 0.1666667 0.08333333
#> 
```
