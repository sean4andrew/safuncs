# Create Flextables from Prevalence Data

Converts the dataframe(s) from
[`Path_Summary()`](https://sean4andrew.github.io/safuncs/reference/Path_Summary.md)
into flextable objects in wide and long formats.

## Usage

``` r
Path_Table(path_sum, digits = 0, prev_append = c("se_plus-minus"))
```

## Arguments

- path_sum:

  A list or dataframe from
  [`Path_Summary()`](https://sean4andrew.github.io/safuncs/reference/Path_Summary.md).

- digits:

  A number specifying the decimal places for the prevalence values in
  the output flextable.

- prev_append:

  A string specifying the suffix to append to prevalence values in the
  output table. Suffix may be the standard error or counts of presence.
  Suffix may be presented inside brackets or after a +/- sign. Choose
  one of the following arguments: "se_brackets", "se_plus-minus",
  "presence_brackets". Defaults to "se_plus-minus".

## Value

A list that contains, at minimum, 2 flextable objects of the summary
statistics (from
[`Path_Summary()`](https://sean4andrew.github.io/safuncs/reference/Path_Summary.md))
in wide and long format. If a contrast dataframe from
[`Path_Summary()`](https://sean4andrew.github.io/safuncs/reference/Path_Summary.md)
is also included in the `path_sum` argument, then the output list
contains additionally a flextable for the contrast dataframe.

## Examples

``` r
#placeholder
```
