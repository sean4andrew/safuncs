# Create Flextables from Prevalence Data

Converts the pathology data summary from
[`Path_Summary()`](https://sean4andrew.github.io/safuncs/reference/Path_Summary.md)
into a list of flextable objects to print out to Word. Summary data are
converted to flextables in wide and long formats. Contrast data
(pairwise comparison p-values) are converted to their standalone
flextable.

## Usage

``` r
Path_Table(path_sum, digits = 0, prev_append = c("se_plus-minus"))
```

## Arguments

- path_sum:

  A list of dataframes or a single dataframe from
  [`Path_Summary()`](https://sean4andrew.github.io/safuncs/reference/Path_Summary.md).

- digits:

  A number specifying the decimal places for the prevalence values in
  the output flextable. Defaults to 0.

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
#Generate summarized pathology dataframe:
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

#Generate pathology table:
path_tables = Path_Table(path_sum = path_sum,
                         digits = 1, #1 decimal places for prevalence values
                         prev_append = c("se_brackets")) #append standard errors (in brackets) to prevalence values

#View pathology table (wide format):
path_tables$wide


.cl-2568e5ea{table-layout:auto;width:100%;}.cl-25610c3a{font-family:'DejaVu Sans';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-25610c4e{font-family:'DejaVu Sans';font-size:6.6pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;vertical-align:super;}.cl-256461aa{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-256486da{background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-256486e4{background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-256486ee{background-color:transparent;vertical-align: middle;border-bottom: 1pt dashed rgba(211, 211, 211, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-256486f8{background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1pt dashed rgba(211, 211, 211, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-256486f9{background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(255, 255, 255, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}

Table 1. Pathological sign prevalence (%) and standard errors (in brackets) across treatments.

Trt.ID
```
