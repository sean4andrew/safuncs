# Create Flextables from Prevalence Data

Converts the dataframe(s) from
[`Path_Summary()`](https://sean4andrew.github.io/safuncs/reference/Path_Summary.md)
into flextable objects in wide and long formats. Flextable can then be
printed to word using `officer` package.

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

path_tables = Path_Table(path_sum = path_sum,
                         digits = 1, #1 decimal places for prevalence values
                         prev_append = c("se_brackets")) #append standard errors (in brackets) to prevalence values
#> Error in loadNamespace(x): there is no package called ‘ftExtra’

#View pathology table (wide format):
path_tables$wide
#> Error: object 'path_tables' not found

#View pathology table (long format):
path_tables$long
#> Error: object 'path_tables' not found

#View pathology contrast table:
path_tables$contrast
#> Error: object 'path_tables' not found

#Print results in word docx:
library(officer)
library(flextable)
results_doc = read_docx() |>
  body_add_flextable(value = path_tables$wide, align = "left", topcaption = TRUE) |>
  body_add_break() |>
  body_add_flextable(value = path_tables$long, align = "left", topcaption = TRUE) |>
  body_add_break() |>
  body_add_flextable(value = path_tables$contrast, align = "left", topcaption = TRUE)
#> Error: object 'path_tables' not found
print(results_doc, target = paste0("ONDA Pathology Analyses ", format(Sys.Date(), "%d%b%Y"),".docx"))
#> Error: object 'results_doc' not found
```
