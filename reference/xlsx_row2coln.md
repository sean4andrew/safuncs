# Set First Row as Column Headers

Transform first row of the data into column names and subsequently
remove that row.

## Usage

``` r
xlsx_row2coln(x)
```

## Arguments

- x:

  A tibble / dataframe object. Initially designed for the output of
  [`readxl::read_xlsx()`](https://readxl.tidyverse.org/reference/read_excel.html).

## Value

Returns a dataframe where the previous first rows are now column names
