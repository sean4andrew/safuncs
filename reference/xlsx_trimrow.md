# Trim Rows Based on Non-NA Values

Remove rows after the last non-NA value in a selected column. Select
column based on the `coli` argument.

## Usage

``` r
xlsx_trimrow(x, coli = 1)
```

## Arguments

- x:

  A dataframe.

- coli:

  A number indicating the index of the column to base the trimming on.

## Value

Returns a dataframe object without the "extra" NA values on the selected
rows.
