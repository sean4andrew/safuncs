# Generates Cleaned Pathology Data

Prepares a standardized pathology dataframe for further analysis using
follow-up functions in the safuncs package (e.g. `Path_Sum()`). Filters
for relevant columns and standardizes contents of pathology columns by
converting values to either NA, "Present", or "Absent".

## Usage

``` r
Path_Gen(
  path_db,
  rel_cols,
  path_cols,
  to_na = c(" ", "N/AP", "", "N/R", "N/A", "N/Ap", "N/ap", "N/Av", "ME", ""),
  to_present = c("Y", "Yes", "yes", "y"),
  to_absent = c("N", "No", "no", "n")
)
```

## Arguments

- path_db:

  A dataframe containing at least one column for the relevant factor
  (e.g. Trt.ID) and other columns for pathological signs (one for each
  sign).

- rel_cols:

  A numeric vector representing the column indices of non-pathology
  columns that are to be carried forward to the output/return of this
  function.

- path_cols:

  A numeric vector specifying the column indices of pathology columns
  which values are to be standardized using the "to\_" arguments of this
  function.

- to_na:

  A character vector specifying what values in `path_cols` are to be
  converted to NA which represents missing data.

- to_present:

  A character vector specifying what values in `path_cols` are to be
  converted to "Present".

- to_absent:

  A character vector specifying what values in `path_cols` are to be
  converted to "Absent".

## Value

A dataframe containing only the relevant columns – all of `rel_cols` and
`path_cols`. Values in the pathology columns standardized according to
function "to\_" arguments.

## Examples
