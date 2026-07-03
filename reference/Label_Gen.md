# Generate Texts for Labels

Combines texts specified in a list. List should contain multiple
variables that holds the texts in a value or vector format. All
combinations of texts across variables are computed. Output text
(combinations) is sorted in order of variables given in the list
(default behavior) or as specified using `sort_by` argument. The output
combinations are tabulated and saved in a .csv in your working
directory.

## Usage

``` r
Label_Gen(
  input_list,
  sort_by = NULL,
  n_col = 6,
  fill_by_row = TRUE,
  save_name = NULL
)
```

## Arguments

- input_list:

  A list of named variables, each containing one or more text/number(s).
  See **Examples** for examples.

- sort_by:

  A value or vector representing the variable(s) to sort the output by.
  For each variable, sorts according to the order of text in the
  variable. When multiple variables is given, prioritizes sorting based
  on the order of variables; leftmost = highest priority. Defaults to
  NULL where sorting is based on `input_list` orders.

- n_col:

  The number of columns in the output table. It should match the number
  of columns in the label paper. Defaults to 6.

- fill_by_row:

  Whether combinations should fill the output table by row (otherwise
  column). Defaults to TRUE.

- save_name:

  Name of the saved .csv. Defaults to NULL where the file name is
  "Label_Gen" and today's date (YYYY-MM-DD).

## Value

A .csv containing all possible combinations. Additionally, a printout
describing the .csv file name and location. Another printout describing
the total number of labels / combinations created.

## See also

[Link](https://sean4andrew.github.io/safuncs/reference/Label_Gen.html)
for web documentation.

## Examples

``` r
# Summarize the input variables in a list
input_variables = list(Time = c("Baseline", "1wpv"),
                       Animal = c("Oysters", "Lobsters"),
                       Tissue = c("Meat", "Shell", "Water", "Head"),
                       Replic_num = 1:3)

# Run Label_Gen() using the input variables.
Label_Gen(input_list = input_variables,
          sort_by = c("Time", "Animal", "Tissue"),
          n_col = 6,
          fill_by_row = TRUE,
          save_name = NULL)
#> [1] "You have 48 total labels"
#> [1] "File saved as Label_Gen 2026-07-03.csv in /home/runner/work/safuncs/safuncs/docs/reference"
```
