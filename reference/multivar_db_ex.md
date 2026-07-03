# Example Multivariate Data

An example dataset representing the mucus chemistry of fish exposed to
different treatments. The dataframe can be inputed into
[`MultiVar()`](https://sean4andrew.github.io/safuncs/reference/MultiVar.md).
Mucus chemistry is described by several variables from columns 2 to 8.
In column 9, a fake factor "Fruits" was added for use as an example
second factor. In column 4, row 10, the available value was replaced by
NA to use as an example missing value.

## Usage

``` r
data(multivar_db_ex)
View(multivar_db_ex)
```

## Format

An object of class `data.frame` with 48 rows and 9 columns.
