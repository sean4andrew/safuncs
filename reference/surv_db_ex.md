# Example Survival Data

A complete survival dataset with survivors, based of `mort_db_ex`. The
dataframe can be inputed into
[`Surv_Plots()`](https://sean4andrew.github.io/safuncs/reference/Surv_Plots.md).

## Usage

``` r
data(surv_db_ex)
View(surv_db_ex)
```

## Format

A data frame containing 1200 rows and 4 columns:

|  |  |  |
|----|----|----|
| `Tank.ID` |  | Unique labels representing the different tanks in the study |
| `Trt.ID` |  | Unique labels representing the different treatments in the study |
| `TTE` |  | Time to Event. In this dataset, TTE = days post challenge |
| `Status` |  | Value indicating what happened at TTE. In this dataset, Status = 1 or 0 indicating death or survival, respectively |
