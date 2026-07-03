# Example Mort Data

An example mortality dataframe that can be accepted by
[`Surv_Gen()`](https://sean4andrew.github.io/safuncs/reference/Surv_Gen.md).

## Usage

``` r
data(mort_db_ex)
View(mort_db_ex)
```

## Format

A data frame containing 399 rows and 4 columns:

|  |  |  |
|----|----|----|
| `Tank.ID` |  | Unique labels for the different tanks in the study |
| `Trt.ID` |  | Unique labels for the different treatments in the study |
| `TTE` |  | Time to Event. In this dataset, TTE = days post challenge |
| `Status` |  | Value indicating what happened at TTE. In this dataset, Status = 1 indicating all events are death |
