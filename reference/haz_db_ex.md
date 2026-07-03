# Example Hazard Data

A reference hazard dataframe created using
`Surv_Plots(..., dailybin = TRUE, data_out = TRUE)$Hazard_DB` which uses
[`bshazard::bshazard()`](https://rdrr.io/pkg/bshazard/man/bshazard.html)
inside. Contains hazard rates over time.

## Usage

``` r
data(haz_db_ex)
View(haz_db_ex)
```

## Format

A data frame containing 54 rows and 3 columns:

|  |  |  |
|----|----|----|
| `Trt.ID` |  | Unique labels representing the different treatment groups used in creating this reference hazard dataframe |
| `Hazard` |  | Hazard values (rates) |
| `Time` |  | Time / TTE in days post challenge. |
