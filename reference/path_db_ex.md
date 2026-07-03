# Example Pathology Data

A dataframe containing Trt.ID, Tank.ID, and various (13) pathological
signs containing binary data (Y/N). Contains "dirty" entries (e.g. N/AP,
ME).

## Usage

``` r
data(path_db_ex)
View(path_db_ex)
```

## Format

A data frame containing 591 rows and 12 columns:

|  |  |  |
|----|----|----|
| `Trt.ID` |  | Unique label representing the treatment group associated with a row of data |
| `Tank.ID` |  | Unique label representing each different tank |
| `misc...` |  | Various pathological signs (columns 3-12) |
