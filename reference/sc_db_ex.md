# Example Lesion Scores Data

A dataframe of lesion scores that can be accepted by the template script
for scores analyses. Except for *Tank.ID*, the columns described below
are necessary.

## Usage

``` r
data(sc_db_ex)
View(sc_db_ex)
```

## Format

A data frame containing 1012 rows (each representing one fish) and 9
columns:

|  |  |  |
|----|----|----|
| `Timepoint` |  | Unique label representing the timepoint associated with a row of data (i.e. one fish) |
| `Trt.ID` |  | Unique label representing the treatment group associated with a row of data |
| `Tank.ID` |  | Unique label representing each different tank |
| `Mouth` |  | Category (must be 0-4) of severity for the fish based on its mouth lesion(s) |
| `Gill` |  | Category (must be 0-2) of severity for the fish based on its gill lesion(s) |
| `Skin1` |  | Count of skin lesions Category 1 |
| `Skin2` |  | Count of skin lesions Category 2 |
| `Skin3` |  | Count of skin lesions Category 3 |
| `Skin3R` |  | Count of skin lesions Category 3R |
