# Generates Cleaned Pathology Data

Prepares a standardized pathology dataframe for further analysis using
follow-up functions in the safuncs package (e.g.
[`Path_Summary()`](https://sean4andrew.github.io/safuncs/reference/Path_Summary.md)).
Filters for relevant columns and standardizes contents of pathology
columns by converting values to either NA, "Present", or "Absent".

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

``` r
#Load example pathology dataframe and view unique values per column
data(path_db_ex)
lapply(path_db_ex, unique)
#> $Trt.ID
#> [1] "D" "A" "C" "B" "E"
#> 
#> $Tank.ID
#>  [1] "E21" "E18" "E08" "E12" "E16" "E17" "E31" "E10" "E22" "E20" "E04" "E05"
#> [13] "E06" "E07" "E03" "E02" "E13" "E15" "E23" "E24" "E26" "E28" "E35" "E38"
#> [25] "E39" "E40" "E41" "E42" "E43" "E44" "E01" "E09" "E11" "E14" "E27" "E30"
#> [37] "E32" "E37" "E19" "E36"
#> 
#> $Dark.body
#> [1] "N"    "N/AP" "Y"   
#> 
#> $Exopthlamia
#> [1] "N"    "Y"    "N/AP"
#> 
#> $Mouth.Haemorrhage
#> [1] "Y"    "N"    "N/AP"
#> 
#> $Skin.Haemorrhage
#> [1] "Y"    "N"    "N/AP"
#> 
#> $Eye.Haemorrhage
#> [1] "N"    "Y"    "N/AP"
#> 
#> $Ascites
#> [1] "N"    "N/AP" "Y"   
#> 
#> $Swollen.Spleen
#> [1] "N"    "Y"    "N/AP"
#> 
#> $Visceral.Haemorrhage
#> [1] "Y"    "N"    "N/AP"
#> 
#> $Grey.Kidney
#> [1] "Y"    "N"    "N/AP"
#> 
#> $Swollen.Kidney
#> [1] "Y"    "N/AP" "N"   
#> 

#Prepare dataset for further analyses
path_db = Path_Gen(path_db = path_db_ex,
                   rel_cols = 1,
                   path_cols = 3:12)
#Again, view the unique values per column, note "N/AP" values converted to "NA",
#"Y" converted to "Present", and "N" converted to "Absent".
lapply(path_db, unique)
#> $Trt.ID
#> [1] "D" "A" "C" "B" "E"
#> 
#> $Dark.body
#> [1] "Absent"  "Present"
#> 
#> $Exopthlamia
#> [1] "Absent"  "Present"
#> 
#> $Mouth.Haemorrhage
#> [1] "Present" "Absent" 
#> 
#> $Skin.Haemorrhage
#> [1] "Present" "Absent" 
#> 
#> $Eye.Haemorrhage
#> [1] "Absent"  "Present"
#> 
#> $Ascites
#> [1] "Absent"  NA        "Present"
#> 
#> $Swollen.Spleen
#> [1] "Absent"  "Present" NA       
#> 
#> $Visceral.Haemorrhage
#> [1] "Present" "Absent"  NA       
#> 
#> $Grey.Kidney
#> [1] "Present" "Absent"  NA       
#> 
#> $Swollen.Kidney
#> [1] "Present" NA        "Absent" 
#> 
```
