# Generate Survivor Data

Produces survival data that includes rows for every surviving fish based
on the starting number of fish and mortality data. To generate survivor
data for tanks absent in the input mortality dataframe, specify the
arguments `tank_without_mort` and `trt_without_mort`. To generate
survivor data with tank specific starting numbers of fish, input a
dataframe into the argument `starting_fish_count` instead of a single
value; details in **Arguments**.

## Usage

``` r
Surv_Gen(
  mort_db,
  mort_fish_count = FALSE,
  starting_fish_count,
  last_tte,
  add_factor = NULL,
  add_sampled = NULL,
  tank_without_mort = NULL,
  trt_without_mort = NULL,
  output_prism = FALSE,
  output_prism_date = NULL
)
```

## Arguments

- mort_db:

  A mort dataframe as described in **Details**.

- mort_fish_count:

  Whether to replicate data (by row) based on the number of fish
  specified under the column 'fish_count'.

- starting_fish_count:

  Value representing the starting number of fish for every tank.
  Alternatively, a dataframe containing the columns "Trt.ID", "Tank.ID",
  and "starting_fish_count" to allow for different fish starting numbers
  per tank.

- last_tte:

  Value representing the time-to-event the fish survived to, assigned to
  every row of survivor data generated.

- add_factor:

  A string or character vector representing the name(s) of column(s) in
  `mort_db` to be carried over to the generated survival data for
  further analysis (e.g. as facet factor in
  [`safuncs::Surv_Plots()`](https://sean4andrew.github.io/safuncs/reference/Surv_Plots.md)).
  Column must also be present in `starting_fish_count` dataframe.
  Defaults to NULL.

- add_sampled:

  A dataframe containing the column names "sampled_per_tank" and
  "sampled_tte" to indicate the amounts and times sampled. Each row of
  the dataframe is correlated (i.e. a specific time for specific
  sampling per tank). Defaults to NULL.

- tank_without_mort:

  A vector of strings specifying the tanks absent from `mort_db`; used
  to generate survivor data for those tanks. Argument ignored if
  `starting_fish_count` is a dataframe.

- trt_without_mort:

  A vector of strings corresponding to `tank_without_mort`. Keep their
  order the same. Argument ignored if `starting_fish_count` is a
  dataframe.

- output_prism:

  Whether to generate and save a prism ready survival csv. Defaults to
  FALSE.

- output_prism_date:

  The starting date to be used in the prism file. Please specify date in
  "dd-Mmm-yyyy" syntax (e.g. "08-Aug-2024").

## Value

A dataframe produced by combining the input mort data and generated rows
of survivor data.

## Details

The mort dataframe supplied as input should consist of the following 4
columns at minimum:

- "Trt.ID" = Labels for treatment groups in the study.

- "Tank.ID" = Labels for tanks in the study (each tank must have a
  unique label).

- "TTE" = Time to Event. Event could be fish death or being sampled and
  removed depending on "Status".

- "Status" = Value indicating what happened at TTE. 1 for dead fish, 0
  for those sampled and removed.

Each row should represent one fish.

For an example dataframe, execute `data(mort_db_ex)` and view.

## See also

[Link](https://sean4andrew.github.io/safuncs/reference/Surv_Gen.html)
for executed **Examples** which includes any figure outputs.

## Examples

``` r
# First, we load an example mortality database available from the safuncs package
data(mort_db_ex)

# Next, we input this data into Surv_Gen() as well as the study details to generate
# entries (rows) for survivors in the output - a "complete" dataframe for further
# survival analysis and data visualization.
Surv_Data_Output = Surv_Gen(mort_db = mort_db_ex,
                            starting_fish_count = 100,
                            last_tte = 54,
                            tank_without_mort = c("C99", "C100"),
                            trt_without_mort = c("A", "B"))
#> [1] "Your total number of tanks is: 14"
#> [1] "Your total number of treatment groups is: 4"
#> [1] "Your total number of fish in the output data is: 1400"

# Below, the bottom 5 rows of the output is displayed to show the rows of survivor
# data generated.
tail(Surv_Data_Output, n = 5)
#>      Tank.ID Trt.ID TTE Status
#> 1396    C100      B  54      0
#> 1397    C100      B  54      0
#> 1398    C100      B  54      0
#> 1399    C100      B  54      0
#> 1400    C100      B  54      0

# Below is another example, this time showing how to specify tank-specific fish
# numbers. First, create the database containing information on the starting fish
# counts for the different tanks. You must include all tanks that are in your mort
# database to get a proper output. For the example below, we will later trim down
# the mort database to only 4 tanks for simplicity. Use these 3 column names:
# starting_fish_count, Tank.ID and Trt.ID.
count_db = data.frame(starting_fish_count = c(100, 100, 120, 120),
                      Tank.ID = c("C1", "C6", "C5", "C8"),
                      Trt.ID = c("B", "B", "D", "D"))

filtered_mort_db = mort_db_ex[mort_db_ex$Tank.ID %in% c("C1", "C6", "C5", "C8"),]

# We then use 'count_db' as input to the argument 'starting_fish_count' in Surv_Gen():
Surv_Data_Output = Surv_Gen(mort_db = filtered_mort_db,
                            starting_fish_count = count_db,
                            last_tte = 54)
#> [1] "Your total number of tanks is: 4"
#> [1] "Your total number of treatment groups is: 2"
#> [1] "Your total number of fish in the output data is: 440"

# Surv_Gen() is also able to generate data for sampled fish (Status = 0) at specified
# TTEs, using the argument add_sampled:
Surv_Data_Output = Surv_Gen(mort_db = filtered_mort_db,
                            starting_fish_count = count_db,
                            last_tte = 54,
                            add_sampled = data.frame(sampled_per_tank = 5,
                                                     sampled_tte = 30))
#> [1] "Your total number of tanks is: 4"
#> [1] "Your total number of treatment groups is: 2"
#> [1] "Your total number of fish in the output data is: 440"
```
