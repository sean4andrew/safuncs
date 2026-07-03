# Get Default Colours by ggplot

Not my function but it is useful so here it is!
[Origin](https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette)

## Usage

``` r
GG_Colour_Hue(n)
```

## Arguments

- n:

  Number of colour groups

## Value

Returns a vector representing the default colour codes assigned to each
group by ggplot.

## See also

[Link](https://sean4andrew.github.io/safuncs/reference/GG_Colour_Hue.html)
for executed **Examples** which includes any figure outputs.

## Examples

``` r
# Get colour codes used for 6 categorical groups
GG_Colour_Hue(6)
#> [1] "#F8766D" "#B79F00" "#00BA38" "#00BFC4" "#619CFF" "#F564E3"
```
