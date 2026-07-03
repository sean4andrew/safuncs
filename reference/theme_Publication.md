# Publication theme for ggplot2

A theme function to add to a ggplot2 object for publication style plots.
Function adapted from
[HanjoStudy/quotidieR](https://rdrr.io/github/HanjoStudy/quotidieR/src/R/theme_publication.R).

## Usage

``` r
theme_Publication(base_size = 14)
```

## Arguments

- base_size:

  size of text in graph

## Value

theme function to add to a ggplot2 object

## See also

[Link](file:///C:/Users/sean4/Documents/GitHub/safuncs/docs/reference/theme_Publication.md)
for executed **Examples** which includes any figure outputs.

## Examples

``` r
# Load an example dataset
data(iris)

# Create a ggplot modified with theme_Publication()
library(ggplot2)
ggplot(data = iris, aes(x = Species, colour = Species, y = Petal.Length)) +
   geom_boxplot() +
   theme_Publication()
```
