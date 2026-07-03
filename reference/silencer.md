# Silence Code Output

Hide output from R console by redirecting output using
`sink(tempfile())` and subsequently
[`sink()`](https://rdrr.io/r/base/sink.html).

## Usage

``` r
silencer(x)
```

## Arguments

- x:

  Code which output is to be directed to the sink.

## Value

Code output without the directed outputs, e.g. output from
[`cat()`](https://rdrr.io/r/base/cat.html).
