
# Your Package Name

This R package provides functions for simulating repeated opponent
processes using the behavioral posology paradigm and for using Bode
plots to analyze the healthy hormetic limits of behaviors.

## Installation

You can install this package from CRAN using the following R command:

``` r
# install.packages("bpos")

## REPLACE WITH GITHUB NAME
```

## Usage

### Simulating Opponent Processes

To simulate opponent processes using the opponentprocess() function, you
can use the following code:

``` r
# library(your_package_name)
# 
# # Simulate opponent processes with default parameters
# opponentprocess()
# 
# # Simulate opponent processes with custom parameters
# opponentprocess(ii = 10, sim_length = 4000, addl = 10000, plot_utility = FALSE)
```

## Creating Bode Plots

You can create Bode plots using the bode_plot() function as follows:

``` r
# # library(your_package_name)
# 
# # Create a Bode plot with default parameters
# bode_plot()
# 
# # Create a Bode plot with custom parameters
# bode_plot(freq_interval = 0.0002, multiply = 40)
```

## Function Details

### opponentprocess()

The opponentprocess() function simulates opponent processes. It takes
various parameters, including ii (dosing interval), sim_length (time
length of simulation), and others. You can customize the simulation by
adjusting these parameters.

### bode_plot()

The bode_plot() function creates a Bode plot using the opponentprocess()
function for a range of dose frequencies. You can specify the
freq_interval and multiply parameters to control the range of
frequencies to plot.

For more details on function parameters and usage, please refer to the
package documentation by using the ‘help()’ function in R.

For more usage examples and detailed explanations, please check the
package vignettes and documentation.

## License

This package is released under the MIT License.

<!-- # Testing out JUST PK alterations. Can STILL ACHIEVE HORMESIS?? -->
<!-- ```{r} -->
<!-- #| fig-width: 9 -->
<!-- #| fig-height: 12 -->
<!-- #| warning: false -->
<!-- #| results: 'hide' -->
<!-- #| fig-dpi: 600 -->
<!-- bode_plot(EC50_b = c(0.8, 0.9, 1, 1.1), k_bpk = 0.02, seq_2 = 0.08, plot_2 = c(0.0015, 0.08)) -->
<!-- ``` -->
<!-- # Work on a combined figure leading up to final hormetic eigenvector -->
<!-- ```{r} -->
<!-- #| fig-width: 9 -->
<!-- #| fig-height: 12 -->
<!-- #| warning: false -->
<!-- #| results: 'hide' -->
<!-- #| fig-dpi: 300 -->
<!-- bode_plot(EC50_b=c(2.1, 2.8, 3.5, 4.2), seq_1=0.00025, sim_length=4000, gg_ylim=-500, meta_curve=TRUE) -->
<!-- ``` -->
<!-- # Repeat, but show that modifying a-process does little to shift hormetic curve -->
<!-- ```{r} -->
<!-- #| fig-width: 9 -->
<!-- #| fig-height: 12 -->
<!-- #| warning: false -->
<!-- #| results: 'hide' -->
<!-- #| fig-dpi: 300 -->
<!-- bode_plot(EC50_a=c(0.6, 0.8, 1.1, 1.5), seq_1=0.00025, sim_length=4000, gg_ylim=-500, meta_curve=TRUE) -->
<!-- ``` -->
