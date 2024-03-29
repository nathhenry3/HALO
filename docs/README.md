
# HALO

#### Nathan Henry, 2023

NOTE: This prototype package has been deprecated as it relies on a
minor violation of utility theory. We have since replaced this with
our PK/PD model of opponent processes, which is less prescriptive and
thus more accurate and flexible. However, I have left this here as an
example of my package-coding ability. 

This package allows you to run behavioral posology simulations using the
HALO paradigm (Hormetic Alignment via Opponent Processes), in order to
simulate the healthy limits of behaviors that have opponent process
dynamics.

The updated HALO paradigm combines the concepts of allostasis, hormesis,
opponent processes, and utility theory. For more information, refer to
the [paper by Henry et al.](https://arxiv.org/abs/2402.07462)

## Installation

You can install the HALO package by running the following code in R or
RStudio:

``` r
# Uncomment the line below to install the devtools package if needed
# install.packages('devtools')

# Load devtools package
library(devtools)

# Install and load the HALO package from GitHub
install_github('nathhenry3/HALO')
library(HALO)
```

## Usage

### Performing a BFRA (Behavioral Frequency Response Analysis)

The HALO paradigm allows you to detect the hormetic (healthy) limits of
behaviors, by performing a Behavioral Frequency Response Analysis (BFRA)
of the opponent processes generated by the behaviors.

The main function to use for simulations is bode_plot(). A Bode plot
allows you to see how changing the frequency of a behavior leads to
different hedonic outcomes. You can create Bode plots using the
bode_plot() function as follows:

``` r
bode_plot()
```

![](README_files/figure-commonmark/unnamed-chunk-2-1.png)

The top graph shows the utility function, which describes the
relationship between the pharmacokinetic values (x-axis) and the
pharmacodynamic effects (y-axis).

The pharmacodynamic effects can be observed in the middle two graphs,
which are simulations of the opponent processes generated by performing
behaviors at different frequencies.

The final graph is the Bode plot, which shows the total hedonic outcome
at each frequency (calculated by taking the integral of the H
compartment over the course of the opponent process simulations). In
this case, performing the behavior at a frequency greater than 0.01
min^(-1) will lead to negative hedonic outcomes.

You can also plot the graphs individually:

``` r
bode_plot(join_plots=FALSE)
```

![](README_files/figure-commonmark/unnamed-chunk-3-1.png)

![](README_files/figure-commonmark/unnamed-chunk-3-2.png)

![](README_files/figure-commonmark/unnamed-chunk-3-3.png)

![](README_files/figure-commonmark/unnamed-chunk-3-4.png)

You can test how varying parameters in the utility function affects the
opponent processes, and in turn affects the hormetic curve. For example,
you can pass a vector to gamma_b to simulate increasing the b-process
curvature in the utility function:

``` r
bode_plot(gamma_b=c(0.7, 0.9, 1.1))
```

![](README_files/figure-commonmark/unnamed-chunk-4-1.png)

Or increase the magnitude of the loss side of the utility function:

``` r
bode_plot(lambda_b=c(2, 2.5, 3))
```

![](README_files/figure-commonmark/unnamed-chunk-5-1.png)

You can also modify the a-process curvature in the utility function:

``` r
bode_plot(gamma_a=c(0.5, 0.7, 0.9), colorscheme=2)
```

![](README_files/figure-commonmark/unnamed-chunk-6-1.png)

Or increase the magnitude of the gain side of the utility function:

``` r
bode_plot(lambda_b=c(1, 1.5, 2), colorscheme=2)
```

![](README_files/figure-commonmark/unnamed-chunk-7-1.png)

You can modify the pharmacodynamic decay constant for the a-process:

``` r
bode_plot(k_apd=c(0.5, 1, 1.5), colorscheme=3)
```

![](README_files/figure-commonmark/unnamed-chunk-8-1.png)

Or modify the pharmacodynamic decay constant for the b-process:

``` r
bode_plot(k_bpd=c(1.5, 2, 2.5), colorscheme=3)
```

![](README_files/figure-commonmark/unnamed-chunk-9-1.png)

Likewise, you can modify the pharmacokinetic decay constant for the
a-process:

``` r
bode_plot(k_apk=c(0.01, 0.02, 0.03), colorscheme=4)
```

![](README_files/figure-commonmark/unnamed-chunk-10-1.png)

Or modify the pharmacokinetic decay constant for the b-process:

``` r
bode_plot(k_bpk=c(0.01, 0.02, 0.03), colorscheme=4)
```

![](README_files/figure-commonmark/unnamed-chunk-11-1.png)

You can extend the frequency range of the Bode plot:

``` r
bode_plot(multiply=300)
```

![](README_files/figure-commonmark/unnamed-chunk-12-1.png)

You can also cut off behaviors at a certain point, so that they stay
within the hormetic limit:

``` r
bode_plot(addl=40)
```

![](README_files/figure-commonmark/unnamed-chunk-13-1.png)

For more usage examples and detailed explanations, please refer to the
package documentation by using the ‘help()’ function in R.

## License

This package is released under the MIT License, and is free to use and
share.
