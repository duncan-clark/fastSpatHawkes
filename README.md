README
================
Duncan Clark
January 16, 2025

# Spatio Temporal Hawkes Process Fitting

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`spatiotemporalHawkes` is a lightweight R package designed to fit
**spatio temporal Hawkes processes**. It is simple and particularly
efficient for small processes with \< 20000 points. Hawkes processes are
a class of self-exciting point process models widely used in seismology,
crime analysis, epidemiology, and other domains where events tend to
cluster in space and time.

## Installation

``` r
devtools::install_github("fastsSpatHawkes")
```

## Example

``` r
# Set some parameters
PARAMS <- list(mu = 5,alpha = 1.5,beta = 5,K = .75)
TIME<- 100
OMEGA <- c(0,10^2,0,10^2)
set.seed(274)

# Generate some data:
data <- sim_hawkes(params = PARAMS,
                   windowT = c(0, TIME),
                   windowS = OMEGA)
# Cleanup into dataframe for fitting
data <- as.data.frame(data)
data <- data[-c(1,5)]

t<-proc.time()[3]
params_fitted <- fit_hawkes(params_init = list(mu = runif(1),
# Fit MLE
                                               alpha = runif(1),
                                               beta = runif(1),
                                               K = runif(1)),
                            realiz = data,
                            windowT = c(0, TIME),
                            windowS = as.owin(OMEGA)
                            )
t<-proc.time()[3]-t
print(params_fitted$par)
```

    ##        mu     alpha      beta         K 
    ## 5.0145437 1.5062312 4.8048690 0.7573824

``` r
print(paste0("Fitting hawkes process to data with ",nrow(data)," points took ",t," seconds."))
```

    ## [1] "Fitting hawkes process to data with 2068 points took 8.833 seconds."
