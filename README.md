# Loglinear-BART
loglinear-BART is an R package designed for advanced statistical analysis, particularly suited for dealing with multinomial logistic and count regression models. This package efficiently implements Bayesian Additive Regression Trees (BART) for log-linear models, providing a powerful tool for analyzing categorical and count data. 

## Installation

To install the latest version of `loglinear-BART` directly from GitHub, use the following R command:

```R
# Install the devtools package if it's not already installed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install loglinear-BART from GitHub
devtools::install_github("NingyaW/Loglinear-BART")
