# svyspline

`svyspline` is an R package designed to support spline analyses of survival outcomes in complex survey data.

It provides a wrapper for fitting survey-weighted Cox proportional hazards models with natural cubic splines, estimating hazard ratios and 95% confidence intervals across the exposure distribution, and testing for non-linearity.

--------------------------------

**Please note:**
- This package is under active development and functionality may change.
- It is shared to support reproducible research and open science.
- Users are encouraged to review and validate results for their specific application.
- Feedback, bug reports, and suggestions are welcome through the GitHub Issues page.

--------------------------------

## What is `svyspline`?

`svyspline` is designed to simplify spline analyses using complex survey data in R.

The package uses survey-weighted quantiles to define exposure trimming, spline knots, and reference values. It fits survey-weighted Cox proportional hazards models and returns estimates that can be used to visualize the exposure-outcome association.

The package also provides a test for non-linearity and fits a corresponding linear model for comparison.

## Getting started

To install `svyspline` from GitHub:

```r
install.packages("remotes")
remotes::install_github("MEB-PASS/svyspline")
```

The package can then be loaded using:

```r
library(svyspline)
```

## How can `svyspline` be used?

The primary function is:

```r
fit_svyspline()
```

For example:

```r
res <- fit_svyspline(
  design = survey_design,
  exposure_var = "bmi",
  time_var = "fu_months",
  event_var = "mortality",
  covariates = c("age", "factor(sex)"),
  trim_probs = c(0.01, 0.99),
  knot_probs = c(0.05, 0.50, 0.95),
  ref_prob = 0.10,
  grid_by = 1,
  round_to = 1
)
```

The function returns the fitted spline and linear models, predicted hazard ratios and 95% confidence intervals, knot locations, reference value, and test for non-linearity.

## Vignette

A vignette provides an example of fitting and visualizing a survey-weighted spline model.

After installing the package with vignettes, it can be viewed using:

```r
vignette("svyspline")
```

## Reporting issues

Please report comments, suggestions, errors, or bugs through the GitHub Issues page.
