# svycoxspline

`svycoxspline` is an R package for fitting and visualizing spline associations between continuous exposures and survival outcomes in complex survey data.

The main function, `fit_svycoxspline()`, fits survey-weighted Cox models with natural cubic splines and returns predicted hazard ratios (HRs), 95% confidence intervals, and a test for non-linearity.

## Getting started

Install from GitHub:

``` r
install.packages("remotes")
remotes::install_github("MEB-PASS/svycoxspline")
```

Load the package:

``` r
library(svycoxspline)
```

## Basic use

`fit_svycoxspline()` requires an existing `survey` design object:

``` r
res <- fit_svycoxspline(
  design = survey_design,
  exposure_var = "bmi",
  time_var = "fu_months",
  event_var = "mortality",
  covariates = c("age", "factor(sex)")
)
```

The function can optionally:

-   trim extreme exposure values using survey-weighted quantiles
-   place spline knots using survey-weighted quantiles
-   define the HR reference value using a survey-weighted quantile
-   generate predictions across the exposure range
-   compare the spline association with a linear model
-   test for non-linearity

## Main arguments

| Argument | Description | Default |
|----------------|----------------------------------------|----------------|
| `design` | Existing `survey` design object | \- |
| `exposure_var` | Continuous exposure | \- |
| `time_var` | Follow-up time | \- |
| `event_var` | Event indicator | \- |
| `covariates` | Covariates in the Cox model | `NULL` |
| `trim_probs` | Lower and upper quantiles used to trim the exposure | `c(0.01, 0.99)` |
| `knot_probs` | Quantiles used for spline knots and boundaries | `c(0.05, 0.50, 0.95)` |
| `ref_prob` | Quantile used for the HR = 1 reference | `0.10` |
| `grid_by` | Prediction grid spacing | `1` |
| `round_to` | Rounding unit for the reference and grid limits | `1` |

For example:

``` r
res <- fit_svycoxspline(
  design = survey_design,
  exposure_var = "bmi",
  time_var = "fu_months",
  event_var = "mortality",
  covariates = c("age", "factor(sex)"),
  trim_probs = NULL,
  knot_probs = c(0.05, 0.50, 0.95),
  ref_prob = 0.10
)
```

Set `trim_probs = NULL` to use the full observed exposure range.

## Returned results

The fitted object includes:

-   `fit_spline` - survey-weighted Cox model with spline terms
-   `fit_linear` - corresponding linear Cox model
-   `pred_spline` - predicted HRs and 95% CIs from the spline model
-   `pred_linear` - predicted HRs and 95% CIs from the linear model
-   `knots` - spline knot locations
-   `ref_val` - exposure value corresponding to HR = 1
-   `nonlinearity_p` - p-value for the test of non-linearity
-   `nonlinearity_p_label` - formatted non-linearity p-value
-   `design_clipped` - survey design used for model fitting

For example:

``` r
res$knots
res$ref_val
res$nonlinearity_p_label
head(res$pred_spline)
```

## Survey design and subsetting

Create the survey design and define survey weights before applying analytic or subgroup restrictions.

``` r
survey_design <- svydesign(
  data = dat,
  id = ~psu,
  strata = ~strata,
  weights = ~weight,
  nest = TRUE
)

survey_subset <- subset(
  survey_design,
  age >= 40
)
```

Pass the resulting survey design to `fit_svycoxspline()`.

## Vignette

The package vignette provides a complete example of fitting and visualizing a survey-weighted spline model:

``` r
vignette("svycoxspline")
```

## Reporting issues

Please report comments, suggestions, errors, and bugs through the GitHub Issues page.
