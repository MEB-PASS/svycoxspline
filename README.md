# svycoxspline

`svycoxspline` is an R package for fitting and visualizing associations between continuous exposures and survival outcomes using survey-weighted Cox proportional hazards models with natural cubic splines.

The main function, `fit_svycoxspline()`, fits survey-weighted Cox models with natural cubic splines and returns predicted hazard ratios (HRs), 95% confidence intervals, and p-values for the overall spline association, departure from linearity, and linear trend.

## Getting started

Install from GitHub:

```r
install.packages("remotes")
remotes::install_github("MEB-PASS/svycoxspline")
```

Load the package:

```r
library(svycoxspline)
```

## Survey design and subsetting

Create the survey design and define survey weights before applying analytic or subgroup restrictions.

This example uses generic survey design variables. Complex surveys may use different sampling designs and design variables.

Users should construct the survey design object using the appropriate specification for their data before passing it to `fit_svycoxspline()`.

```r
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

## Basic use

`fit_svycoxspline()` requires an existing `survey` design object:

```r
res <- fit_svycoxspline(
  design = survey_design,
  exposure_var = "bmi",
  time_var = "fu_months",
  event_var = "mortality",
  covariates = c("age", "factor(sex)")
)
```

The function can optionally:

- trim extreme exposure values using survey-weighted quantiles
- place spline knots using survey-weighted quantiles
- define the HR reference value using a survey-weighted quantile
- generate predictions across the exposure range
- return HRs and 95% CIs at user-specified exposure values
- compare the spline association with a linear model
- test the overall spline association
- test for departure from linearity
- test for linear trend

## Main arguments and spline settings

The default settings provide a starting point, but exposure trimming, knot placement, and the reference value should be selected based on the exposure distribution and scientific question. These choices can affect the fitted association and should not be treated as universal defaults (Desquilbet & Mariotti, 2010).

| Argument | Purpose | Default | Example |
|:--------------|:--------------------------|:--------------|:--------------|
| `trim_probs` | Lower and upper survey-weighted quantiles defining the exposure range used for modeling | `c(0.01, 0.99)` | `c(0.01, 0.99)` or `NULL` |
| `knot_probs` | Survey-weighted quantiles used to place spline knots and define the boundary knots | `c(0.05, 0.50, 0.95)` | `c(0.05, 0.50, 0.95)` |
| `ref_prob` | Survey-weighted quantile defining the HR = 1 reference value | `0.10` | `0.10` |
| `report_points` | Exposure values at which HRs and 95% CIs are returned | `NULL` | `c(25, 30)` |
| `grid_by` | Spacing of exposure values used for prediction | `1` | `1` |
| `round_to` | Increment used to round the reference value and prediction limits; `1` rounds to whole units, `0.1` to tenths, and `10` to multiples of 10 | `1` | `1` |

For example, set `trim_probs = NULL` to use the full observed exposure range:

```r
res <- fit_svycoxspline(
  design = survey_design,
  exposure_var = "bmi",
  time_var = "fu_months",
  event_var = "mortality",
  covariates = c("age", "factor(sex)"),
  trim_probs = NULL
)
```

For example, to specify trimming, knots, and exposure values of interest:

```r
res <- fit_svycoxspline(
  design = survey_design,
  exposure_var = "bmi",
  time_var = "fu_months",
  event_var = "mortality",
  covariates = c("age", "factor(sex)"),
  trim_probs = c(0.01, 0.99),
  knot_probs = c(0.05, 0.50, 0.95),
  report_points = c(20, 30),
  ref_prob = 0.10
)
```

## Test for non-linearity

Departure from linearity is assessed using a Wald test. The fitted natural spline is reparameterized so that the exposure is included explicitly as the linear component, with the remaining spline terms representing departures from linearity. The reparameterized model spans the same natural spline space as the model used to estimate the exposure-response curve.

A joint Wald test of the nonlinear terms is used to test the null hypothesis that the exposure-response association is linear. A smaller p-value provides evidence that the association departs from a linear relationship, while a larger p-value indicates insufficient evidence to reject linearity.

The test should be interpreted alongside the estimated spline curve, confidence intervals, exposure distribution, and model assumptions.

## Returned results

The fitted object includes:

- `fit_spline` - survey-weighted Cox model with natural spline terms
- `fit_linear` - corresponding survey-weighted Cox model with a linear exposure term
- `fit_nonlinearity` - equivalent reparameterization of the spline model used for the test of non-linearity
- `pred_spline` - predicted HRs and 95% CIs from the spline model
- `pred_linear` - predicted HRs and 95% CIs from the linear model
- `pred_points` - predicted HRs and 95% CIs at user-specified `report_points`
- `knots` - spline knot locations
- `ref_val` - exposure value corresponding to HR = 1
- `overall_p` - numeric p-value for the overall spline association
- `overall_p_label` - formatted p-value for the overall spline association
- `nonlinearity_p` - numeric p-value for the Wald test of departure from linearity
- `nonlinearity_p_label` - formatted non-linearity p-value
- `nonlinear_terms` - spline terms jointly tested for departure from linearity
- `trend_p` - numeric p-value for the linear trend test
- `trend_p_label` - formatted p-value for the linear trend test
- `design_clipped` - survey design used for model fitting after optional trimming

The three association tests answer different questions:

- `overall_p` tests whether the spline terms are jointly associated with the outcome.
- `nonlinearity_p` tests whether the spline association departs from linearity using a joint Wald test of the nonlinear spline terms.
- `trend_p` tests the continuous exposure term in the corresponding linear Cox model.

The fitted object also contains internal diagnostics that can be used to verify the spline reparameterization:

- `x_reconstruction_error`
- `spline_space_error`
- `reparameterization_lp_difference`

These values should be approximately zero when the reparameterization successfully reproduces the original natural spline model.

For example:

```r
res$knots
res$ref_val

res$overall_p
res$nonlinearity_p
res$trend_p

res$overall_p_label
res$nonlinearity_p_label
res$trend_p_label

head(res$pred_spline)
res$pred_points
```

The numeric p-values can be used directly when creating or exporting results tables, while the corresponding `_label` objects provide formatted values for display.

## Vignette

The package vignette provides a complete example of fitting and visualizing a survey-weighted spline model:

```r
vignette("vignette-svycoxspline")
```

## Reporting issues

Please report comments, suggestions, errors, and bugs through the GitHub Issues page.

## Acknowledgements

Special thanks to Hayden Hayes, Sang Kyu Lee, and Barry Graubard for their contributions to the development and testing of this package.
