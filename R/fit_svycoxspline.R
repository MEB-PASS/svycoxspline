#' Fit a survey-weighted Cox proportional hazards model with natural cubic splines
#'
#' @export
fit_svycoxspline <- function(design,
                             exposure_var,
                             time_var,
                             event_var,
                             covariates,
                             trim_probs = c(0.01, 0.99),
                             knot_probs = c(0.05, 0.50, 0.95),
                             ref_prob = 0.10,
                             grid_by = 1,
                             round_to = 1,
                             report_points = NULL) {

  # ------------------------------------------------------------
  # 0. Basic checks
  # ------------------------------------------------------------

  if (!is.null(trim_probs)) {

    if (length(trim_probs) != 2 ||
        trim_probs[1] >= trim_probs[2]) {
      stop("trim_probs must contain two increasing probabilities.")
    }

    if (any(trim_probs <= 0 | trim_probs >= 1)) {
      stop("trim_probs must be between 0 and 1.")
    }

  }

  if (length(knot_probs) < 3) {
    stop(
      "knot_probs must contain at least three probabilities: ",
      "two boundary probabilities and at least one interior knot."
    )
  }

  if (any(knot_probs <= 0 | knot_probs >= 1) ||
      any(diff(knot_probs) <= 0)) {
    stop("knot_probs must be strictly increasing and between 0 and 1.")
  }

  if (ref_prob <= 0 || ref_prob >= 1) {
    stop("ref_prob must be between 0 and 1.")
  }

  if (grid_by <= 0) {
    stop("grid_by must be greater than 0.")
  }

  if (round_to <= 0) {
    stop("round_to must be greater than 0.")
  }

  if (!is.null(report_points)) {

    if (!is.numeric(report_points) ||
        any(!is.finite(report_points))) {
      stop("report_points must be a numeric vector.")
    }

  }

  # ------------------------------------------------------------
  # 1. Trim exposure extremes
  # ------------------------------------------------------------

  if (is.null(trim_probs)) {

    # No trimming
    lo_b <- min(
      design$variables[[exposure_var]],
      na.rm = TRUE
    )

    hi_b <- max(
      design$variables[[exposure_var]],
      na.rm = TRUE
    )

  } else {

    bounds <- survey::svyquantile(
      as.formula(paste0("~", exposure_var)),
      design = design,
      quantiles = trim_probs,
      method = "linear",
      ci = FALSE
    )

    lo_b <- as.numeric(bounds[[1]][1])
    hi_b <- as.numeric(bounds[[1]][2])

    # Subset the existing survey design rather than recreating it.
    # This preserves the survey design structure and weights.
    design <- subset(
      design,
      get(exposure_var) >= lo_b &
        get(exposure_var) <= hi_b
    )
  }

  # Check report points against the modeled exposure range
  if (!is.null(report_points)) {

    if (any(report_points < lo_b | report_points > hi_b)) {
      stop(
        "All report_points must fall within the modeled exposure range: ",
        round(lo_b, 3),
        " to ",
        round(hi_b, 3),
        "."
      )
    }

  }

  # Number of participants and outcome events after trimming
  n_used <- nrow(design$variables)

  n_events <- sum(
    design$variables[[event_var]] == TRUE,
    na.rm = TRUE
  )

  # ------------------------------------------------------------
  # 2. Calculate spline knots and reference value
  # ------------------------------------------------------------

  knots <- as.numeric(
    survey::svyquantile(
      as.formula(paste0("~", exposure_var)),
      design = design,
      quantiles = knot_probs,
      method = "linear",
      ci = FALSE
    )[[1]]
  )

  # Check for duplicate knots
  if (anyDuplicated(knots)) {
    stop(
      "Spline knots are not unique. ",
      "Try different knot probabilities or check the exposure distribution."
    )
  }

  ref_raw <- as.numeric(
    survey::svyquantile(
      as.formula(paste0("~", exposure_var)),
      design = design,
      quantiles = ref_prob,
      method = "linear",
      ci = FALSE
    )[[1]]
  )

  ref_val <- round(ref_raw / round_to) * round_to

  # ------------------------------------------------------------
  # 3. Build spline basis and add to survey design
  # ------------------------------------------------------------

  x <- design$variables[[exposure_var]]

  X_ns <- splines::ns(
    x,
    knots = knots[-c(1, length(knots))],
    Boundary.knots = range(knots)
  )

  nc <- ncol(X_ns)

  colnames(X_ns) <- paste0("s", seq_len(nc))

  update_args <- as.list(as.data.frame(X_ns))

  names(update_args) <- colnames(X_ns)

  design <- do.call(
    update,
    c(list(design), update_args)
  )

  # ------------------------------------------------------------
  # 4. Fit spline Cox model
  # ------------------------------------------------------------

  spline_terms <- paste0(
    "s",
    seq_len(nc),
    collapse = " + "
  )

  covariate_terms <- paste(
    covariates,
    collapse = " + "
  )

  if (length(covariates) == 0) {

    form_spline <- as.formula(
      paste0(
        "survival::Surv(",
        time_var,
        ", ",
        event_var,
        " == TRUE) ~ ",
        spline_terms
      )
    )

  } else {

    form_spline <- as.formula(
      paste0(
        "survival::Surv(",
        time_var,
        ", ",
        event_var,
        " == TRUE) ~ ",
        spline_terms,
        " + ",
        covariate_terms
      )
    )
  }

  fit_spline <- survey::svycoxph(
    form_spline,
    design = design
  )

  # ------------------------------------------------------------
  # 5. Prediction grid
  # ------------------------------------------------------------

  grid <- seq(
    floor(lo_b / round_to) * round_to,
    ceiling(hi_b / round_to) * round_to,
    by = grid_by
  )

  # Spline basis for prediction grid
  X_grid <- as.matrix(
    splines::ns(
      grid,
      knots = knots[-c(1, length(knots))],
      Boundary.knots = range(knots)
    )
  )

  # Spline basis for reference value
  X_ref <- as.numeric(
    splines::ns(
      ref_val,
      knots = knots[-c(1, length(knots))],
      Boundary.knots = range(knots)
    )
  )

  # Contrast each grid value against reference
  L <- t(X_grid) -
    matrix(
      X_ref,
      nrow = nc,
      ncol = nrow(X_grid)
    )

  # Spline coefficients and covariance matrix
  spline_names <- paste0(
    "s",
    seq_len(nc)
  )

  b <- coef(fit_spline)[spline_names]

  V <- vcov(fit_spline)[
    spline_names,
    spline_names,
    drop = FALSE
  ]

  # Log HR and standard error
  lnHR <- as.numeric(
    t(L) %*% b
  )

  SE <- sqrt(
    diag(t(L) %*% V %*% L)
  )

  # Prediction data frame
  pred_df <- data.frame(
    x = grid,
    lnHR = lnHR,
    SE = SE
  )

  pred_df$HR <- exp(pred_df$lnHR)

  pred_df$lo <- exp(
    pred_df$lnHR -
      1.96 * pred_df$SE
  )

  pred_df$hi <- exp(
    pred_df$lnHR +
      1.96 * pred_df$SE
  )

  names(pred_df)[1] <- exposure_var

  # ------------------------------------------------------------
  # 5a. Predictions at user-specified reporting points
  # ------------------------------------------------------------

  if (!is.null(report_points)) {

    X_points <- as.matrix(
      splines::ns(
        report_points,
        knots = knots[-c(1, length(knots))],
        Boundary.knots = range(knots)
      )
    )

    # Contrast each reporting point against reference
    L_points <- t(X_points) -
      matrix(
        X_ref,
        nrow = nc,
        ncol = length(report_points)
      )

    # Log HR and standard error
    lnHR_points <- as.numeric(
      t(L_points) %*% b
    )

    SE_points <- sqrt(
      diag(
        t(L_points) %*% V %*% L_points
      )
    )

    # Reporting point results
    pred_points <- data.frame(
      x = report_points,
      lnHR = lnHR_points,
      SE = SE_points
    )

    pred_points$HR <- exp(
      pred_points$lnHR
    )

    pred_points$lo <- exp(
      pred_points$lnHR -
        1.96 * pred_points$SE
    )

    pred_points$hi <- exp(
      pred_points$lnHR +
        1.96 * pred_points$SE
    )

    names(pred_points)[1] <- exposure_var

  } else {

    pred_points <- NULL

  }

  # ------------------------------------------------------------
  # 6. Fit linear comparison model
  # ------------------------------------------------------------

  if (length(covariates) == 0) {

    form_lin <- as.formula(
      paste0(
        "survival::Surv(",
        time_var,
        ", ",
        event_var,
        " == TRUE) ~ ",
        exposure_var
      )
    )

  } else {

    form_lin <- as.formula(
      paste0(
        "survival::Surv(",
        time_var,
        ", ",
        event_var,
        " == TRUE) ~ ",
        exposure_var,
        " + ",
        covariate_terms
      )
    )
  }

  fit_lin <- survey::svycoxph(
    form_lin,
    design = design
  )

  beta <- coef(fit_lin)[exposure_var]

  vbeta <- vcov(fit_lin)[
    exposure_var,
    exposure_var
  ]

  # Linear prediction data frame
  pred_linear <- data.frame(
    x = grid,
    lnHR = beta * (grid - ref_val)
  )

  pred_linear$SE <- sqrt(vbeta) *
    abs(grid - ref_val)

  pred_linear$HR <- exp(
    pred_linear$lnHR
  )

  pred_linear$lo <- exp(
    pred_linear$lnHR -
      1.96 * pred_linear$SE
  )

  pred_linear$hi <- exp(
    pred_linear$lnHR +
      1.96 * pred_linear$SE
  )

  names(pred_linear)[1] <- exposure_var

  # ------------------------------------------------------------
  # 7. Test for non-linearity
  # ------------------------------------------------------------

  # NOTE:
  # This uses s2:nc as the non-linear terms. This should be
  # validated before treating the resulting p-value as a
  # general-purpose test of non-linearity.
  if (nc > 1) {

    nonlin_terms <- paste0(
      "s",
      2:nc,
      collapse = " + "
    )

    wald_p <- tryCatch({

      p <- survey::regTermTest(
        fit_spline,
        as.formula(
          paste("~", nonlin_terms)
        )
      )$p

      if (is.null(p)) {
        NA_real_
      } else {
        as.numeric(p)
      }

    }, error = function(e) {
      NA_real_
    })

  } else {

    wald_p <- NA_real_

  }

  # ------------------------------------------------------------
  # 8. Format p-value
  # ------------------------------------------------------------

  p_label <- if (is.na(wald_p)) {

    NA_character_

  } else if (wald_p < 0.005) {

    "<0.005"

  } else {

    formatC(
      wald_p,
      format = "f",
      digits = 3
    )
  }

  # ------------------------------------------------------------
  # 9. Interpret p-value
  # ------------------------------------------------------------

  interpretation <- if (is.na(wald_p)) {

    "Non-linearity could not be assessed."

  } else if (wald_p < 0.05) {

    "Evidence of departure from linearity."

  } else {

    "No evidence of departure from linearity."
  }

  # ------------------------------------------------------------
  # 10. Print model summary
  # ------------------------------------------------------------

  message(
    "Spline model fitted.\n",
    "Participants included = ",
    format(n_used, big.mark = ","),
    "\n",
    "Outcome events = ",
    format(n_events, big.mark = ","),
    "\n",
    "P for non-linearity = ",
    p_label,
    "\n",
    interpretation
  )

  # ------------------------------------------------------------
  # 11. Return results
  # ------------------------------------------------------------

  list(
    fit_spline = fit_spline,
    fit_linear = fit_lin,
    pred_spline = pred_df,
    pred_linear = pred_linear,
    pred_points = pred_points,
    knots = knots,
    ref_val = ref_val,
    nonlinearity_p = wald_p,
    nonlinearity_p_label = p_label,
    design_clipped = design
  )
}
