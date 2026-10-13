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

    # Subset existing survey design so survey structure is retained
    design <- subset(
      design,
      get(exposure_var) >= lo_b &
        get(exposure_var) <= hi_b
    )
  }

  # Check reporting points
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

  # Sample size and events
  n_used <- nrow(design$variables)

  n_events <- sum(
    design$variables[[event_var]] == TRUE,
    na.rm = TRUE
  )

  # ------------------------------------------------------------
  # 2. Calculate survey-weighted knots and reference value
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
  # 3. Build natural spline basis
  # ------------------------------------------------------------

  x <- design$variables[[exposure_var]]

  X_ns <- splines::ns(
    x,
    knots = knots[-c(1, length(knots))],
    Boundary.knots = range(knots)
  )

  nc <- ncol(X_ns)

  colnames(X_ns) <- paste0(
    "s",
    seq_len(nc)
  )

  update_args <- as.list(
    as.data.frame(X_ns)
  )

  names(update_args) <- colnames(X_ns)

  design <- do.call(
    stats::update,
    c(list(design), update_args)
  )

  # ------------------------------------------------------------
  # 4. Fit survey-weighted spline Cox model
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
  # 4a. Test for overall spline association
  # ------------------------------------------------------------

  overall_p <- tryCatch({

    p <- survey::regTermTest(
      fit_spline,
      as.formula(
        paste("~", spline_terms)
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

  # ------------------------------------------------------------
  # 5. Prediction grid
  # ------------------------------------------------------------

  grid <- sort(unique(c(
    seq(
      floor(lo_b / round_to) * round_to,
      ceiling(hi_b / round_to) * round_to,
      by = grid_by
    ),
    ref_val
  )))

  X_grid <- as.matrix(
    splines::ns(
      grid,
      knots = knots[-c(1, length(knots))],
      Boundary.knots = range(knots)
    )
  )

  X_ref <- as.numeric(
    splines::ns(
      ref_val,
      knots = knots[-c(1, length(knots))],
      Boundary.knots = range(knots)
    )
  )

  # Contrasts relative to reference
  L <- t(X_grid) -
    matrix(
      X_ref,
      nrow = nc,
      ncol = nrow(X_grid)
    )

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

  lnHR <- as.numeric(
    t(L) %*% b
  )

  SE <- sqrt(
    pmax(
      0,
      diag(t(L) %*% V %*% L)
    )
  )

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

    L_points <- t(X_points) -
      matrix(
        X_ref,
        nrow = nc,
        ncol = length(report_points)
      )

    lnHR_points <- as.numeric(
      t(L_points) %*% b
    )

    SE_points <- sqrt(
      pmax(
        0,
        diag(
          t(L_points) %*% V %*% L_points
        )
      )
    )

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
  # 6. Fit linear comparison model and test linear trend
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

  # Linear coefficient and variance
  beta <- coef(fit_lin)[exposure_var]

  vbeta <- vcov(fit_lin)[
    exposure_var,
    exposure_var
  ]

  # P for linear trend
  trend_p <- summary(fit_lin)$coefficients[
    exposure_var,
    "Pr(>|z|)"
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
  # 7. Formal test for departure from linearity
  # ------------------------------------------------------------
  #
  # Reparameterize the SAME natural spline space so that the
  # original exposure x is explicitly included as the linear
  # component.
  #
  # The remaining basis terms then represent departures from
  # linearity.
  #
  # This avoids assuming that s1 from splines::ns() is itself
  # the linear exposure term.

  if (nc > 1) {

    B <- as.matrix(X_ns)

    # ----------------------------------------------------------
    # 7a. Confirm that x belongs to the natural spline space
    # ----------------------------------------------------------

    basis_fit <- lm.fit(
      x = cbind(Intercept = 1, B),
      y = x
    )

    x_reconstructed <- as.numeric(
      cbind(Intercept = 1, B) %*%
        basis_fit$coefficients
    )

    x_reconstruction_error <- max(
      abs(x - x_reconstructed),
      na.rm = TRUE
    )

    if (!is.finite(x_reconstruction_error) ||
        x_reconstruction_error > 1e-7) {

      stop(
        "The linear exposure could not be reconstructed from ",
        "the natural spline space. Non-linearity test aborted."
      )
    }

    # ----------------------------------------------------------
    # 7b. Replace one spline basis column with x
    # ----------------------------------------------------------

    basis_coef <- basis_fit$coefficients[-1]

    drop_col <- which.max(
      abs(basis_coef)
    )

    keep_cols <- setdiff(
      seq_len(nc),
      drop_col
    )

    X_nonlin <- B[
      ,
      keep_cols,
      drop = FALSE
    ]

    colnames(X_nonlin) <- paste0(
      "nl",
      seq_len(ncol(X_nonlin))
    )

    # ----------------------------------------------------------
    # 7c. Verify that the reparameterized basis spans exactly
    #     the same space as the original natural spline
    # ----------------------------------------------------------

    X_old <- cbind(
      Intercept = 1,
      B
    )

    X_new <- cbind(
      Intercept = 1,
      exposure = x,
      X_nonlin
    )

    reconstruction <- qr.resid(
      qr(X_new),
      X_old
    )

    spline_space_error <- max(
      abs(reconstruction),
      na.rm = TRUE
    )

    if (!is.finite(spline_space_error) ||
        spline_space_error > 1e-7) {

      stop(
        "Reparameterized spline basis does not reproduce ",
        "the original natural spline space. ",
        "Non-linearity test aborted."
      )
    }

    # ----------------------------------------------------------
    # 7d. Add nonlinear terms to survey design
    # ----------------------------------------------------------

    update_nonlin <- as.list(
      as.data.frame(X_nonlin)
    )

    names(update_nonlin) <- colnames(X_nonlin)

    design_nonlin <- do.call(
      stats::update,
      c(list(design), update_nonlin)
    )

    nonlinear_terms <- colnames(
      X_nonlin
    )

    # ----------------------------------------------------------
    # 7e. Fit equivalent reparameterized survey-weighted model
    # ----------------------------------------------------------

    if (length(covariates) == 0) {

      form_nonlin <- as.formula(
        paste0(
          "survival::Surv(",
          time_var,
          ", ",
          event_var,
          " == TRUE) ~ ",
          exposure_var,
          " + ",
          paste(
            nonlinear_terms,
            collapse = " + "
          )
        )
      )

    } else {

      form_nonlin <- as.formula(
        paste0(
          "survival::Surv(",
          time_var,
          ", ",
          event_var,
          " == TRUE) ~ ",
          exposure_var,
          " + ",
          paste(
            nonlinear_terms,
            collapse = " + "
          ),
          " + ",
          covariate_terms
        )
      )
    }

    fit_nonlin <- survey::svycoxph(
      form_nonlin,
      design = design_nonlin
    )

    # ----------------------------------------------------------
    # 7f. Verify fitted model equivalence
    # ----------------------------------------------------------

    lp_original <- as.numeric(
      predict(
        fit_spline,
        type = "lp"
      )
    )

    lp_reparam <- as.numeric(
      predict(
        fit_nonlin,
        type = "lp"
      )
    )

    # Cox linear predictors can differ by an arbitrary constant,
    # so center both before comparing.
    lp_original_centered <-
      lp_original -
      mean(lp_original, na.rm = TRUE)

    lp_reparam_centered <-
      lp_reparam -
      mean(lp_reparam, na.rm = TRUE)

    max_lp_difference <- max(
      abs(
        lp_original_centered -
          lp_reparam_centered
      ),
      na.rm = TRUE
    )

    if (!is.finite(max_lp_difference) ||
        max_lp_difference > 1e-6) {

      stop(
        "Original and reparameterized spline models do not ",
        "produce equivalent linear predictors. ",
        "Non-linearity test aborted."
      )
    }

    # ----------------------------------------------------------
    # 7g. Joint Wald test of nonlinear component(s)
    # ----------------------------------------------------------

    nonlinearity_p <- tryCatch({

      test_result <- survey::regTermTest(
        fit_nonlin,
        as.formula(
          paste(
            "~",
            paste(
              nonlinear_terms,
              collapse = " + "
            )
          )
        )
      )

      if (is.null(test_result$p)) {
        NA_real_
      } else {
        as.numeric(test_result$p)
      }

    }, error = function(e) {

      NA_real_

    })

  } else {

    fit_nonlin <- NULL
    nonlinear_terms <- character(0)
    x_reconstruction_error <- NA_real_
    spline_space_error <- NA_real_
    max_lp_difference <- NA_real_
    nonlinearity_p <- NA_real_
  }

  # ------------------------------------------------------------
  # 8. Format p-values
  # ------------------------------------------------------------

  format_p <- function(p) {

    if (is.na(p)) {

      NA_character_

    } else if (p < 0.005) {

      "<0.005"

    } else {

      formatC(
        p,
        format = "f",
        digits = 3
      )
    }
  }

  overall_p_label <- format_p(overall_p)

  nonlinearity_p_label <- format_p(
    nonlinearity_p
  )

  trend_p_label <- format_p(trend_p)

  # ------------------------------------------------------------
  # 9. Interpret p-value for non-linearity
  # ------------------------------------------------------------

  interpretation <- if (is.na(nonlinearity_p)) {

    "Non-linearity could not be assessed."

  } else if (nonlinearity_p < 0.05) {

    paste0(
      "Wald P-value is statistically significant - ",
      "evidence of departure from linearity."
    )

  } else {

    paste0(
      "Wald P-value is not statistically significant - ",
      "no evidence of departure from linearity."
    )
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
    "P for overall association = ",
    overall_p_label,
    "\n",
    "P for non-linearity = ",
    nonlinearity_p_label,
    "\n",
    "P for linear trend = ",
    trend_p_label,
    "\n",
    interpretation,
    if (nc > 1) {
      paste0(
        "\nMaximum LP difference after reparameterization = ",
        format(
          max_lp_difference,
          scientific = TRUE,
          digits = 3
        )
      )
    } else {
      ""
    }
  )

  # ------------------------------------------------------------
  # 11. Return results
  # ------------------------------------------------------------

  list(
    fit_spline = fit_spline,
    fit_linear = fit_lin,

    # Reparameterized version of same spline model
    fit_nonlinearity = fit_nonlin,

    pred_spline = pred_df,
    pred_linear = pred_linear,
    pred_points = pred_points,

    knots = knots,
    ref_val = ref_val,

    overall_p = overall_p,
    overall_p_label = overall_p_label,

    nonlinearity_p = nonlinearity_p,
    nonlinearity_p_label = nonlinearity_p_label,
    nonlinear_terms = nonlinear_terms,

    trend_p = trend_p,
    trend_p_label = trend_p_label,

    # Diagnostics confirming equivalent parameterization
    x_reconstruction_error = x_reconstruction_error,
    spline_space_error = spline_space_error,
    reparameterization_lp_difference = max_lp_difference,

    design_clipped = design
  )
}
