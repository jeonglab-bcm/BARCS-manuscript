# Beta-binomial regression with t-based coefficient tests
#
# This is a dependency-free research implementation.  It fits the mean model
#
#   logit(mu_i) = x_i' beta
#
# using the beta-binomial variance
#
#   Var(K_i) = n_i mu_i (1 - mu_i) {1 + (n_i - 1) rho}.
#
# The coefficient estimates are obtained by feasible IRLS.  rho is estimated
# from the Pearson estimating equation, and coefficient tests use a Student
# t reference distribution with residual degrees of freedom.

.bb_stop <- function(message) {
  stop(message, call. = FALSE)
}

.bb_validate_response <- function(count, total) {
  if (!is.numeric(count) || !is.numeric(total)) {
    .bb_stop("`count` and `total` must be numeric.")
  }
  if (length(count) != length(total) || length(count) < 2L) {
    .bb_stop("`count` and `total` must have the same length (at least two).")
  }
  if (anyNA(count) || anyNA(total) ||
      any(!is.finite(count)) || any(!is.finite(total))) {
    .bb_stop("`count` and `total` cannot contain missing or non-finite values.")
  }
  if (any(total <= 0) || any(count < 0) || any(count > total)) {
    .bb_stop("Each observation must satisfy 0 <= `count` <= `total`, with `total` > 0.")
  }
  integer_like <- function(x) all(abs(x - round(x)) < sqrt(.Machine$double.eps))
  if (!integer_like(count) || !integer_like(total)) {
    .bb_stop("`count` and `total` must contain integer-valued counts.")
  }
  invisible(TRUE)
}

.bb_make_design <- function(formula, data, n) {
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    .bb_stop("`formula` must be a one-sided formula, for example `~ dose + batch`.")
  }
  if (!is.data.frame(data) || nrow(data) != n) {
    .bb_stop("`data` must be a data frame with one row per observation.")
  }
  mf <- model.frame(formula, data = data, na.action = na.fail)
  x <- model.matrix(formula, data = mf)
  qr_x <- qr(x)
  if (qr_x$rank < ncol(x)) {
    .bb_stop("The design matrix is not full rank; remove aliased covariates.")
  }
  if (nrow(x) <= ncol(x)) {
    .bb_stop("The model needs more observations than fitted coefficients.")
  }
  list(x = x, terms = terms(mf), contrasts = attr(x, "contrasts"))
}

.bb_inverse <- function(a) {
  chol2inv(chol(a))
}

# Weighted least squares, the inner step of the IRLS loop.
#
# Both helpers below look for an `*_cpp` function first. Those are the
# RcppArmadillo kernels in the CB2 package: when this file is sourced
# alongside a loaded CB2 they are found and used, and when it is sourced on
# its own -- which is how the benchmark scripts run it, so a reviewer can read
# a benchmark without installing anything -- the R fallback runs instead.
#
# The two paths must stay numerically identical, since results committed to
# this repository may have been produced under either. `tests/run_tests.R`
# covers the R path; the package's own suite covers the C++ one.
.bb_wls_system <- function(x, weight, response) {
  if (exists("bb_wls_system_cpp", mode = "function", inherits = TRUE)) {
    return(bb_wls_system_cpp(x, weight, response))
  }
  list(
    information = crossprod(x, weight * x),
    score_target = crossprod(x, weight * response)
  )
}

.bb_wls_solve <- function(x, weight, response, covariance = FALSE) {
  if (exists("bb_wls_solve_cpp", mode = "function", inherits = TRUE)) {
    result <- bb_wls_solve_cpp(x, weight, response, covariance)
    result$coefficient <- drop(result$coefficient)
    return(result)
  }
  system <- .bb_wls_system(x, weight, response)
  coefficient <- drop(solve(system$information, system$score_target))
  if (!covariance) {
    return(list(coefficient = coefficient))
  }
  list(
    coefficient = coefficient,
    covariance = .bb_inverse(system$information)
  )
}

.bb_estimate_rho <- function(count, total, mu, df_residual,
                             upper = 1 - 1e-8) {
  binomial_variance <- total * mu * (1 - mu)
  pearson <- function(rho) {
    sum((count - total * mu)^2 /
          (binomial_variance * (1 + (total - 1) * rho)))
  }
  q0 <- pearson(0)
  # `pearson_null` is the Pearson chi-square under the binomial model. Unlike
  # `rho` it is not truncated at the lower boundary, so `pearson_null /
  # df_residual` is an unconstrained estimate of the guide's variance
  # inflation and is the quantity `bb_moderate_dispersion()` shrinks.
  if (!is.finite(q0)) {
    return(list(rho = 0, scale = 1, pearson = q0, pearson_null = q0,
                boundary = TRUE))
  }
  if (q0 <= df_residual) {
    return(list(rho = 0, scale = 1, pearson = q0, pearson_null = q0,
                boundary = FALSE))
  }
  qu <- pearson(upper)
  if (qu >= df_residual) {
    return(list(
      rho = upper,
      scale = max(1, qu / df_residual),
      pearson = qu,
      pearson_null = q0,
      boundary = TRUE
    ))
  }
  rho <- uniroot(
    function(r) pearson(r) - df_residual,
    interval = c(0, upper),
    tol = 1e-10
  )$root
  list(rho = rho, scale = 1, pearson = pearson(rho), pearson_null = q0,
       boundary = FALSE)
}

#' Fit beta-binomial regression for one guide
#'
#' @param count Guide counts, one value per sample.
#' @param total Total library counts, one value per sample.
#' @param formula One-sided mean-model formula, such as `~ dose + batch`.
#' @param data Sample-level covariates.
#' @param maxit Maximum feasible-IRLS iterations.
#' @param tolerance Relative convergence tolerance.
#' @param mu_bound Numerical bound applied to fitted proportions.
#' @return An object of class `bbreg`.
#' @export
bbreg <- function(count, total, formula, data, maxit = 100L,
                  tolerance = 1e-8, mu_bound = 1e-8) {
  .bb_validate_response(count, total)
  count <- as.numeric(count)
  total <- as.numeric(total)
  design <- .bb_make_design(formula, data, length(count))
  x <- design$x
  rank <- ncol(x)
  df_residual <- nrow(x) - rank

  initial <- suppressWarnings(glm.fit(
    x = x,
    y = count / total,
    weights = total,
    family = binomial(link = "logit"),
    control = glm.control(maxit = 50L, epsilon = tolerance)
  ))
  beta <- initial$coefficients
  if (any(!is.finite(beta))) {
    pooled <- (sum(count) + 0.5) / (sum(total) + 1)
    beta <- numeric(rank)
    beta[1L] <- qlogis(pooled)
  }

  converged <- FALSE
  rho_fit <- list(rho = 0, scale = 1, pearson = NA_real_,
                  boundary = FALSE)
  for (iteration in seq_len(maxit)) {
    eta <- drop(x %*% beta)
    mu <- pmin(pmax(plogis(eta), mu_bound), 1 - mu_bound)
    rho_fit <- .bb_estimate_rho(count, total, mu, df_residual)
    rho <- rho_fit$rho

    working_response <- eta + (count / total - mu) / (mu * (1 - mu))
    working_weight <- total * mu * (1 - mu) /
      (1 + (total - 1) * rho)
    beta_new <- tryCatch(
      .bb_wls_solve(
        x, working_weight, working_response, covariance = FALSE
      )$coefficient,
      error = function(e) rep(NA_real_, rank)
    )
    if (any(!is.finite(beta_new))) {
      .bb_stop("The IRLS update was singular; inspect sparse counts and the design.")
    }

    change <- max(abs(beta_new - beta) / pmax(1, abs(beta)))
    beta <- beta_new
    if (change < tolerance) {
      converged <- TRUE
      break
    }
  }

  eta <- drop(x %*% beta)
  mu <- pmin(pmax(plogis(eta), mu_bound), 1 - mu_bound)
  rho_fit <- .bb_estimate_rho(count, total, mu, df_residual)
  rho <- rho_fit$rho
  working_weight <- total * mu * (1 - mu) /
    (1 + (total - 1) * rho)
  final_wls <- .bb_wls_solve(
    x, working_weight, eta, covariance = TRUE
  )
  unscaled_covariance <- final_wls$covariance
  covariance <- rho_fit$scale * unscaled_covariance
  standard_error <- sqrt(diag(covariance))
  statistic <- beta / standard_error
  p_value <- 2 * pt(-abs(statistic), df = df_residual)
  coefficient_table <- cbind(
    estimate = beta,
    std_error = standard_error,
    t_value = statistic,
    df = rep(df_residual, rank),
    p_value = p_value
  )
  rownames(coefficient_table) <- colnames(x)

  structure(list(
    coefficients = setNames(beta, colnames(x)),
    coefficient_table = coefficient_table,
    covariance = covariance,
    fitted.values = mu,
    linear.predictors = eta,
    residuals = count / total - mu,
    pearson = rho_fit$pearson,
    pearson_null = rho_fit$pearson_null,
    rho = rho,
    scale = rho_fit$scale,
    dispersion_boundary = rho_fit$boundary,
    df.residual = df_residual,
    rank = rank,
    count = count,
    total = total,
    formula = formula,
    data = data,
    design = x,
    weights = working_weight,
    converged = converged,
    iterations = iteration,
    call = match.call()
  ), class = "bbreg")
}

#' @export
coef.bbreg <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.bbreg <- function(object, ...) {
  object$covariance
}

#' @export
fitted.bbreg <- function(object, ...) {
  object$fitted.values
}

#' @export
residuals.bbreg <- function(object, type = c("response", "pearson"), ...) {
  type <- match.arg(type)
  if (type == "response") {
    return(object$residuals)
  }
  variance <- object$fitted.values * (1 - object$fitted.values) *
    (1 + (object$total - 1) * object$rho) / object$total
  object$residuals / sqrt(variance)
}

#' @export
summary.bbreg <- function(object, ...) {
  ans <- list(
    call = object$call,
    coefficients = object$coefficient_table,
    rho = object$rho,
    pearson = object$pearson,
    df.residual = object$df.residual,
    converged = object$converged,
    iterations = object$iterations
  )
  class(ans) <- "summary.bbreg"
  ans
}

#' @export
print.summary.bbreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, digits = digits, P.values = TRUE,
               has.Pvalue = TRUE)
  cat("\nBeta-binomial intraclass correlation (rho):",
      formatC(x$rho, digits = digits), "\n")
  cat("Pearson statistic / residual df:",
      formatC(x$pearson, digits = digits), "/",
      x$df.residual, "\n")
  cat("Converged:", x$converged, "after", x$iterations, "iterations\n")
  invisible(x)
}

#' Test a linear contrast of beta-binomial regression coefficients
#'
#' @param object A fitted `bbreg` object.
#' @param contrast Numeric contrast vector in coefficient order, or a named
#'   numeric vector whose names identify coefficients.
#' @param null Null value for the contrast.
#' @return A one-row data frame with estimate, standard error, t statistic,
#'   degrees of freedom, and two-sided p-value.
#' @export
bb_contrast <- function(object, contrast, null = 0) {
  if (!inherits(object, "bbreg")) {
    .bb_stop("`object` must be a fitted `bbreg` object.")
  }
  p <- length(object$coefficients)
  if (!is.numeric(contrast) || anyNA(contrast)) {
    .bb_stop("`contrast` must be a numeric vector without missing values.")
  }
  if (!is.null(names(contrast))) {
    if (!all(names(contrast) %in% names(object$coefficients))) {
      .bb_stop("A named contrast contains an unknown coefficient.")
    }
    complete <- setNames(numeric(p), names(object$coefficients))
    complete[names(contrast)] <- contrast
    contrast <- complete
  }
  if (length(contrast) != p) {
    .bb_stop("An unnamed contrast must have one value per coefficient.")
  }
  estimate <- drop(crossprod(contrast, object$coefficients))
  standard_error <- sqrt(drop(t(contrast) %*% object$covariance %*% contrast))
  statistic <- (estimate - null) / standard_error
  data.frame(
    estimate = estimate,
    std_error = standard_error,
    t_value = statistic,
    df = object$df.residual,
    p_value = 2 * pt(-abs(statistic), df = object$df.residual),
    row.names = NULL
  )
}

#' Apply beta-binomial regression guide by guide
#'
#' @param counts Guide-by-sample count matrix.
#' @param data Sample-level covariate data frame.
#' @param formula One-sided regression formula.
#' @param term Coefficient name to report.
#' @param totals Optional library-size vector. Defaults to column sums.
#' @param guide Optional guide identifiers. Defaults to row names.
#' @param gene Optional gene identifier per guide.
#' @param min_total_count Guides below this total count receive missing results.
#' @param ncores Number of forked workers on Unix-like systems. Windows uses
#'   one worker.
#' @param ... Additional arguments passed to `bbreg`.
#' @return A data frame with one row per guide and Benjamini--Hochberg FDR.
#' @export
bb_screen <- function(counts, data, formula, term, totals = NULL,
                      guide = rownames(counts), gene = NULL,
                      min_total_count = 10, ncores = 1L, ...) {
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    .bb_stop("`counts` must be a numeric matrix or data frame.")
  }
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
    .bb_stop("`counts` must contain finite, non-negative values.")
  }
  if (nrow(data) != ncol(counts)) {
    .bb_stop("`data` must have one row per count-matrix column.")
  }
  if (is.null(totals)) {
    totals <- colSums(counts)
  }
  if (length(totals) != ncol(counts)) {
    .bb_stop("`totals` must have one value per count-matrix column.")
  }
  if (any(counts > rep(totals, each = nrow(counts)))) {
    .bb_stop("A guide count cannot exceed its sample's `total`.")
  }
  # Checked here as well as per guide because the per-guide failure is silent:
  # every fit would return an all-NA row with `converged = FALSE`, which reads
  # as a modelling failure rather than a malformed argument. Size-factor
  # normalization is the usual way to arrive with non-integer totals.
  if (any(abs(totals - round(totals)) >= sqrt(.Machine$double.eps))) {
    .bb_stop(paste0(
      "`totals` must be integer-valued library sizes; round them first. ",
      "A beta-binomial denominator counts sequenced reads, so a fractional ",
      "total has no likelihood."
    ))
  }
  if (is.null(guide)) {
    guide <- sprintf("guide_%d", seq_len(nrow(counts)))
  }
  if (length(guide) != nrow(counts) || anyDuplicated(guide)) {
    .bb_stop("`guide` must uniquely identify every row of `counts`.")
  }
  if (!is.null(gene) && length(gene) != nrow(counts)) {
    .bb_stop("`gene` must have one value per guide.")
  }
  if (length(ncores) != 1L || !is.finite(ncores) || ncores < 1) {
    .bb_stop("`ncores` must be one positive integer.")
  }
  ncores <- as.integer(ncores)

  design <- .bb_make_design(formula, data, ncol(counts))
  if (!term %in% colnames(design$x)) {
    .bb_stop(sprintf(
      "`term` must be one model-matrix coefficient: %s",
      paste(colnames(design$x), collapse = ", ")
    ))
  }

  one_guide <- function(i) {
    if (sum(counts[i, ]) < min_total_count) {
      return(c(estimate = NA_real_, std_error = NA_real_,
               t_value = NA_real_, df = NA_real_, p_value = NA_real_,
               rho = NA_real_, pearson_null = NA_real_,
               mean_cpm = mean(counts[i, ] / totals * 1e6),
               converged = 0))
    }
    fit <- tryCatch(
      bbreg(counts[i, ], totals, formula, data, ...),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(c(estimate = NA_real_, std_error = NA_real_,
               t_value = NA_real_, df = NA_real_, p_value = NA_real_,
               rho = NA_real_, pearson_null = NA_real_,
               mean_cpm = mean(counts[i, ] / totals * 1e6),
               converged = 0))
    }
    tab <- fit$coefficient_table[term, ]
    c(
      estimate = tab[["estimate"]],
      std_error = tab[["std_error"]],
      t_value = tab[["t_value"]],
      df = tab[["df"]],
      p_value = tab[["p_value"]],
      rho = fit$rho,
      pearson_null = fit$pearson_null,
      mean_cpm = mean(counts[i, ] / totals * 1e6),
      converged = as.numeric(fit$converged)
    )
  }
  if (ncores > 1L && .Platform$OS.type == "unix") {
    pieces <- parallel::mclapply(
      seq_len(nrow(counts)), one_guide, mc.cores = ncores
    )
    statistics <- do.call(rbind, pieces)
  } else {
    statistics <- t(vapply(
      seq_len(nrow(counts)), one_guide, numeric(9L)
    ))
  }
  result <- data.frame(
    guide = guide,
    statistics,
    row.names = NULL,
    check.names = FALSE
  )
  result$converged <- as.logical(result$converged)
  if (!is.null(gene)) {
    result <- cbind(gene = gene, result)
  }
  result$fdr <- p.adjust(result$p_value, method = "BH")
  result
}

#' Calibrate guide-level t tests with negative-control guides
#'
#' Estimates a one-parameter empirical-null scale from the absolute
#' negative-control t statistics.  The scale matches their empirical
#' `(1 - alpha)` quantile to the corresponding two-sided Student t cutoff.
#' Scales below `min_scale` are truncated, so calibration need not make an
#' already conservative analysis more liberal.
#'
#' @param result A guide-level data frame returned by [bb_screen()].
#' @param control Logical vector identifying negative-control guides.
#' @param alpha Tail probability at which to estimate the null scale.
#' @param min_controls Minimum number of finite negative controls required.
#' @param min_scale Lower bound for the estimated scale.
#' @param method Null-scale estimator. `"tail_quantile"` (default) matches the
#'   empirical `(1 - alpha)` quantile of the absolute control statistics to the
#'   corresponding two-sided Student t cutoff. `"qq_slope"` instead fits the
#'   slope of the control quantile-quantile plot against the t reference over
#'   the 0.50 to 0.95 probability band, which uses many order statistics rather
#'   than one and is markedly less variable.
#' @return `result` with recalibrated standard errors, t statistics, p-values,
#'   and FDR. Original inferential columns are retained with a `raw_` prefix;
#'   the scale and alpha are stored as attributes.
#' @export
bb_calibrate_controls <- function(result, control, alpha = 0.05,
                                  min_controls = 20L, min_scale = 1,
                                  method = c("tail_quantile", "qq_slope")) {
  method <- match.arg(method)
  required <- c("estimate", "std_error", "t_value", "df", "p_value", "fdr")
  if (!is.data.frame(result) || !all(required %in% names(result))) {
    .bb_stop(
      "`result` must be a guide-level result from `bb_screen()`."
    )
  }
  if (!is.logical(control) || length(control) != nrow(result) ||
      anyNA(control)) {
    .bb_stop("`control` must be a non-missing logical vector, one per guide.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 0.5) {
    .bb_stop("`alpha` must be one finite number between 0 and 0.5.")
  }
  if (length(min_controls) != 1L || !is.finite(min_controls) ||
      min_controls < 2) {
    .bb_stop("`min_controls` must be at least two.")
  }
  if (length(min_scale) != 1L || !is.finite(min_scale) ||
      min_scale <= 0) {
    .bb_stop("`min_scale` must be positive.")
  }

  valid <- control & is.finite(result$t_value) & is.finite(result$df)
  if (sum(valid) < as.integer(min_controls)) {
    .bb_stop(sprintf(
      "At least %d finite negative-control statistics are required.",
      as.integer(min_controls)
    ))
  }
  control_df <- unique(result$df[valid])
  if (length(control_df) != 1L || control_df <= 0) {
    .bb_stop("Finite negative-control guides must share one positive `df`.")
  }
  if (method == "tail_quantile") {
    # One order statistic in the far tail. Simple, but with a few hundred
    # controls it is a high-variance estimate of the null scale.
    empirical_cutoff <- as.numeric(stats::quantile(
      abs(result$t_value[valid]),
      probs = 1 - alpha,
      names = FALSE,
      type = 8
    ))
    reference_cutoff <- stats::qt(1 - alpha / 2, df = control_df)
    ratio <- empirical_cutoff / reference_cutoff
  } else {
    # Slope of the control quantile-quantile plot against the t reference,
    # through the origin, over the probability band where a few hundred
    # controls still supply many order statistics. This uses the same null
    # information but averages over the band instead of trusting one point.
    probabilities <- seq(0.50, 0.95, by = 0.01)
    empirical <- as.numeric(stats::quantile(
      abs(result$t_value[valid]),
      probs = probabilities,
      names = FALSE,
      type = 8
    ))
    reference <- stats::qt((1 + probabilities) / 2, df = control_df)
    ratio <- sum(empirical * reference) / sum(reference^2)
  }
  scale <- max(min_scale, ratio)

  result$raw_std_error <- result$std_error
  result$raw_t_value <- result$t_value
  result$raw_p_value <- result$p_value
  result$raw_fdr <- result$fdr
  result$std_error <- result$std_error * scale
  result$t_value <- result$t_value / scale
  result$p_value <- 2 * pt(-abs(result$t_value), df = result$df)
  result$fdr <- p.adjust(result$p_value, method = "BH")
  attr(result, "control_scale") <- scale
  attr(result, "control_alpha") <- alpha
  result
}

#' Invert the trigamma function
#'
#' Newton iteration solving `trigamma(y) = x`, used by the scaled-F moment
#' estimator in [bb_moderate_dispersion()].
#'
#' @param x One positive number.
#' @return The value `y` with `trigamma(y)` equal to `x`, or `Inf`.
#' @keywords internal
# Solve trigamma(y) = x for y.
#
# This is the step that turns an observed spread of log-dispersions into a
# prior degrees of freedom. Smyth's (2004) scaled-F moment estimator gives the
# variance of the log guide dispersions about the fitted trend; under the
# model that variance equals trigamma(prior_df / 2) plus a known term, so
# recovering prior_df means inverting trigamma. There is no closed form, hence
# this. It mirrors `limma::trigammaInverse`.
#
# Newton's method, but on 1/trigamma rather than on trigamma itself. Because
# trigamma is convex and decreasing, iterating on the reciprocal converges
# monotonically from the starting value below, which is what makes a fixed
# 50-iteration cap safe: it is a bound that is never reached in practice, not
# a tolerance.
#
# The three early returns are the boundaries where the iteration is either
# unnecessary or numerically unhelpful:
#
#   x <= 0 or non-finite   No solution -- trigamma is strictly positive. Inf
#                          means "infinite prior df", which the caller reads
#                          as complete shrinkage to the trend. This is the
#                          right answer for a screen whose guide dispersions
#                          show no more spread than sampling alone explains.
#   x > 1e7                trigamma(y) -> 1/y^2 as y -> 0, so y ~ 1/sqrt(x).
#   x < 1e-6               trigamma(y) -> 1/y as y -> Inf, so y ~ 1/x. Here
#                          the Newton step would difference two nearly equal
#                          quantities and lose most of its precision.
.bb_trigamma_inverse <- function(x) {
  if (!is.finite(x) || x <= 0) {
    return(Inf)
  }
  if (x > 1e7) {
    return(1 / sqrt(x))
  }
  if (x < 1e-6) {
    return(1 / x)
  }
  # Starting value chosen to sit below the root, so the monotone convergence
  # above applies.
  y <- 0.5 + 1 / x
  for (iteration in seq_len(50L)) {
    tri <- trigamma(y)
    # Newton step for the reciprocal form; psigamma(y, deriv = 2) is
    # trigamma'(y), the tetragamma function.
    step <- tri * (1 - tri / x) / psigamma(y, deriv = 2L)
    y <- y + step
    if (abs(step / y) < 1e-8) {
      break
    }
  }
  y
}

#' Moderate guide dispersion toward a library-wide trend
#'
#' [bb_screen()] estimates one beta-binomial dispersion per guide from that
#' guide's own residuals.  A pooled screen leaves few residual degrees of
#' freedom per guide, so that estimate is noisy: it inflates some standard
#' errors and deflates others, and the noise carries through to the gene
#' ranking.  This function applies the standard empirical-Bayes remedy.  It
#' shrinks each guide's variance inflation toward a trend fitted across the
#' whole library and returns a moderated `t` statistic on
#' `df_residual + prior_df` degrees of freedom.
#'
#' The shrinkage target and the prior degrees of freedom are estimated from
#' the data by the scaled-F moment method of Smyth (2004), so the function
#' carries no tuning constant that could be fitted to an outcome.
#'
#' Moderation improves power on screen data for a specific reason.  A handful
#' of guides whose own residuals understate their variance drive the
#' negative-control tail, and [bb_calibrate_controls()] can only answer that
#' tail with a blanket penalty applied to every guide in the screen.
#' Correcting those guides individually removes the need for the blanket
#' penalty, so most guides end up tested more sharply even though the noisy
#' ones are tested more strictly.
#'
#' `one_way` and `borrow_df` provide a strictly conservative variant, under
#' which the moderated guide test is never more liberal than the unmoderated
#' one.  It is not the default: on the held-out scenarios it gave up most of
#' the power gain without being needed for calibration.
#'
#' Guide variance inflation is taken as `pearson_null / df_residual`, the
#' Pearson dispersion under the binomial model.  Unlike `rho` it is not
#' truncated at zero, so it behaves like a scaled chi-square statistic and is
#' the right quantity to shrink.  Rescaling the stored standard error by the
#' square root of the moderated-to-fitted inflation ratio is exact when the
#' sample library totals are equal and is a close approximation otherwise.
#'
#' This is guide-level variance moderation.  It is independent of, and can be
#' combined with, the gene-level moderation in [bb_gene_eb_moderate()].
#'
#' @param result Guide-level result returned by [bb_screen()].
#' @param trend Shrink toward an abundance trend in `mean_cpm` (default) or
#'   toward one library-wide constant.
#' @param one_way Only ever raise a guide's variance toward the library
#'   estimate, never lower it.  Off by default.  Lowering the variance of a
#'   genuinely noisy guide is the mechanism that makes shared-dispersion count
#'   models anti-conservative, so this is available as a conservative guard,
#'   but on the held-out CRISPulator scenarios it costs real power and was not
#'   needed to hold the realized false-discovery proportion below nominal.
#' @param borrow_df Claim the prior degrees of freedom in the reference
#'   distribution, as ordinary moderated-`t` inference does.  On by default.
#'   Turning it off uses the moderated variance while keeping the residual
#'   degrees of freedom the design alone supports.
#' @param span Smoother span used for the abundance trend.
#' @param min_guides Minimum number of usable guides required.
#' @return `result` with moderated standard errors, `t` statistics, degrees of
#'   freedom, p-values, and FDR.  The unmoderated columns are retained with an
#'   `unmoderated_` prefix, and the estimated `prior_df` is stored as an
#'   attribute.
#' @export
bb_moderate_dispersion <- function(result, trend = TRUE, one_way = FALSE,
                                   borrow_df = TRUE, span = 0.5,
                                   min_guides = 50L) {
  required <- c("estimate", "std_error", "t_value", "df", "p_value",
                "pearson_null")
  if (!is.data.frame(result) || !all(required %in% names(result))) {
    .bb_stop(
      "`result` must be a guide-level result from `bb_screen()` carrying `pearson_null`."
    )
  }
  if (length(min_guides) != 1L || !is.finite(min_guides) || min_guides < 2) {
    .bb_stop("`min_guides` must be at least two.")
  }

  usable <- is.finite(result$pearson_null) & result$pearson_null >= 0 &
    is.finite(result$std_error) & result$std_error > 0 &
    is.finite(result$estimate) & is.finite(result$df) & result$df > 0
  if ("converged" %in% names(result)) {
    usable <- usable & !is.na(result$converged) & result$converged
  }
  if (sum(usable) < as.integer(min_guides)) {
    .bb_stop(sprintf(
      "At least %d usable guide fits are required to moderate dispersion.",
      as.integer(min_guides)
    ))
  }
  df_residual <- unique(result$df[usable])
  if (length(df_residual) != 1L) {
    .bb_stop("Usable guides must share one residual degrees-of-freedom value.")
  }

  inflation <- pmax(result$pearson_null[usable] / df_residual, 1e-8)
  log_inflation <- log(inflation)

  if (isTRUE(trend) && "mean_cpm" %in% names(result)) {
    abundance <- log(pmax(result$mean_cpm[usable], .Machine$double.eps))
    ordering <- order(abundance)
    smooth <- stats::lowess(
      abundance[ordering], log_inflation[ordering], f = span
    )
    trend_log <- stats::approx(
      smooth$x, smooth$y, xout = abundance, rule = 2
    )$y
  } else {
    trend_log <- rep(mean(log_inflation), length(log_inflation))
  }

  # Scaled-F moment estimator (Smyth 2004) on the log dispersion scale.
  residual_log <- log_inflation - trend_log
  centred <- residual_log - digamma(df_residual / 2) + log(df_residual / 2)
  n_guides <- length(centred)
  excess_variance <- n_guides / (n_guides - 1) *
    mean((centred - mean(centred))^2) - trigamma(df_residual / 2)
  if (excess_variance > 0) {
    prior_df <- 2 * .bb_trigamma_inverse(excess_variance)
    location <- mean(centred) + digamma(prior_df / 2) - log(prior_df / 2)
  } else {
    prior_df <- Inf
    location <- mean(centred)
  }
  prior_inflation <- exp(trend_log + location)

  prior_df_used <- min(prior_df, 1e6)
  moderated <- (df_residual * inflation + prior_df_used * prior_inflation) /
    (df_residual + prior_df_used)

  # `bb_screen()` already scaled the standard error by the fitted inflation,
  # which is truncated at one because `rho` cannot be negative.
  fitted_inflation <- pmax(1, inflation)
  if (isTRUE(one_way)) {
    # Raise a guide's variance toward the library when its own residuals look
    # implausibly quiet, but never lower it below what its own data support.
    # Shrinking a genuinely noisy guide down is what makes shared-dispersion
    # count models anti-conservative on screen data, and this screens that
    # failure mode out while keeping the correction that matters.
    moderated <- pmax(moderated, fitted_inflation)
  }
  df_total <- if (isTRUE(borrow_df)) {
    df_residual + prior_df_used
  } else {
    # Use the better variance estimate without also claiming the degrees of
    # freedom the prior would buy. The extra df is the part of moderation that
    # depends on the prior being exactly right.
    df_residual
  }
  rescale <- sqrt(moderated / fitted_inflation)

  moderated_result <- result
  moderated_result$unmoderated_std_error <- result$std_error
  moderated_result$unmoderated_t_value <- result$t_value
  moderated_result$unmoderated_df <- result$df
  moderated_result$unmoderated_p_value <- result$p_value
  moderated_result$variance_inflation <- NA_real_
  moderated_result$variance_inflation[usable] <- inflation
  moderated_result$moderated_inflation <- NA_real_
  moderated_result$moderated_inflation[usable] <- moderated

  moderated_result$std_error[!usable] <- NA_real_
  moderated_result$t_value[!usable] <- NA_real_
  moderated_result$df[!usable] <- NA_real_
  moderated_result$p_value[!usable] <- NA_real_
  moderated_result$std_error[usable] <- result$std_error[usable] * rescale
  moderated_result$t_value[usable] <- result$estimate[usable] /
    moderated_result$std_error[usable]
  moderated_result$df[usable] <- df_total
  moderated_result$p_value[usable] <- 2 * stats::pt(
    -abs(moderated_result$t_value[usable]), df = df_total
  )
  moderated_result$fdr <- stats::p.adjust(
    moderated_result$p_value, method = "BH"
  )
  attr(moderated_result, "prior_df") <- prior_df
  attr(moderated_result, "df_total") <- df_total
  moderated_result
}

#' Test a shared guide effect against an empirical gene-level null
#'
#' This function is intended for exploratory screens in which several
#' independently designed guides target each gene but biological replication
#' is too limited for reliable guide-level reference distributions. It does
#' not treat guides as biological replicates. Instead, it estimates one shared
#' gene effect by inverse-variance weighting of guide coefficients, forms its
#' model-based Wald statistic, then calibrates the gene-statistic distribution
#' with a robust empirical null. It does not combine guide p-values by Fisher's
#' or Stouffer's method.
#'
#' For gene \(g\), let \(w_{gj}=\operatorname{SE}(\widehat\beta_{gj})^{-2}\).
#' The shared effect and raw statistic are
#' \deqn{\widehat\beta_g =
#' \frac{\sum_j w_{gj}\widehat\beta_{gj}}{\sum_j w_{gj}},\qquad
#' T_g = \widehat\beta_g\sqrt{\sum_jw_{gj}}.}
#' Its null center is the median score among control genes when enough are
#' supplied, and otherwise the median across all genes. The null scale is the
#' largest of `min_scale`, the all-gene MAD, and the control-gene tail scale.
#' The all-gene MAD assumes that fewer than half of genes are active.
#'
#' @param result A guide-level result from [bb_screen()]. If control
#'   calibration has already been applied, the retained `raw_std_error`
#'   column is used automatically.
#' @param control Optional non-missing logical vector identifying control
#'   guides. A control gene must contain only control guides.
#' @param min_guides Minimum number of finite guide scores required per gene.
#' @param alpha Tail probability used to estimate the control-gene scale.
#' @param min_control_genes Minimum number of valid control genes needed to
#'   use their median and tail scale.
#' @param min_scale Lower bound for the empirical-null scale.
#' @return A data frame with one row per testable gene. `statistic` is the
#'   empirical-null standardized gene score and `p_value` uses a standard
#'   normal reference. Null parameters are also stored as attributes.
#' @export
bb_gene_consistency <- function(result, control = NULL, min_guides = 3L,
                                alpha = 0.05, min_control_genes = 10L,
                                min_scale = 1) {
  required <- c("gene", "estimate", "std_error")
  if (!is.data.frame(result) || !all(required %in% names(result))) {
    .bb_stop(
      "`result` must contain guide-level `gene`, `estimate`, and `std_error` columns."
    )
  }
  if (anyNA(result$gene)) {
    .bb_stop("`result$gene` cannot contain missing values.")
  }
  if (length(min_guides) != 1L || !is.finite(min_guides) ||
      min_guides < 2) {
    .bb_stop("`min_guides` must be an integer of at least two.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 0.5) {
    .bb_stop("`alpha` must be one finite number between 0 and 0.5.")
  }
  if (length(min_control_genes) != 1L ||
      !is.finite(min_control_genes) || min_control_genes < 2) {
    .bb_stop("`min_control_genes` must be an integer of at least two.")
  }
  if (length(min_scale) != 1L || !is.finite(min_scale) ||
      min_scale <= 0) {
    .bb_stop("`min_scale` must be positive.")
  }
  min_guides <- as.integer(min_guides)
  min_control_genes <- as.integer(min_control_genes)

  if (is.null(control)) {
    control <- rep(FALSE, nrow(result))
  } else if (!is.logical(control) || length(control) != nrow(result) ||
             anyNA(control)) {
    .bb_stop("`control` must be a non-missing logical vector, one per guide.")
  }
  control_by_gene <- split(control, result$gene)
  mixed_control <- vapply(
    control_by_gene,
    function(value) any(value) && !all(value),
    logical(1L)
  )
  if (any(mixed_control)) {
    .bb_stop("A gene cannot mix control and non-control guides.")
  }

  standard_error <- if ("raw_std_error" %in% names(result)) {
    result$raw_std_error
  } else {
    result$std_error
  }
  valid <- is.finite(result$estimate) &
    is.finite(standard_error) &
    standard_error > 0
  groups <- split(seq_len(nrow(result)), result$gene)
  pieces <- lapply(names(groups), function(gene_name) {
    index <- groups[[gene_name]]
    index <- index[valid[index]]
    if (length(index) < min_guides) {
      return(NULL)
    }
    guide_weight <- 1 / standard_error[index]^2
    gene_estimate <- sum(
      guide_weight * result$estimate[index]
    ) / sum(guide_weight)
    gene_standard_error <- sqrt(1 / sum(guide_weight))
    nonzero <- result$estimate[index] != 0
    agreement <- if (gene_estimate == 0 || !any(nonzero)) {
      NA_real_
    } else {
      mean(
        sign(result$estimate[index][nonzero]) == sign(gene_estimate)
      )
    }
    data.frame(
      gene = gene_name,
      n_guides = length(index),
      estimate = gene_estimate,
      std_error = gene_standard_error,
      raw_statistic = gene_estimate / gene_standard_error,
      guide_direction_agreement = agreement,
      converged_fraction = if ("converged" %in% names(result)) {
        mean(result$converged[index], na.rm = TRUE)
      } else {
        NA_real_
      },
      control_gene = all(control[index]),
      row.names = NULL
    )
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (length(pieces) < 2L) {
    .bb_stop("At least two genes must have enough finite guide scores.")
  }
  gene_result <- do.call(rbind, pieces)
  rownames(gene_result) <- NULL

  global_center <- stats::median(gene_result$raw_statistic)
  global_scale <- stats::mad(
    gene_result$raw_statistic,
    center = global_center,
    constant = 1 / stats::qnorm(0.75)
  )
  if (!is.finite(global_scale) || global_scale <= 0) {
    global_scale <- 1
  }

  control_statistic <- gene_result$raw_statistic[
    gene_result$control_gene
  ]
  enough_controls <- length(control_statistic) >= min_control_genes
  null_center <- if (enough_controls) {
    stats::median(control_statistic)
  } else {
    global_center
  }
  control_scale <- if (enough_controls) {
    as.numeric(stats::quantile(
      abs(control_statistic - null_center),
      probs = 1 - alpha,
      names = FALSE,
      type = 8
    )) / stats::qnorm(1 - alpha / 2)
  } else {
    NA_real_
  }
  scale_candidates <- c(min_scale, global_scale, control_scale)
  null_scale <- max(scale_candidates[is.finite(scale_candidates)])

  gene_result$statistic <-
    (gene_result$raw_statistic - null_center) / null_scale
  gene_result$p_value <- 2 * stats::pnorm(-abs(gene_result$statistic))
  gene_result$fdr <- stats::p.adjust(gene_result$p_value, method = "BH")
  attr(gene_result, "null_center") <- null_center
  attr(gene_result, "null_scale") <- null_scale
  attr(gene_result, "global_scale") <- global_scale
  attr(gene_result, "control_scale") <- control_scale
  attr(gene_result, "control_genes") <- length(control_statistic)
  attr(gene_result, "null_assumption") <-
    paste(
      "Shared-effect guide-consistency empirical null;",
      "not biological-replicate inference."
    )
  gene_result
}

.bb_gene_pooling_inputs <- function(result, control, min_guides) {
  required <- c("gene", "estimate", "std_error")
  if (!is.data.frame(result) || !all(required %in% names(result))) {
    .bb_stop(
      "`result` must contain guide-level `gene`, `estimate`, and `std_error` columns."
    )
  }
  if (anyNA(result$gene) || any(!nzchar(as.character(result$gene)))) {
    .bb_stop("`result$gene` must contain non-missing gene identifiers.")
  }
  if (length(min_guides) != 1L || !is.finite(min_guides) ||
      min_guides < 2) {
    .bb_stop("`min_guides` must be an integer of at least two.")
  }
  if (is.null(control)) {
    control <- rep(FALSE, nrow(result))
  } else if (!is.logical(control) || length(control) != nrow(result) ||
             anyNA(control)) {
    .bb_stop("`control` must be a non-missing logical vector, one per guide.")
  }
  control_by_gene <- split(control, result$gene)
  mixed_control <- vapply(
    control_by_gene,
    function(value) any(value) && !all(value),
    logical(1L)
  )
  if (any(mixed_control)) {
    .bb_stop("A gene cannot mix control and non-control guides.")
  }
  standard_error <- if ("raw_std_error" %in% names(result)) {
    result$raw_std_error
  } else {
    result$std_error
  }
  valid <- is.finite(result$estimate) &
    is.finite(standard_error) &
    standard_error > 0
  if ("converged" %in% names(result)) {
    valid <- valid & !is.na(result$converged) & result$converged
  }
  list(
    result = result,
    control = control,
    standard_error = standard_error,
    valid = valid,
    min_guides = as.integer(min_guides)
  )
}

#' Reproduce the original BARCS guide-to-gene statistic
#'
#' Converts each two-sided guide p-value to a signed standard-normal score and
#' sums those scores within genes. This function exists to keep the historical
#' BARCS benchmark calculation explicit while newer effect-pooling methods are
#' evaluated beside it.
#'
#' @param result Guide-level result returned by [bb_screen()], optionally
#'   calibrated by [bb_calibrate_controls()].
#' @param min_guides Minimum number of finite guide results required per gene.
#' @return One row per testable gene with the historical signed-z statistic.
#' @export
bb_gene_original <- function(result, min_guides = 1L) {
  required <- c("gene", "estimate", "p_value")
  if (!is.data.frame(result) || !all(required %in% names(result))) {
    .bb_stop(
      "`result` must contain guide-level `gene`, `estimate`, and `p_value` columns."
    )
  }
  if (anyNA(result$gene) || any(!nzchar(as.character(result$gene)))) {
    .bb_stop("`result$gene` must contain non-missing gene identifiers.")
  }
  if (length(min_guides) != 1L || !is.finite(min_guides) ||
      min_guides < 1) {
    .bb_stop("`min_guides` must be one positive integer.")
  }
  min_guides <- as.integer(min_guides)
  valid <- is.finite(result$estimate) &
    is.finite(result$p_value) &
    result$p_value >= 0 &
    result$p_value <= 1
  if ("converged" %in% names(result)) {
    valid <- valid & !is.na(result$converged) & result$converged
  }
  groups <- split(seq_len(nrow(result)), result$gene)
  pieces <- lapply(names(groups), function(gene_name) {
    all_index <- groups[[gene_name]]
    index <- all_index[valid[all_index]]
    if (length(index) < min_guides) {
      return(NULL)
    }
    signed_z <- sign(result$estimate[index]) * stats::qnorm(
      pmax(result$p_value[index] / 2, .Machine$double.xmin),
      lower.tail = FALSE
    )
    combined_z <- sum(signed_z) / sqrt(length(signed_z))
    gene_estimate <- stats::median(result$estimate[index])
    data.frame(
      gene = gene_name,
      n_guides = length(index),
      estimate = gene_estimate,
      statistic = combined_z,
      guide_direction_agreement = if (gene_estimate == 0) {
        NA_real_
      } else {
        mean(sign(result$estimate[index]) == sign(gene_estimate))
      },
      effect_statistic_sign_agreement =
        gene_estimate == 0 || sign(gene_estimate) == sign(combined_z),
      p_value = 2 * stats::pnorm(-abs(combined_z)),
      converged_fraction = if ("converged" %in% names(result)) {
        mean(result$converged[all_index], na.rm = TRUE)
      } else {
        NA_real_
      },
      method = "original",
      row.names = NULL
    )
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) {
    .bb_stop("No gene has enough finite guide results.")
  }
  gene_result <- do.call(rbind, pieces)
  rownames(gene_result) <- NULL
  gene_result$fdr <- stats::p.adjust(gene_result$p_value, method = "BH")
  gene_result
}

#' Test whether normally distributed guide coefficients have nonzero mean
#'
#' Treats the fitted guide coefficients within each gene as exchangeable
#' observations
#' \(\widehat\beta_{gj} \sim N(\mu_g, \sigma_g^2)\). The gene effect is the
#' arithmetic mean, the guide-effect standard deviation is estimated by the
#' ordinary sample standard deviation, and the standard error of the mean is
#' \(s_g/\sqrt{m_g}\). Because \(\sigma_g\) is estimated rather than known,
#' the default p-value uses a Student t reference with \(m_g-1\) degrees of
#' freedom. A plug-in standard-normal p-value is returned separately and can
#' be selected explicitly for sensitivity analysis.
#'
#' This is guide-consistency inference, not biological-replicate inference.
#' It deliberately gives every valid guide equal weight and does not use the
#' guide-specific regression standard errors.
#'
#' @param result Guide-level result returned by [bb_screen()].
#' @param min_guides Minimum number of finite guide coefficients per gene.
#' @param reference Reference distribution for the reported `p_value`:
#'   `"student_t"` (default) or the plug-in `"normal"` approximation.
#' @return One row per testable gene with the estimated normal mean and
#'   standard deviation, standard error of the mean, t/z statistics, both
#'   reference p-values, selected p-value, and FDR.
#' @export
bb_gene_normal <- function(result, min_guides = 3L,
                           reference = c("student_t", "normal")) {
  required <- c("gene", "estimate")
  if (!is.data.frame(result) || !all(required %in% names(result))) {
    .bb_stop("`result` must contain guide-level `gene` and `estimate` columns.")
  }
  if (anyNA(result$gene) || any(!nzchar(as.character(result$gene)))) {
    .bb_stop("`result$gene` must contain non-missing gene identifiers.")
  }
  if (length(min_guides) != 1L || !is.finite(min_guides) ||
      min_guides < 2) {
    .bb_stop("`min_guides` must be an integer of at least two.")
  }
  reference <- match.arg(reference)
  min_guides <- as.integer(min_guides)
  valid <- is.finite(result$estimate)
  if ("converged" %in% names(result)) {
    valid <- valid & !is.na(result$converged) & result$converged
  }

  groups <- split(seq_len(nrow(result)), result$gene)
  pieces <- lapply(names(groups), function(gene_name) {
    all_index <- groups[[gene_name]]
    index <- all_index[valid[all_index]]
    n_guides <- length(index)
    if (n_guides < min_guides) {
      return(NULL)
    }
    guide_beta <- result$estimate[index]
    mu <- mean(guide_beta)
    sigma <- stats::sd(guide_beta)
    standard_error <- sigma / sqrt(n_guides)
    statistic <- if (standard_error > 0) {
      mu / standard_error
    } else if (mu == 0) {
      0
    } else {
      sign(mu) * Inf
    }
    degrees_freedom <- n_guides - 1L
    student_p_value <- 2 * stats::pt(
      -abs(statistic), df = degrees_freedom
    )
    normal_p_value <- 2 * stats::pnorm(-abs(statistic))
    data.frame(
      gene = gene_name,
      n_guides = n_guides,
      estimate = mu,
      sigma = sigma,
      std_error = standard_error,
      statistic = statistic,
      df = degrees_freedom,
      student_p_value = student_p_value,
      normal_p_value = normal_p_value,
      p_value = if (reference == "student_t") {
        student_p_value
      } else {
        normal_p_value
      },
      guide_direction_agreement = if (mu == 0) {
        NA_real_
      } else {
        mean(sign(guide_beta) == sign(mu))
      },
      converged_fraction = if ("converged" %in% names(result)) {
        mean(result$converged[all_index], na.rm = TRUE)
      } else {
        NA_real_
      },
      method = paste0("normal_beta_", reference),
      row.names = NULL
    )
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) {
    .bb_stop("No gene has enough finite guide coefficients.")
  }
  gene_result <- do.call(rbind, pieces)
  rownames(gene_result) <- NULL
  gene_result$fdr <- stats::p.adjust(gene_result$p_value, method = "BH")
  attr(gene_result, "reference") <- reference
  attr(gene_result, "null_assumption") <- paste(
    "Exchangeable normally distributed guide coefficients;",
    "not biological-replicate inference."
  )
  gene_result
}

.bb_gene_pooling_components <- function(inputs) {
  result <- inputs$result
  standard_error <- inputs$standard_error
  valid <- inputs$valid
  control <- inputs$control
  groups <- split(seq_len(nrow(result)), result$gene)
  pieces <- lapply(names(groups), function(gene_name) {
    all_index <- groups[[gene_name]]
    index <- all_index[valid[all_index]]
    if (length(index) < inputs$min_guides) {
      return(NULL)
    }
    estimate <- result$estimate[index]
    variance <- standard_error[index]^2
    fixed_weight <- 1 / variance
    fixed_estimate <- sum(fixed_weight * estimate) / sum(fixed_weight)
    q <- sum(fixed_weight * (estimate - fixed_estimate)^2)
    heterogeneity_df <- length(index) - 1
    c_value <- sum(fixed_weight) -
      sum(fixed_weight^2) / sum(fixed_weight)
    tau2 <- if (is.finite(c_value) && c_value > 0) {
      max(0, (q - heterogeneity_df) / c_value)
    } else {
      0
    }
    direction_agreement <- if (fixed_estimate == 0) {
      NA_real_
    } else {
      mean(sign(estimate) == sign(fixed_estimate))
    }
    data.frame(
      gene = gene_name,
      n_guides = length(index),
      fixed_estimate = fixed_estimate,
      fixed_information = sum(fixed_weight),
      heterogeneity_q = q,
      heterogeneity_df = heterogeneity_df,
      heterogeneity_c = c_value,
      raw_tau2 = tau2,
      i_squared = if (q > 0) {
        max(0, (q - heterogeneity_df) / q)
      } else {
        0
      },
      guide_direction_agreement = direction_agreement,
      converged_fraction = if ("converged" %in% names(result)) {
        mean(result$converged[all_index], na.rm = TRUE)
      } else {
        NA_real_
      },
      control_gene = all(control[all_index]),
      row.names = NULL
    )
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (length(pieces) < 2L) {
    .bb_stop("At least two genes must have enough finite guide scores.")
  }
  components <- do.call(rbind, pieces)
  rownames(components) <- NULL
  components
}

.bb_calibrate_gene_statistics <- function(result, alpha,
                                          min_control_genes,
                                          min_scale) {
  if (length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 0.5) {
    .bb_stop("`alpha` must be one finite number between 0 and 0.5.")
  }
  if (length(min_control_genes) != 1L ||
      !is.finite(min_control_genes) || min_control_genes < 2) {
    .bb_stop("`min_control_genes` must be an integer of at least two.")
  }
  if (length(min_scale) != 1L || !is.finite(min_scale) ||
      min_scale <= 0) {
    .bb_stop("`min_scale` must be positive.")
  }
  global_center <- stats::median(result$raw_statistic)
  global_scale <- stats::mad(
    result$raw_statistic,
    center = global_center,
    constant = 1 / stats::qnorm(0.75)
  )
  if (!is.finite(global_scale) || global_scale <= 0) {
    global_scale <- 1
  }
  control_statistic <- result$raw_statistic[result$control_gene]
  enough_controls <-
    length(control_statistic) >= as.integer(min_control_genes)
  null_center <- if (enough_controls) {
    stats::median(control_statistic)
  } else {
    global_center
  }
  control_scale <- if (enough_controls) {
    as.numeric(stats::quantile(
      abs(control_statistic - null_center),
      probs = 1 - alpha,
      names = FALSE,
      type = 8
    )) / stats::qnorm(1 - alpha / 2)
  } else {
    NA_real_
  }
  scale_candidates <- c(min_scale, global_scale, control_scale)
  null_scale <- max(scale_candidates[is.finite(scale_candidates)])
  result$statistic <-
    (result$raw_statistic - null_center) / null_scale
  result$p_value <- 2 * stats::pnorm(-abs(result$statistic))
  result$fdr <- stats::p.adjust(result$p_value, method = "BH")
  attr(result, "null_center") <- null_center
  attr(result, "null_scale") <- null_scale
  attr(result, "global_scale") <- global_scale
  attr(result, "control_scale") <- control_scale
  attr(result, "control_genes") <- length(control_statistic)
  result
}

.bb_finish_gene_pooling <- function(inputs, components, tau2,
                                    method, alpha,
                                    min_control_genes, min_scale) {
  result <- inputs$result
  standard_error <- inputs$standard_error
  valid <- inputs$valid
  groups <- split(seq_len(nrow(result)), result$gene)
  rows <- lapply(seq_len(nrow(components)), function(row_index) {
    gene_name <- components$gene[row_index]
    all_index <- groups[[gene_name]]
    index <- all_index[valid[all_index]]
    estimate <- result$estimate[index]
    variance <- standard_error[index]^2
    gene_tau2 <- tau2[row_index]
    weight <- 1 / (variance + gene_tau2)
    gene_estimate <- sum(weight * estimate) / sum(weight)
    gene_standard_error <- sqrt(1 / sum(weight))
    posterior_weight <- gene_tau2 / (gene_tau2 + variance)
    shrunken_guide_effect <-
      gene_estimate + posterior_weight * (estimate - gene_estimate)
    leave_one_out <- vapply(seq_along(index), function(drop_index) {
      keep <- seq_along(index) != drop_index
      sum(weight[keep] * estimate[keep]) / sum(weight[keep])
    }, numeric(1L))
    data.frame(
      gene = gene_name,
      n_guides = length(index),
      estimate = gene_estimate,
      std_error = gene_standard_error,
      raw_statistic = gene_estimate / gene_standard_error,
      tau2 = gene_tau2,
      raw_tau2 = components$raw_tau2[row_index],
      heterogeneity_q = components$heterogeneity_q[row_index],
      heterogeneity_df = components$heterogeneity_df[row_index],
      i_squared = components$i_squared[row_index],
      guide_direction_agreement =
        components$guide_direction_agreement[row_index],
      max_weight_fraction = max(weight / sum(weight)),
      max_shrinkage = max(abs(shrunken_guide_effect - estimate)),
      leave_one_out_max_change =
        max(abs(leave_one_out - gene_estimate)),
      leave_one_out_sign_stable = if (gene_estimate == 0) {
        NA
      } else {
        all(sign(leave_one_out) == sign(gene_estimate))
      },
      converged_fraction = components$converged_fraction[row_index],
      control_gene = components$control_gene[row_index],
      method = method,
      row.names = NULL
    )
  })
  pooled <- do.call(rbind, rows)
  rownames(pooled) <- NULL
  .bb_calibrate_gene_statistics(
    pooled,
    alpha = alpha,
    min_control_genes = min_control_genes,
    min_scale = min_scale
  )
}

#' Combine guide effects by random-effects partial pooling
#'
#' Preserves every guide-level BARCS coefficient and standard error. For each
#' gene, a DerSimonian--Laird guide-heterogeneity variance is estimated and
#' guide effects are combined with weights
#' \(1 / \{\operatorname{SE}_{gj}^2 + \tau_g^2\}\). Guide count is not used as
#' biological residual degrees of freedom. The resulting working statistic
#' is calibrated with a robust gene-level empirical null.
#'
#' @param result Guide-level result returned by [bb_screen()].
#' @param control Optional logical vector identifying negative-control guides.
#' @param min_guides Minimum number of finite guides required per gene.
#' @param alpha Tail probability used for empirical-null calibration.
#' @param min_control_genes Minimum valid control genes required to use the
#'   control-only null.
#' @param min_scale Lower bound for empirical-null scale.
#' @return One row per testable gene with the partially pooled effect,
#'   uncertainty, heterogeneity, influence diagnostics, p-value, and FDR.
#' @export
bb_gene_partial_pool <- function(result, control = NULL, min_guides = 2L,
                                 alpha = 0.05,
                                 min_control_genes = 10L,
                                 min_scale = 1) {
  inputs <- .bb_gene_pooling_inputs(result, control, min_guides)
  components <- .bb_gene_pooling_components(inputs)
  pooled <- .bb_finish_gene_pooling(
    inputs,
    components,
    tau2 = components$raw_tau2,
    method = "partial_pooling",
    alpha = alpha,
    min_control_genes = min_control_genes,
    min_scale = min_scale
  )
  attr(pooled, "heterogeneity_estimator") <- "DerSimonian-Laird"
  attr(pooled, "null_assumption") <-
    "Guide random-effects model with a robust gene-level empirical null."
  pooled
}

#' Combine guide effects with empirical-Bayes heterogeneity moderation
#'
#' Starts from the same random-effects guide model as
#' [bb_gene_partial_pool()]. A screen-wide heterogeneity variance is estimated
#' by pooling excess Cochran Q across genes. Each noisy gene-specific
#' heterogeneity estimate is then shrunk toward that screen-wide value:
#' \deqn{\widetilde\tau_g^2 =
#' \frac{d_g\widehat\tau_g^2 + d_0\tau_0^2}{d_g+d_0}.}
#' The prior scale is learned without truth labels; `prior_df` controls its
#' prespecified strength.
#'
#' @inheritParams bb_gene_partial_pool
#' @param prior_df Positive prior degrees of freedom controlling moderation.
#' @return One row per testable gene with moderated heterogeneity and the same
#'   guide-agreement and influence diagnostics as partial pooling.
#' @export
bb_gene_eb_moderate <- function(result, control = NULL, min_guides = 2L,
                                prior_df = 4,
                                alpha = 0.05,
                                min_control_genes = 10L,
                                min_scale = 1) {
  if (length(prior_df) != 1L || !is.finite(prior_df) || prior_df <= 0) {
    .bb_stop("`prior_df` must be one positive finite number.")
  }
  inputs <- .bb_gene_pooling_inputs(result, control, min_guides)
  components <- .bb_gene_pooling_components(inputs)
  total_c <- sum(components$heterogeneity_c)
  prior_tau2 <- if (is.finite(total_c) && total_c > 0) {
    max(
      0,
      sum(
        components$heterogeneity_q - components$heterogeneity_df
      ) / total_c
    )
  } else {
    0
  }
  moderated_tau2 <- (
    components$heterogeneity_df * components$raw_tau2 +
      prior_df * prior_tau2
  ) / (components$heterogeneity_df + prior_df)
  pooled <- .bb_finish_gene_pooling(
    inputs,
    components,
    tau2 = moderated_tau2,
    method = "empirical_bayes",
    alpha = alpha,
    min_control_genes = min_control_genes,
    min_scale = min_scale
  )
  attr(pooled, "prior_tau2") <- prior_tau2
  attr(pooled, "prior_df") <- prior_df
  attr(pooled, "heterogeneity_estimator") <-
    "DerSimonian-Laird moderated toward pooled excess-Q prior"
  attr(pooled, "null_assumption") <- paste(
    "Empirical-Bayes moderated guide random-effects model;",
    "not biological-replicate inference."
  )
  pooled
}
