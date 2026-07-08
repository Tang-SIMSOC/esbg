#' Equal-Size Binary Grouping for polarization measurement
#'
#' Estimate an Equal-Size Binary Grouping (ESBG) partition for numeric opinion
#' data and compute within-group dispersion and between-group separation.
#'
#' @param data A data frame containing individual-level observations.
#' @param vars Optional character vector of numeric columns defining the opinion
#'   space. If `NULL`, numeric columns are inferred from `data`, excluding common
#'   identifier columns such as `id`.
#' @param exclude_vars Character vector of columns to exclude when `vars = NULL`.
#'   Ignored when `vars` is supplied.
#' @param n_starts Number of random restarts. For `cluster_method = "em"`, each
#'   restart begins from a random equal-size partition. For
#'   `cluster_method = "centroid"`, the input order is randomly shuffled before
#'   each call to `anticlust::balanced_clustering()`. Across starts, ESBG keeps
#'   the partition with the smallest within-group sum of squares, equivalent to
#'   the smallest `W`. The starts are not averaged.
#' @param max_iter Maximum number of EM iterations per random start. Used only
#'   when `cluster_method = "em"`; ignored when `cluster_method = "centroid"`.
#' @param tol Convergence tolerance for the within-group sum of squares. Used
#'   only when `cluster_method = "em"`; ignored when
#'   `cluster_method = "centroid"`.
#' @param delta Positive normalizing parameter used to compute the polarization
#'   index `P = (1 / delta) * B / (W + 1)`. The default is `1`. For comparing
#'   multiple datasets, provide the same `delta` for each call.
#' @param cluster_method Clustering backend. `"em"` uses the constrained
#'   k-means/EM implementation described in Tang et al. (2022). `"centroid"`
#'   uses `anticlust::balanced_clustering(..., method = "centroid")`.
#' @param scale Logical. If `TRUE`, selected variables are standardized before
#'   ESBG is estimated.
#' @param na_rm Logical. If `TRUE`, rows with missing values in `vars` are
#'   omitted from estimation. If `FALSE`, missing values trigger an error.
#' @param odd_method How to handle an odd number of complete observations.
#'   `"error"` stops with an error. `"remove"` removes the observation closest
#'   to the component-wise median and reports ESBG for the remaining even sample.
#'   `"average"` removes that observation, estimates ESBG on the remaining even
#'   sample, evaluates both ways of assigning the median observation to one of
#'   the two groups, and averages the resulting scores. `"min"` reports the
#'   assignment of the median observation that gives the smaller within-group sum
#'   of squares. The default is `"average"`.
#'
#' @return An object of class `"esbg"` containing group assignments, centroids,
#'   scores, and estimation metadata.
#' @export
#'
#' @examples
#' data <- generate_blobs(n = 100)
#' fit <- esbg(data, vars = c("x", "y"), n_starts = 20)
#' fit_all_numeric <- esbg(data, n_starts = 20)
#' summary(fit)
#' plot(fit)
esbg <- function(data, vars = NULL, n_starts = 50, max_iter = 100, tol = 1e-8,
                 delta = 1, scale = FALSE, na_rm = FALSE, exclude_vars = c("id"),
                 cluster_method = c("em", "centroid"),
                 odd_method = "average") {
  cluster_method <- match.arg(cluster_method)
  odd_method <- match.arg(odd_method, c("error", "remove", "average", "min"))
  prepared <- prepare_esbg_data(
    data,
    vars,
    scale = scale,
    na_rm = na_rm,
    exclude_vars = exclude_vars
  )
  coords <- prepared$coords
  vars <- prepared$vars

  partition <- run_esbg_partition(
    coords,
    n_starts = n_starts,
    max_iter = max_iter,
    tol = tol,
    cluster_method = cluster_method,
    odd_method = odd_method
  )
  if (!is.null(partition$odd$median_index)) {
    partition$odd$median_row <- prepared$row_index[partition$odd$median_index]
  }
  groups_complete <- partition$groups
  scores <- evaluate_esbg_partition(coords, partition)
  if (!is.null(scores$selected_assignment)) {
    groups_complete[partition$odd$median_index] <- scores$selected_assignment
    partition$groups <- groups_complete
    partition$odd$selected_assignment <- scores$selected_assignment
  }

  groups <- rep(NA_integer_, nrow(data))
  groups[prepared$row_index] <- as.integer(groups_complete)
  groups <- factor(groups, levels = c(1, 2), labels = c("1", "2"))
  delta <- resolve_delta(delta, coords)
  scores$P <- polarization_index(W = scores$W, B = scores$B, delta = delta)

  result <- list(
    groups = groups,
    groups_complete = groups_complete,
    centroids = scores$centroids,
    scores = list(w1 = scores$w1, w2 = scores$w2, W = scores$W, B = scores$B, P = scores$P),
    delta = delta,
    data = data,
    coords = coords,
    vars = vars,
    row_index = prepared$row_index,
    scale = scale,
    center = prepared$center,
    scale_values = prepared$scale_values,
    n_starts = n_starts,
    max_iter = max_iter,
    tol = tol,
    cluster_method = cluster_method,
    odd = partition$odd,
    call = match.call()
  )

  class(result) <- "esbg"
  result
}

#' Print an ESBG fit
#'
#' @param x An object returned by [esbg()].
#' @param ... Unused.
#' @export
print.esbg <- function(x, ...) {
  cat("Equal-Size Binary Grouping (ESBG)\n")
  cat("Variables:", paste(x$vars, collapse = ", "), "\n")
  cat("Complete observations:", length(x$groups_complete), "\n")
  cat("Group sizes:", paste(as.integer(table(x$groups_complete)), collapse = " / "), "\n")
  if (!is.null(x$odd$median_index)) {
    cat("Odd-n handling:", odd_method_label(x$odd$method), "\n")
  }
  cat(sprintf("w1: %.4f\n", x$scores$w1))
  cat(sprintf("w2: %.4f\n", x$scores$w2))
  cat(sprintf("W: %.4f\n", x$scores$W))
  cat(sprintf("B: %.4f\n", x$scores$B))
  cat(sprintf("P: %.4f (delta = %.4f)\n", x$scores$P, x$delta))
  invisible(x)
}

#' Summarize an ESBG fit
#'
#' @param object An object returned by [esbg()].
#' @param ... Unused.
#' @export
summary.esbg <- function(object, ...) {
  out <- list(
    call = object$call,
    vars = object$vars,
    n = length(object$groups_complete),
    group_sizes = table(object$groups_complete),
    centroids = object$centroids,
    scores = object$scores,
    delta = object$delta,
    odd = object$odd
  )
  class(out) <- "summary.esbg"
  out
}

#' Print an ESBG summary
#'
#' @param x An object returned by [summary.esbg()].
#' @param ... Unused.
#' @export
print.summary.esbg <- function(x, ...) {
  cat("Summary of Equal-Size Binary Grouping (ESBG)\n\n")
  cat("Variables:", paste(x$vars, collapse = ", "), "\n")
  cat("Observations:", x$n, "\n")
  cat("Group sizes:", paste(as.integer(x$group_sizes), collapse = " / "), "\n\n")
  if (!is.null(x$odd$median_index)) {
    cat("Odd-n handling:\n")
    cat(" ", odd_method_label(x$odd$method), "\n\n")
  }
  cat("Scores:\n")
  cat(sprintf("  w1, group 1 within-group heterogeneity: %.4f\n", x$scores$w1))
  cat(sprintf("  w2, group 2 within-group heterogeneity: %.4f\n", x$scores$w2))
  cat(sprintf("  W, overall within-group heterogeneity: %.4f\n", x$scores$W))
  cat(sprintf("  B, between-group heterogeneity: %.4f\n", x$scores$B))
  cat(sprintf("  P, polarization index: %.4f (delta = %.4f)\n\n", x$scores$P, x$delta))
  cat("Centroids:\n")
  print(x$centroids)
  invisible(x)
}

#' Plot an ESBG fit
#'
#' The plot uses a Figure 11-inspired style: ESBG groups are shown in blue and
#' orange, and group centroids are shown as lighter triangles.
#'
#' @param x An object returned by [esbg()].
#' @param ... Unused.
#'
#' @return A ggplot object.
#' @export
plot.esbg <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  if (length(x$vars) != 2) {
    stop(
      "`plot()` is available only for ESBG fits with exactly two variables. ",
      "This fit has ", length(x$vars), " variables.",
      call. = FALSE
    )
  }

  plot_data <- x$data[x$row_index, , drop = FALSE]
  plot_data$.esbg_group <- x$groups_complete
  x_var <- x$vars[1]
  y_var <- x$vars[2]
  plot_data$.esbg_x <- plot_data[[x_var]]
  plot_data$.esbg_y <- plot_data[[y_var]]

  centroids <- x$centroids
  if (isTRUE(x$scale)) {
    centroids <- sweep(centroids, 2, x$scale_values, FUN = "*")
    centroids <- sweep(centroids, 2, x$center, FUN = "+")
  }

  centroid_data <- data.frame(
    .esbg_x = centroids[, x_var],
    .esbg_y = centroids[, y_var],
    .esbg_group = factor(rownames(centroids), levels = c("1", "2")),
    row.names = NULL
  )

  display_order <- rownames(centroids)[do.call(order, as.data.frame(centroids[, c(x_var, y_var), drop = FALSE]))]
  display_labels <- stats::setNames(c("1", "2"), display_order)
  plot_data$.esbg_display_group <- factor(
    display_labels[as.character(plot_data$.esbg_group)],
    levels = c("1", "2")
  )
  centroid_data$.esbg_display_group <- factor(
    display_labels[as.character(centroid_data$.esbg_group)],
    levels = c("1", "2")
  )

  group_colors <- c("1" = "#1f77b4", "2" = "#ff7f0e")
  centroid_colors <- c("1" = "#8fc4e6", "2" = "#ffbd73")
  centroid_labels <- c(
    "1" = "Centroid (group 1)",
    "2" = "Centroid (group 2)"
  )
  internal_w <- c("1" = x$scores$w1, "2" = x$scores$w2)
  display_w <- stats::setNames(internal_w[display_order], c("1", "2"))

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .esbg_x, y = .esbg_y)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = .esbg_display_group),
      shape = 21,
      color = "black",
      stroke = 0.35,
      alpha = 0.90,
      size = 2.2
    ) +
    ggplot2::geom_point(
      data = centroid_data[centroid_data$.esbg_display_group == "1", , drop = FALSE],
      ggplot2::aes(x = .esbg_x, y = .esbg_y, shape = centroid_labels[["1"]]),
      fill = centroid_colors[["1"]],
      color = "black",
      stroke = 0.35,
      size = 2.9,
      show.legend = TRUE
    ) +
    ggplot2::geom_point(
      data = centroid_data[centroid_data$.esbg_display_group == "2", , drop = FALSE],
      ggplot2::aes(x = .esbg_x, y = .esbg_y, shape = centroid_labels[["2"]]),
      fill = centroid_colors[["2"]],
      color = "black",
      stroke = 0.35,
      size = 2.9,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      name = "Group",
      values = group_colors,
      breaks = c("1", "2"),
      na.value = "grey70",
      guide = ggplot2::guide_legend(
        override.aes = list(
          shape = 21,
          fill = unname(group_colors),
          color = "black",
          size = 3
        )
      )
    ) +
    ggplot2::scale_shape_manual(
      name = "Centroids",
      values = c("Centroid (group 1)" = 24, "Centroid (group 2)" = 24),
      guide = ggplot2::guide_legend(
        override.aes = list(
          fill = unname(centroid_colors),
          color = "black",
          size = 2.9
        )
      )
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.4),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(
      title = "Equal-Size Binary Grouping (ESBG)",
      subtitle = sprintf("w1: %.3f, w2: %.3f, W: %.3f, B: %.3f, P: %.3f",
                         display_w[["1"]], display_w[["2"]], x$scores$W, x$scores$B, x$scores$P),
      x = x_var,
      y = y_var,
      shape = "Centroids"
    )
}

#' Compare two binary partitions
#'
#' Checks whether two binary group assignments represent the same partition,
#' allowing the group labels to be swapped. Missing values are ignored by default,
#' which is useful when `odd_method = "average"` or `odd_method = "remove"` leaves
#' one median observation unassigned.
#'
#' @param x A vector of group assignments, or an object returned by [esbg()].
#' @param y A vector of group assignments, or an object returned by [esbg()].
#' @param na_rm Logical. If `TRUE`, compare only observations assigned in both
#'   partitions. If `FALSE`, missingness must also match.
#'
#' @return `TRUE` if the partitions are identical up to label switching;
#'   otherwise `FALSE`.
#' @export
#'
#' @examples
#' same_partition(c(1, 1, 2, 2), c(2, 2, 1, 1))
same_partition <- function(x, y, na_rm = TRUE) {
  x <- extract_groups(x)
  y <- extract_groups(y)

  if (length(x) != length(y)) {
    stop("`x` and `y` must have the same length.", call. = FALSE)
  }

  if (na_rm) {
    keep <- !is.na(x) & !is.na(y)
    x <- x[keep]
    y <- y[keep]
  } else if (!identical(is.na(x), is.na(y))) {
    return(FALSE)
  } else {
    keep <- !is.na(x)
    x <- x[keep]
    y <- y[keep]
  }

  if (length(x) == 0) {
    return(TRUE)
  }

  x <- factor(x)
  y <- factor(y)

  if (length(levels(droplevels(x))) != 2 || length(levels(droplevels(y))) != 2) {
    stop("`x` and `y` must each contain exactly two non-missing groups.", call. = FALSE)
  }

  rel_x <- outer(as.character(x), as.character(x), FUN = "==")
  rel_y <- outer(as.character(y), as.character(y), FUN = "==")
  identical(rel_x, rel_y)
}

extract_groups <- function(x) {
  if (inherits(x, "esbg")) {
    return(x$groups)
  }
  x
}

polarization_index <- function(W, B = NULL, delta = NULL) {
  if (inherits(W, "esbg")) {
    fit <- W
    W <- fit$scores$W
    B <- fit$scores$B
    if (is.null(delta)) {
      delta <- fit$delta
    }
  } else {
    if (is.null(B)) {
      stop("`B` must be supplied when `W` is numeric.", call. = FALSE)
    }
  }

  if (!is.numeric(W) || length(W) != 1 || is.na(W) || W < 0) {
    stop("`W` must be a non-negative numeric value.", call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1 || is.na(B) || B < 0) {
    stop("`B` must be a non-negative numeric value.", call. = FALSE)
  }
  if (!is.numeric(delta) || length(delta) != 1 || is.na(delta) || delta <= 0) {
    stop("`delta` must be a positive numeric value.", call. = FALSE)
  }

  (1 / delta) * (B / (W + 1))
}

resolve_delta <- function(delta = 1, coords) {
  if (!is.numeric(delta) || length(delta) != 1 || is.na(delta) || delta <= 0) {
    stop("`delta` must be a positive numeric value.", call. = FALSE)
  }

  delta
}

max_squared_distance <- function(x) {
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("`x` must be numeric.", call. = FALSE)
  }
  if (nrow(x) < 2) {
    stop("`x` must contain at least two rows.", call. = FALSE)
  }

  max(stats::dist(x)^2)
}

#' Generate synthetic blob data
#'
#' Creates a two-dimensional synthetic dataset with one or more natural blobs.
#' This is useful for demonstrating ESBG on data with visible spatial structure.
#'
#' @param n Number of observations.
#' @param n_blobs Number of blobs.
#' @param sizes_blobs Optional blob sizes. If `NULL`, blob sizes are as equal as
#'   possible. If `sizes_blobs` is a positive integer vector that sums to `n`, it
#'   is used as exact blob counts. Otherwise, `sizes_blobs` is interpreted as
#'   relative weights. For example, with `n = 300`, `sizes_blobs = c(1, 2, 3)`
#'   produces blob sizes
#'   `c(50, 100, 150)`.
#' @param spreads_blobs Optional blob spreads, used as the standard deviation
#'   for each blob. If `NULL`, spreads increase evenly from `0.4` to `0.8`.
#' @param seed Optional random seed.
#' @param include_id Logical. If `TRUE`, include an `id` column.
#' @param rescale Logical. If `TRUE`, generate unbounded normal blobs and
#'   min-max rescale `x` and `y` to `[0, 1]`. If `FALSE`, generate bounded
#'   normal blobs directly inside `[0, 1]`.
#'
#' @return A data frame with columns `x` and `y`, plus `id` when
#'   `include_id = TRUE`.
#' @export
generate_blobs <- function(n = 100, n_blobs = 2, sizes_blobs = NULL, spreads_blobs = NULL, seed = NULL,
                           include_id = FALSE, rescale = TRUE) {
  if (!is.logical(include_id) || length(include_id) != 1 || is.na(include_id)) {
    stop("`include_id` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (!is.logical(rescale) || length(rescale) != 1 || is.na(rescale)) {
    stop("`rescale` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL
    on.exit({
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  blob_sizes <- resolve_blob_sizes(n = n, n_blobs = n_blobs, sizes_blobs = sizes_blobs)
  centers <- blob_centers(n_blobs, unit_square = !rescale)
  blob_sds <- resolve_blob_spread(n_blobs = n_blobs, spreads_blobs = spreads_blobs)

  x <- numeric(n)
  y <- numeric(n)
  start <- 1

  for (j in seq_len(n_blobs)) {
    end <- start + blob_sizes[j] - 1
    idx <- start:end
    if (rescale) {
      x[idx] <- stats::rnorm(blob_sizes[j], mean = centers[j, 1], sd = blob_sds[j])
      y[idx] <- stats::rnorm(blob_sizes[j], mean = centers[j, 2], sd = blob_sds[j])
    } else {
      x[idx] <- bounded_normal_01(blob_sizes[j], mean = centers[j, 1], sd = blob_sds[j])
      y[idx] <- bounded_normal_01(blob_sizes[j], mean = centers[j, 2], sd = blob_sds[j])
    }
    start <- end + 1
  }

  if (rescale) {
    x <- min_max_scale(x)
    y <- min_max_scale(y)
  }

  out <- data.frame(
    x = x,
    y = y
  )

  if (include_id) {
    out <- data.frame(id = seq_len(n), out)
  }

  out
}

resolve_blob_sizes <- function(n, n_blobs, sizes_blobs = NULL) {
  if (length(n) != 1 || is.na(n) || n < 2 || n != as.integer(n)) {
    stop("`n` must be an integer greater than or equal to 2.", call. = FALSE)
  }
  if (length(n_blobs) != 1 || is.na(n_blobs) || n_blobs < 1 || n_blobs != as.integer(n_blobs)) {
    stop("`n_blobs` must be a positive integer.", call. = FALSE)
  }
  if (n_blobs > n) {
    stop("`n_blobs` cannot be larger than `n`.", call. = FALSE)
  }

  if (is.null(sizes_blobs)) {
    base <- rep(floor(n / n_blobs), n_blobs)
    base[seq_len(n %% n_blobs)] <- base[seq_len(n %% n_blobs)] + 1
    return(base)
  }

  if (!is.numeric(sizes_blobs) || length(sizes_blobs) != n_blobs || any(is.na(sizes_blobs)) || any(sizes_blobs <= 0)) {
    stop("`sizes_blobs` must be a positive numeric vector with length equal to `n_blobs`.", call. = FALSE)
  }

  if (all(sizes_blobs == as.integer(sizes_blobs)) && sum(sizes_blobs) == n) {
    return(as.integer(sizes_blobs))
  }

  weights <- sizes_blobs / sum(sizes_blobs)
  raw_sizes <- weights * n
  sizes <- floor(raw_sizes)
  remainder <- n - sum(sizes)
  if (remainder > 0) {
    add_to <- order(raw_sizes - sizes, decreasing = TRUE)[seq_len(remainder)]
    sizes[add_to] <- sizes[add_to] + 1
  }

  if (any(sizes == 0)) {
    stop("`sizes_blobs` creates at least one empty blob. Use fewer blobs or larger weights.", call. = FALSE)
  }

  as.integer(sizes)
}

resolve_blob_spread <- function(n_blobs, spreads_blobs = NULL) {
  if (is.null(spreads_blobs)) {
    return(seq(0.4, 0.8, length.out = n_blobs))
  }

  if (!is.numeric(spreads_blobs) || any(is.na(spreads_blobs)) || any(spreads_blobs <= 0)) {
    stop("`spreads_blobs` must contain positive numeric values.", call. = FALSE)
  }

  if (length(spreads_blobs) == 1) {
    return(rep(spreads_blobs, n_blobs))
  }

  if (length(spreads_blobs) != n_blobs) {
    stop("`spreads_blobs` must have length 1 or length equal to `n_blobs`.", call. = FALSE)
  }

  spreads_blobs
}

blob_centers <- function(n_blobs, unit_square = FALSE) {
  if (unit_square) {
    if (n_blobs == 1) {
      return(matrix(c(0.5, 0.5), ncol = 2))
    }
    if (n_blobs == 2) {
      return(matrix(c(0.8, 0.8, 0.2, 0.2), ncol = 2, byrow = TRUE))
    }

    theta <- seq(0, 2 * pi, length.out = n_blobs + 1)[seq_len(n_blobs)]
    return(cbind(0.5 + 0.35 * cos(theta), 0.5 + 0.35 * sin(theta)))
  }

  if (n_blobs == 1) {
    return(matrix(c(0, 0), ncol = 2))
  }
  if (n_blobs == 2) {
    return(matrix(c(4, 0, 0, -4), ncol = 2, byrow = TRUE))
  }

  theta <- seq(0, 2 * pi, length.out = n_blobs + 1)[seq_len(n_blobs)]
  cbind(4 * cos(theta), 4 * sin(theta))
}

bounded_normal_01 <- function(n, mean, sd) {
  lower_p <- stats::pnorm(0, mean = mean, sd = sd)
  upper_p <- stats::pnorm(1, mean = mean, sd = sd)
  stats::qnorm(stats::runif(n, min = lower_p, max = upper_p), mean = mean, sd = sd)
}

prepare_esbg_data <- function(data, vars = NULL, scale = FALSE, na_rm = FALSE,
                              exclude_vars = c("id")) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  vars <- resolve_esbg_vars(data, vars, exclude_vars = exclude_vars)

  if (is.null(vars)) {
    stop("Internal error: variable selection failed.", call. = FALSE)
  }
  if (!is.character(vars) || length(vars) < 1) {
    stop("`vars` must be a character vector of column names.", call. = FALSE)
  }
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop("These variables are not in `data`: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  non_numeric <- vars[!vapply(data[vars], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("These variables are not numeric: ", paste(non_numeric, collapse = ", "), call. = FALSE)
  }

  complete <- stats::complete.cases(data[vars])
  if (!all(complete) && !na_rm) {
    stop("Missing values found in `vars`. Remove/impute them or set `na_rm = TRUE`.", call. = FALSE)
  }

  row_index <- which(complete)
  coords <- as.matrix(data[row_index, vars, drop = FALSE])
  if (nrow(coords) < 2) {
    stop("ESBG requires at least two complete observations.", call. = FALSE)
  }

  center <- rep(0, ncol(coords))
  scale_values <- rep(1, ncol(coords))
  names(center) <- names(scale_values) <- vars

  if (scale) {
    center <- colMeans(coords)
    scale_values <- apply(coords, 2, stats::sd)
    zero_variance <- names(scale_values)[scale_values == 0 | is.na(scale_values)]
    if (length(zero_variance) > 0) {
      stop("Cannot scale zero-variance variables: ", paste(zero_variance, collapse = ", "), call. = FALSE)
    }
    coords <- sweep(coords, 2, center, FUN = "-")
    coords <- sweep(coords, 2, scale_values, FUN = "/")
  }

  list(
    coords = coords,
    vars = vars,
    row_index = row_index,
    center = center,
    scale_values = scale_values
  )
}

resolve_esbg_vars <- function(data, vars = NULL, exclude_vars = c("id")) {
  if (!is.null(vars)) {
    return(vars)
  }

  numeric_vars <- names(data)[vapply(data, is.numeric, logical(1))]
  excluded_numeric_vars <- intersect(numeric_vars, exclude_vars)
  numeric_vars <- setdiff(numeric_vars, exclude_vars)

  if (length(numeric_vars) < 1) {
    stop(
      "`vars` is NULL, but `data` has no usable numeric columns after excluding identifier columns.",
      call. = FALSE
    )
  }

  msg <- paste0(
    "`vars` is NULL; using numeric columns: ",
    paste(numeric_vars, collapse = ", ")
  )
  if (length(excluded_numeric_vars) > 0) {
    msg <- paste0(
      msg,
      ". Excluded: ",
      paste(excluded_numeric_vars, collapse = ", ")
    )
  }
  message(msg, ". For reproducible analysis, consider setting `vars` explicitly.")

  numeric_vars
}

run_esbg_partition <- function(coords, n_starts = 50, max_iter = 100, tol = 1e-8,
                               cluster_method = c("em", "centroid"),
                               odd_method = "average") {
  cluster_method <- match.arg(cluster_method)
  if (length(n_starts) != 1 || n_starts < 1 || n_starts != as.integer(n_starts)) {
    stop("`n_starts` must be a positive integer.", call. = FALSE)
  }
  if (length(max_iter) != 1 || max_iter < 1 || max_iter != as.integer(max_iter)) {
    stop("`max_iter` must be a positive integer.", call. = FALSE)
  }
  if (length(tol) != 1 || is.na(tol) || tol < 0) {
    stop("`tol` must be a non-negative number.", call. = FALSE)
  }
  odd_method <- match.arg(odd_method, c("error", "remove", "average", "min"))

  if (nrow(coords) %% 2 == 0) {
    groups <- run_even_esbg_partition(
      coords,
      n_starts = n_starts,
      max_iter = max_iter,
      tol = tol,
      cluster_method = cluster_method
    )
    return(list(
      groups = groups,
      alternatives = NULL,
      odd = list(method = "none", median_index = NULL)
    ))
  }

  if (odd_method == "error") {
    stop(
      "ESBG requires equal-sized groups for the selected clustering method. ",
      "Use an even number of complete observations or choose another `odd_method` option.",
      call. = FALSE
    )
  }

  median_index <- closest_to_component_median(coords)
  even_coords <- coords[-median_index, , drop = FALSE]
  even_groups <- run_even_esbg_partition(
    even_coords,
    n_starts = n_starts,
    max_iter = max_iter,
    tol = tol,
    cluster_method = cluster_method
  )

  groups <- rep(NA_character_, nrow(coords))
  groups[-median_index] <- as.character(even_groups)
  groups <- factor(groups, levels = c("1", "2"))

  group_1 <- groups
  group_1[median_index] <- "1"
  group_2 <- groups
  group_2[median_index] <- "2"

  list(
    groups = groups,
    alternatives = list(
      median_to_group_1 = factor(group_1, levels = c("1", "2")),
      median_to_group_2 = factor(group_2, levels = c("1", "2"))
    ),
    odd = list(method = odd_method, median_index = median_index)
  )
}

run_even_esbg_partition <- function(coords, n_starts = 50, max_iter = 100, tol = 1e-8,
                                    cluster_method = c("em", "centroid")) {
  cluster_method <- match.arg(cluster_method)
  if (cluster_method == "centroid") {
    return(run_anticlust_centroid_partition(coords, n_starts = n_starts))
  }

  if (cluster_method == "em") {
    best_groups <- NULL
    best_score <- Inf

    for (i in seq_len(n_starts)) {
      groups <- random_equal_partition(nrow(coords))
      previous_score <- Inf

      for (iter in seq_len(max_iter)) {
        centroids <- compute_centroids(coords, groups)
        updated_groups <- constrained_centroid_assignment(coords, centroids)
        current_score <- within_group_ss(coords, updated_groups)

        if (identical(as.character(updated_groups), as.character(groups))) {
          groups <- updated_groups
          break
        }
        if (abs(previous_score - current_score) <= tol) {
          groups <- updated_groups
          break
        }

        groups <- updated_groups
        previous_score <- current_score
      }

      current_score <- within_group_ss(coords, groups)
      if (current_score < best_score) {
        best_score <- current_score
        best_groups <- groups
      }
    }

    return(best_groups)
  }

  stop("Unknown `cluster_method`: ", cluster_method, call. = FALSE)
}

run_anticlust_centroid_partition <- function(coords, n_starts = 50) {
  if (!requireNamespace("anticlust", quietly = TRUE)) {
    stop(
      "Package 'anticlust' is required for `cluster_method = \"centroid\"`.",
      call. = FALSE
    )
  }

  best_groups <- NULL
  best_score <- Inf

  for (i in seq_len(n_starts)) {
    shuffled_index <- sample(nrow(coords))
    shuffled_coords <- coords[shuffled_index, , drop = FALSE]
    shuffled_groups <- anticlust::balanced_clustering(
      shuffled_coords,
      K = 2,
      method = "centroid"
    )

    groups <- integer(nrow(coords))
    groups[shuffled_index] <- shuffled_groups
    groups <- factor(groups, levels = c(1, 2), labels = c("1", "2"))

    current_score <- within_group_ss(coords, groups)
    if (current_score < best_score) {
      best_score <- current_score
      best_groups <- groups
    }
  }

  best_groups
}

random_equal_partition <- function(n) {
  if (n %% 2 != 0) {
    stop("`n` must be even for an equal-size binary partition.", call. = FALSE)
  }

  shuffled_index <- sample(n)
  groups <- integer(n)
  groups[shuffled_index[seq_len(n / 2)]] <- 1L
  groups[shuffled_index[(n / 2 + 1):n]] <- 2L
  factor(groups, levels = c(1, 2), labels = c("1", "2"))
}

compute_centroids <- function(coords, groups) {
  centroids <- vapply(levels(groups), function(g) {
    colMeans(coords[groups == g, , drop = FALSE])
  }, numeric(ncol(coords)))

  t(centroids)
}

constrained_centroid_assignment <- function(coords, centroids) {
  d1 <- rowSums((sweep(coords, 2, centroids["1", ], FUN = "-"))^2)
  d2 <- rowSums((sweep(coords, 2, centroids["2", ], FUN = "-"))^2)
  delta <- abs(d1 - d2)
  groups <- ifelse(d1 <= d2, "1", "2")

  target_size <- nrow(coords) / 2
  n_group_1 <- sum(groups == "1")

  if (n_group_1 > target_size) {
    move <- order(delta[groups == "1"])[seq_len(n_group_1 - target_size)]
    move_index <- which(groups == "1")[move]
    groups[move_index] <- "2"
  } else if (n_group_1 < target_size) {
    move <- order(delta[groups == "2"])[seq_len(target_size - n_group_1)]
    move_index <- which(groups == "2")[move]
    groups[move_index] <- "1"
  }

  factor(groups, levels = c("1", "2"))
}

evaluate_esbg_partition <- function(coords, partition) {
  if (is.null(partition$alternatives)) {
    return(evaluate_esbg_coords(coords, partition$groups))
  }

  if (partition$odd$method == "remove") {
    kept <- !is.na(partition$groups)
    return(evaluate_esbg_coords(coords[kept, , drop = FALSE], partition$groups[kept]))
  }

  score_1 <- evaluate_esbg_coords(coords, partition$alternatives$median_to_group_1)
  score_2 <- evaluate_esbg_coords(coords, partition$alternatives$median_to_group_2)

  if (partition$odd$method == "min") {
    if (score_1$W <= score_2$W) {
      return(c(score_1, list(selected_assignment = "1")))
    }
    return(c(score_2, list(selected_assignment = "2")))
  }

  list(
    W = mean(c(score_1$W, score_2$W)),
    w1 = mean(c(score_1$w1, score_2$w1)),
    w2 = mean(c(score_1$w2, score_2$w2)),
    B = mean(c(score_1$B, score_2$B)),
    centroids = (score_1$centroids + score_2$centroids) / 2,
    alternatives = list(
      median_to_group_1 = score_1,
      median_to_group_2 = score_2
    )
  )
}

odd_method_label <- function(method) {
  switch(
    method,
    none = "none",
    remove = "median observation removed before scoring",
    average = "median observation evaluated in both groups; scores averaged",
    min = "median observation assigned to the group with smaller within-group sum of squares",
    method
  )
}

closest_to_component_median <- function(coords) {
  med <- apply(coords, 2, stats::median)
  distances <- rowSums((sweep(coords, 2, med, FUN = "-"))^2)
  unname(which.min(distances))
}

evaluate_esbg_coords <- function(coords, groups) {
  groups <- factor(groups)
  group_levels <- levels(groups)

  if (length(group_levels) != 2) {
    stop("ESBG evaluation requires exactly two groups.", call. = FALSE)
  }

  centroids <- matrix(
    NA_real_,
    nrow = 2,
    ncol = ncol(coords),
    dimnames = list(group_levels, colnames(coords))
  )

  ssw <- 0
  group_w <- numeric(length(group_levels))
  names(group_w) <- group_levels
  for (g in group_levels) {
    pts <- coords[groups == g, , drop = FALSE]
    centroid <- colMeans(pts)
    group_ssw <- sum(rowSums((sweep(pts, 2, centroid, FUN = "-"))^2))
    centroids[g, ] <- centroid
    group_w[g] <- group_ssw / nrow(pts)
    ssw <- ssw + group_ssw
  }

  W <- ssw / nrow(coords)
  B <- sum((centroids[1, ] - centroids[2, ])^2)

  list(w1 = unname(group_w[1]), w2 = unname(group_w[2]), W = W, B = B, centroids = centroids)
}

within_group_ss <- function(coords, groups) {
  sum(vapply(levels(groups), function(g) {
    pts <- coords[groups == g, , drop = FALSE]
    centroid <- colMeans(pts)
    sum(rowSums((sweep(pts, 2, centroid, FUN = "-"))^2))
  }, numeric(1)))
}

min_max_scale <- function(x) {
  range_x <- range(x, na.rm = TRUE)
  if (diff(range_x) == 0) {
    return(rep(0, length(x)))
  }
  (x - range_x[1]) / diff(range_x)
}
