test_that("generate_blobs returns requested size and variables", {
  dat <- generate_blobs(n = 300, n_blobs = 3, sizes_blobs = c(1, 2, 3), seed = 1)

  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 300)
  expect_equal(names(dat), c("x", "y"))
  expect_true(all(dat$x >= 0 & dat$x <= 1))
  expect_true(all(dat$y >= 0 & dat$y <= 1))
})

test_that("generate_blobs defaults to 100 observations", {
  dat <- generate_blobs(seed = 1)

  expect_equal(nrow(dat), 100)
})

test_that("generate_blobs can return bounded unscaled coordinates", {
  dat <- generate_blobs(
    n = 100,
    n_blobs = 2,
    sizes_blobs = c(60, 40),
    spreads_blobs = c(0.1, 0.1),
    seed = 123,
    rescale = FALSE
  )

  expect_equal(names(dat), c("x", "y"))
  expect_true(all(dat$x >= 0 & dat$x <= 1))
  expect_true(all(dat$y >= 0 & dat$y <= 1))
})

test_that("EM returns an equal-size partition for even data", {
  dat <- generate_blobs(n = 100, seed = 1)
  fit <- esbg(dat, vars = c("x", "y"), n_starts = 5, cluster_method = "em")

  expect_s3_class(fit, "esbg")
  expect_equal(as.integer(table(fit$groups)), c(50, 50))
  expect_true(fit$scores$w1 >= 0)
  expect_true(fit$scores$w2 >= 0)
  expect_true(fit$scores$W >= 0)
  expect_equal(fit$scores$W, mean(c(fit$scores$w1, fit$scores$w2)))
  expect_true(fit$scores$B >= 0)
  expect_equal(fit$scores$P, fit$scores$B / (fit$scores$W + 1))
})

test_that("known 2D minimum polarization case has P equal to zero", {
  dat <- data.frame(
    x = c(0, 0, 0, 0),
    y = c(0, 0, 0, 0)
  )

  fit <- esbg(dat, vars = c("x", "y"), n_starts = 5, delta = 1)

  expect_equal(as.integer(table(fit$groups)), c(2, 2))
  expect_equal(fit$scores$w1, 0)
  expect_equal(fit$scores$w2, 0)
  expect_equal(fit$scores$W, 0)
  expect_equal(fit$scores$B, 0)
  expect_equal(fit$scores$P, 0)
})

test_that("known 2D perfect separation case has maximum normalized P", {
  dat <- data.frame(
    x = c(0, 0, 1, 1),
    y = c(0, 0, 1, 1)
  )

  fit <- esbg(dat, vars = c("x", "y"), n_starts = 5, delta = 2)

  expect_equal(as.integer(table(fit$groups)), c(2, 2))
  expect_equal(fit$scores$w1, 0)
  expect_equal(fit$scores$w2, 0)
  expect_equal(fit$scores$W, 0)
  expect_equal(fit$scores$B, 2)
  expect_equal(fit$scores$P, 1)
})

test_that("odd_method remove leaves one observation unassigned", {
  dat <- generate_blobs(n = 101, seed = 1)
  fit <- esbg(dat, vars = c("x", "y"), n_starts = 5, odd_method = "remove")

  expect_equal(sum(is.na(fit$groups)), 1)
  expect_equal(as.integer(table(fit$groups)), c(50, 50))
  expect_equal(fit$odd$method, "remove")
})

test_that("odd_method min assigns all observations", {
  dat <- generate_blobs(n = 101, seed = 1)
  fit <- esbg(dat, vars = c("x", "y"), n_starts = 5, odd_method = "min")

  expect_equal(sum(is.na(fit$groups)), 0)
  expect_equal(sort(as.integer(table(fit$groups))), c(50, 51))
  expect_equal(fit$odd$method, "min")
  expect_true(fit$odd$selected_assignment %in% c("1", "2"))
})

test_that("vars are inferred while id is excluded", {
  dat <- generate_blobs(n = 100, seed = 1, include_id = TRUE)
  fit <- esbg(dat, n_starts = 5)

  expect_equal(fit$vars, c("x", "y"))
})

test_that("same_partition ignores label switching", {
  groups_a <- factor(c(1, 1, 2, 2))
  groups_b <- factor(c(2, 2, 1, 1))

  expect_true(same_partition(groups_a, groups_b))
})

test_that("plot works only for two-dimensional fits", {
  dat <- generate_blobs(n = 100, seed = 1)
  fit_2d <- esbg(dat, vars = c("x", "y"), n_starts = 5)

  expect_s3_class(plot(fit_2d), "ggplot")

  dat$z <- dat$x + dat$y
  fit_3d <- esbg(dat, vars = c("x", "y", "z"), n_starts = 5)
  expect_error(plot(fit_3d), "exactly two variables")
})
