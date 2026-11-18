mock_Qq_2d <- function(X, Qq, Qq.mean = rowMeans(Qq)) {
  mesh <- fmesher::fm_mesh_1d(seq(-pi, pi, length.out = 5), boundary = 'cyclic', degree = 1)
  list(X = X, mesh = mesh, Qq = Qq, Qq.mean = list(mean = Qq.mean))
}

test_that('runs clusters reproduce the toy example of the manuscript', {
  cl <- runs_clusters(c(3L, 5L, 6L, 11L, 12L, 16L), m = 3)
  expect_equal(cl, list(c(3L, 5L, 6L), c(11L, 12L), 16L))
  expect_equal(sapply(cl, function(x) diff(range(x)) + 1), c(4, 2, 1))
  expect_equal(runs_clusters(integer(0), m = 3), list())
  expect_equal(runs_clusters(7L, m = 3), list(7L))
})

test_that('radial maxima are exceedances under a direction-dependent threshold', {
  r <- c(3, 5, 2.5)
  w <- c(0, pi / 2, 0)
  X <- cbind(r * cos(w), r * sin(w))
  fQ <- mock_Qq_2d(X, Qq = cbind(c(2, 2, 2, 6)))
  expect_equal(get_exceedances(fQ), list(c(1L, 3L)))
  expect_equal(get_clusters_runs(fQ, m = 3), list(list(c(1L, 3L))))
  expect_equal(get_radial_maxima(fQ, m = 3), list(1L))
})

test_that('one result per posterior sample of Q_q, or one for its posterior mean', {
  set.seed(3)
  w <- runif(200, -pi, pi)
  r <- rexp(200)
  X <- cbind(r * cos(w), r * sin(w))
  Qq <- cbind(rep(0.5, 4), rep(1, 4), rep(2, 4))
  fQ <- mock_Qq_2d(X, Qq, Qq.mean = rep(3, 4))
  exc <- get_exceedances(fQ)
  expect_length(exc, 3)
  expect_equal(exc[[2]], which(r > 1))
  expect_equal(get_exceedances(fQ, use.mean.Qq = TRUE), list(which(r > 3)))
  rm <- get_radial_maxima(fQ, m = 2)
  expect_length(rm, 3)
  for (i in 1:3) expect_true(all(rm[[i]] %in% exc[[i]]))
  expect_equal(get_radial_maxima(fQ, m = 2, use.mean.Qq = TRUE)[[1]],
               vapply(runs_clusters(which(r > 3), 2), function(x) x[which.max(r[x])], integer(1)))
})

test_that('filtered Z indicators follow the runs definition', {
  exc <- c(3L, 5L, 6L, 11L, 12L, 16L)
  expect_equal(filtered_Z_indicators(exc, m = 3), c(0, 0, 1))
  expect_equal(filtered_Z_indicators(exc, m = 1), c(1, 0, 1, 0, 1))
  expect_equal(filtered_Z_indicators(exc, m = 1), as.numeric(diff(exc) > 1)[1:5])
  expect_length(filtered_Z_indicators(c(2L, 4L), m = 3), 0)
})

test_that('fit_theta recovers a known directional extremal index in 2d', {
  set.seed(11)
  theta_true <- function(w) 0.2 + 0.6 * plogis(3 * cos(w - 1))
  make_Z <- function(n) { w <- runif(n, -pi, pi); cbind(rbinom(n, 1, theta_true(w)), w) }
  filtered_Z <- list(make_Z(1500), make_Z(1500))
  opts <- list(mesh.knots.2d = seq(-pi, pi, by = pi / 32), zeta.prior.range = c(0.2, 0.8),
               zeta.prior.sigma = c(10, 0.8), N.theta = 20, seed = 0L)
  th <- suppressMessages(fit_theta(filtered_Z, opts))
  expect_length(th$theta, 2)
  expect_true(all(vapply(th$theta, function(x) is.matrix(x) && ncol(x) == 20 && nrow(x) == th$mesh$n, logical(1))))
  expect_true(all(th$theta[[1]] > 0 & th$theta[[1]] < 1))
  expect_lt(max(abs(rowMeans(th$theta[[1]]) - theta_true(th$mesh$loc))), 0.15)
  tf <- tempfile(fileext = '.pdf'); grDevices::pdf(tf)
  on.exit({grDevices::dev.off(); unlink(tf)}, add = TRUE)
  expect_no_error(plot_theta(filtered_Z, th, sim_conf = FALSE))
  expect_no_error(plot_mean_cluster_size(filtered_Z, th, sim_conf = FALSE))
  th0 <- suppressMessages(fit_theta(filtered_Z[1], modifyList(opts, list(N.theta = 0))))
  expect_s3_class(th0$theta[[1]], 'data.frame')
  expect_error(plot_theta(filtered_Z[1], th0), 'N.theta > 0')
})

test_that('fit_GL fits the radial likelihood to cluster radial maxima', {
  set.seed(21)
  n <- 2000
  e <- mvtnorm::rmvnorm(n, c(0, 0), matrix(c(1, .7, .7, 1), 2))
  Y <- e
  for (t in 2:n) Y[t, ] <- 0.8 * Y[t - 1, ] + 0.6 * e[t, ]
  X <- rmutil::qlaplace(pnorm(Y))
  op <- set.options(X = X, excess.dist.fam = 'E', W.model = 'M1', W.data = 'ExcOnly',
                    use.mean.Qq = FALSE, q = 0.9, N.Qq = 2, N.GW = 5,
                    mesh.knots.2d = seq(-pi, pi, by = pi / 32), seed = 0L)
  cfg <- set.configs(save = FALSE, progress = FALSE)
  fQ <- suppressMessages(fit_Qq(X, op, cfg))
  rm <- get_radial_maxima(fQ, m = 5)

  fQ2 <- fQ; fQ2$options$W.model <- 'M2'
  expect_error(fit_GL(fQ2, cfg, radial.idx = rm), 'M2')
  expect_error(fit_GL(fQ, cfg, radial.idx = rm[1]), 'one element per posterior sample')
  bad <- rm; bad[[1]][1] <- setdiff(seq_len(n), get_exceedances(fQ)[[1]])[1]
  expect_error(suppressMessages(fit_GL(fQ, cfg, radial.idx = bad)), 'indices of exceedances')
  fQ3 <- list(X = cbind(X, 1), options = op)
  expect_error(fit_GL(fQ3, cfg, radial.idx = rm), 'd=3')

  fm <- suppressMessages(fit_GL(fQ, cfg, return_fitted_obj = TRUE))
  fm_all <- suppressMessages(fit_GL(fQ, cfg, return_fitted_obj = TRUE, radial.idx = get_exceedances(fQ)))
  for (i in 1:2)
    expect_equal(fm_all$GL.fit.obj[[i]]$summary.fixed$mean, fm$GL.fit.obj[[i]]$summary.fixed$mean, tolerance = 0.05)
  expect_null(fm$radial.idx)

  fm_rm <- suppressMessages(fit_GL(fQ, cfg, radial.idx = rm))
  expect_identical(fm_rm$radial.idx, rm)
  expect_gt(mean(unlist(fm_rm$G)), mean(unlist(fm$G)))

  fQ4 <- fQ; fQ4$options$W.model <- 'M3'
  expect_warning(suppressMessages(fit_GL(fQ4, cfg, radial.idx = rm)), 'M3')
})
