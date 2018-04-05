N_TEST <- 5 # replications of tests with random input data

test_that("Single-bandwidth Multiscale CPTS bottom up equals MOSUM", {
  for (i in 1:N_TEST) {
    alpha <- runif(1, 0, 1)
    eta <- runif(1, 0, 1)
    epsilon <- runif(1, 0, 1)
    ts <- list(piecewiseStationary_timeSeries(model="blocks"),
               piecewiseStationary_timeSeries(model="fms"),
               piecewiseStationary_timeSeries(model="mix"),
               piecewiseStationary_timeSeries(model="stairs10"),
               piecewiseStationary_timeSeries(model="teeth10"))
    for (x in ts) {
      n <- length(x)
      G <- max(floor(runif(1, 20, 40)), ceiling(runif(1, 0.05*n, 0.2*n)))
      cpts.mosum <- mosum.cpts(x, G, alpha=alpha, criterion="eta",
                               eta=eta)$cpts
      cpts.moba <- multiscale.bottomUp.cpts(x, G=G, eta=eta,
                                            alpha=alpha)$cpts
      names(cpts.moba) <- NULL
      cpts.moba <- as.numeric(cpts.moba)
      expect_equal(cpts.moba, cpts.mosum)
    }
  }
})