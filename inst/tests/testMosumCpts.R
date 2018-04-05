N_TEST <- 5 # replications of tests with random input data

test_that("Set of detected change points increases with decreasing epsilon", {
  for (i in 1:N_TEST) {
    alpha <- runif(1, 0, 1)
    ts <- list(piecewiseStationary_timeSeries(model="blocks"),
               piecewiseStationary_timeSeries(model="fms"),
               piecewiseStationary_timeSeries(model="mix"),
               piecewiseStationary_timeSeries(model="stairs10"),
               piecewiseStationary_timeSeries(model="teeth10"))
    for (x in ts) {
      cpts.prev <- numeric(0)
      cpts.prev.b <- numeric(0) # Include boundary extension
      G <- floor(runif(1, 5, 15))
      for (epsilon in (11:0)/10) { #TODO include 0
        cpts.curr <- mosum.cpts(x, G,
                                alpha=alpha, criterion="epsilon", epsilon=epsilon,
                                boundary.extension=F)$cpts
        cpts.curr.b <- mosum.cpts(x, G, boundary.extension=T,
                                  alpha=alpha, criterion="epsilon", epsilon=epsilon)$cpts
        expect_true(all(cpts.prev %in% cpts.curr))
        expect_true(all(cpts.prev.b %in% cpts.curr.b))
        cpts.prev <- cpts.curr
        cpts.prev.b <- cpts.curr.b
      }
    }
  }
})
test_that("Set of detected change points increases with decreasing eta", {
  for (i in 1:N_TEST) {
    alpha <- runif(1, 0, 1)
    ts <- list(piecewiseStationary_timeSeries(model="blocks"),
               piecewiseStationary_timeSeries(model="fms"),
               piecewiseStationary_timeSeries(model="mix"),
               piecewiseStationary_timeSeries(model="stairs10"),
               piecewiseStationary_timeSeries(model="teeth10"))
    for (x in ts) {
      cpts.prev <- numeric(0)
      cpts.prev.b <- numeric(0) # Include boundary extension
      G <- floor(runif(1, 5, 15))
      for (eta in (20:0)/10) {
        cpts.curr <- mosum.cpts(x, G,
                                alpha=alpha, criterion="eta", eta=eta,
                                boundary.extension=F)$cpts
        cpts.curr.b <- mosum.cpts(x, G, boundary.extension=T,
                                alpha=alpha, criterion="eta", eta=eta)$cpts
        expect_true(all(cpts.prev %in% cpts.curr))
        expect_true(all(cpts.prev.b %in% cpts.curr.b))
        cpts.prev <- cpts.curr
        cpts.prev.b <- cpts.curr.b
      }
    }
  }
})