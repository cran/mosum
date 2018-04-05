N_TEST <- 20

test_that("Custom threshold is consistent to critical value", {
  for (i in 1:N_TEST) {
    ts <- list(piecewiseStationary_timeSeries(model="blocks"),
               piecewiseStationary_timeSeries(model="fms"),
               piecewiseStationary_timeSeries(model="mix"),
               piecewiseStationary_timeSeries(model="stairs10"),
               piecewiseStationary_timeSeries(model="teeth10"))
    for (x in ts) {
      n <- length(x)
      alpha <- runif(1, 0, 1)
      G_left <- max(floor(runif(1, 20, 40)), ceiling(runif(1, 0.05*n, 0.2*n)))
      G_right <- max(floor(runif(1, 20, 40)), ceiling(runif(1, 0.05*n, 0.2*n)))
      
      # Just re-define the critical value, but in 'custom' branch of code
      threshold.custom <- mosum.criticalValue(n, G_left, G_right, alpha)
      threshold.function1 <- function(G) { mosum.criticalValue(n,G,G,alpha) }
      threshold.function2 <- function(G_l, G_r) { mosum.criticalValue(n,G_l,G_r,alpha) }
      
      mcpts1 <- mosum.cpts(x, G_left, G.right=G_right, alpha=alpha)
      mcpts2 <- multiscale.bottomUp.cpts(x, G=c(G_left,G_right), alpha=alpha)
      mcpts3 <- multiscale.cpts(x, G=c(G_left,G_right), alpha=alpha)
      
      mcpts4 <- mosum.cpts(x, G_left, G.right=G_right, 
                           threshold="custom", 
                           threshold.custom=threshold.custom)
      mcpts5 <- multiscale.bottomUp.cpts(x, G=c(G_left,G_right), 
                                         threshold="custom", 
                                         threshold.function=threshold.function1)
      mcpts6 <- multiscale.cpts(x, G=c(G_left,G_right), 
                                threshold="custom",
                                threshold.function=threshold.function2)

      expect_equal(mcpts1$cpts, mcpts4$cpts)
      expect_equal(mcpts2$cpts, mcpts5$cpts)
      expect_equal(mcpts3$cpts, mcpts6$cpts)
    }
  }
})