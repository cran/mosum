N_TEST <- 5 # replications of tests with random input data
N_BOOTSTRAP <- 1000
NUMERICAL_TOL <- 1e-12

test_that("Width of bootstrap intervals is consistent", {
  for (i in 1:N_TEST) {
    alpha <- runif(1, 0, 1)
    # eta <- runif(1, 0, 1)
    # epsilon <- runif(1, 0, 1)
    ts <- list(piecewiseStationary_timeSeries(model="blocks"),
               piecewiseStationary_timeSeries(model="fms"),
               piecewiseStationary_timeSeries(model="mix"),
               piecewiseStationary_timeSeries(model="stairs10"),
               piecewiseStationary_timeSeries(model="teeth10"))
    for (x in ts) {
      n <- length(x)
      G_left <- max(floor(runif(1, 20, 40)), ceiling(runif(1, 0.05*n, 0.1*n)))
      G_right <- max(floor(runif(1, 20, 40)), ceiling(runif(1, 0.05*n, 0.1*n)))
      mcpts1 <- mosum.cpts(x, G_left, G.right=G_right, alpha=alpha)
      mcpts2 <- multiscale.bottomUp.cpts(x, G=c(G_left,G_right), alpha=alpha)
      mcpts3 <- multiscale.cpts(x, G=c(G_left,G_right), alpha=alpha)

      s1 <- sample(1:10000,1)
      s2 <- sample(1:10000,1)
      s3 <- sample(1:10000,1)

      alpha1 <- runif(1, 0, .5)
      set.seed(s1); b1 <- cpts.bootstrap(mcpts1, N_BOOTSTRAP, alpha=alpha1)
      set.seed(s2); b2 <- cpts.bootstrap(mcpts2, N_BOOTSTRAP, alpha=alpha1)
      set.seed(s3); b3 <- cpts.bootstrap(mcpts3, N_BOOTSTRAP, alpha=alpha1)

      # Pointwise intervals contain estimate
      expect_true(all(b1$CI[,2] <= b1$CI[,1] & b1$CI[,3] >= b1$CI[,1]))
      expect_true(all(b2$CI[,2] <= b2$CI[,1] & b2$CI[,3] >= b2$CI[,1]))
      expect_true(all(b3$CI[,2] <= b3$CI[,1] & b3$CI[,3] >= b3$CI[,1]))

      # Uniform intervals contain estimate
      expect_true(all(b1$CI[,4]-NUMERICAL_TOL <= b1$CI[,1] & b1$CI[,5]+NUMERICAL_TOL >= b1$CI[,1]))
      expect_true(all(b2$CI[,4]-NUMERICAL_TOL <= b2$CI[,1] & b2$CI[,5]+NUMERICAL_TOL >= b2$CI[,1]))
      expect_true(all(b3$CI[,4]-NUMERICAL_TOL <= b3$CI[,1] & b3$CI[,5]+NUMERICAL_TOL >= b3$CI[,1]))

      alpha2 <- runif(1, 0, alpha1)
      set.seed(s1); b4 <- cpts.bootstrap(mcpts1, N_BOOTSTRAP, alpha=alpha2)
      set.seed(s2); b5 <- cpts.bootstrap(mcpts2, N_BOOTSTRAP, alpha=alpha2)
      set.seed(s3); b6 <- cpts.bootstrap(mcpts3, N_BOOTSTRAP, alpha=alpha2)

      # Pointwise Intervals grow (with smaller alpha)
      expect_true(all(b1$CI[,2] >= b4$CI[,2] & b1$CI[,3] <= b4$CI[,3]))
      expect_true(all(b2$CI[,2] >= b5$CI[,2] & b2$CI[,3] <= b5$CI[,3]))
      expect_true(all(b3$CI[,2] >= b6$CI[,2] & b3$CI[,3] <= b6$CI[,3]))

      # Uniform Intervals grow (with smaller alpha)
      expect_true(all(b1$CI[,4]+NUMERICAL_TOL >= b4$CI[,4] & b1$CI[,5]-NUMERICAL_TOL <= b4$CI[,5]))
      expect_true(all(b2$CI[,4]+NUMERICAL_TOL >= b5$CI[,4] & b2$CI[,5]-NUMERICAL_TOL <= b5$CI[,5]))
      expect_true(all(b3$CI[,4]+NUMERICAL_TOL >= b6$CI[,4] & b3$CI[,5]-NUMERICAL_TOL <= b6$CI[,5]))
    }
  }
})
