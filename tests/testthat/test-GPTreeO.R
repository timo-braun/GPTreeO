## load the data set beforehand
test_that("training works and mean of node 0", {
  set.seed(42)
  X <- as.matrix(sample(seq(0, 10, length.out = 51)))
  y <- sin(2 * pi * X / 10) + 0.2 * sin(2 * pi * X / 2.5)
  y_variance <- rep(0.1**2, 51)

  gptree = GPTree$new(Nbar = 50, gradual_split = FALSE, theta = 0.05)

  for (i in 1:nrow(X)) {
    gptree$update(X[i,], y[i], y_variance[i])
  }

  ## Select first 50 data points
  X_50 = as.matrix(X[1:50,])

  ## Test that mean is equal
  expect_equal(gptree$nodes$`0`$position_split, median(X_50))
})

test_that("GP gets removed after splitting", {
  set.seed(42)
  X <- as.matrix(sample(seq(0, 10, length.out = 101)))
  y <- sin(2 * pi * X / 10) + 0.2 * sin(2 * pi * X / 2.5)
  y_variance <- rep(0.1**2, 101)

  gptree = GPTree$new(Nbar = 50, retrain_buffer_length = 50)
  
  for (i in 1:51) {
    gptree$update(X[i,], y[i], y_variance[i])
  }

  ## Test that the wrapped_gp was removed
  expect_equal(gptree$nodes$`0`$wrapped_gp, NULL)

})

test_that("node is no longer in leaf_keys after splitting", {
  set.seed(42)
  X <- as.matrix(sample(seq(0, 10, length.out = 101)))
  y <- sin(2 * pi * X / 10) + 0.2 * sin(2 * pi * X / 2.5)
  y_variance <- rep(0.1**2, 101)

  gptree = GPTree$new(Nbar = 50, retrain_buffer_length = 50)
  
  for (i in 1:51) {
    gptree$update(X[i,], y[i], y_variance[i])
  }

  ## Test that the node 0 is no longer in the leaf_keys
  expect_equal("0" %in% gptree$leaf_keys, FALSE)

})


