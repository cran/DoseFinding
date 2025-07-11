test_that("emax", {
  emx1 <- guesst(d=0.3, p=0.8, model="emax")
  expect_equal(unname(emax(0.3,0,1,emx1)),
               0.8,
               tolerance = 0.001)
  
})

test_that("emax local", {
  emx2 <- guesst(d=0.3, p=0.8, model="emax", local = TRUE, Maxd = 1)
  expect_equal(unname(emax(0.3,0,1,emx2)/emax(1,0,1,emx2)),
               0.8,
               tolerance = 0.001)
})

test_that("betaMod", {
  bta <- guesst(d=0.4, p=0.8, model="betaMod", dMax=0.8, scal=1.2, Maxd=1)
  expect_equal(betaMod(c(0.4,0.8), 0, 1, bta[1], bta[2], scal=1.2),
               c(0.8, 1.0),
               tolerance = 0.001)
})

test_that("exponential", {
  expo <- guesst(d = 0.8, p = 0.5, "exponential", Maxd=1)
  expect_equal(unname(exponential(0.8,0,1,expo)/exponential(1,0,1,expo)),
               0.5,
               tolerance = 0.001)
})

test_that("quadratic", {
  quad <- guesst(d = 0.7, p = 1, "quadratic")
  mm <- Mods(quadratic=quad, doses=c(0,0.7,1))
  expect_equal(getResp(mm)[2],
               1,
               tolerance = 0.001)
})

test_that("logistic", {
  lgc1 <- guesst(d = c(0.2, 0.6), p = c(0.2, 0.95), "logistic")
  expect_equal(logistic(c(0.2,0.6), 0, 1, lgc1[1], lgc1[2]),
               c(0.2, 0.95),
               tolerance = 0.001)
})

test_that("logistic local", {
  lgc2 <- guesst(d = c(0.2, 0.6), p = c(0.2, 0.95), "logistic", 
                 local = TRUE, Maxd = 1)
  r0 <- logistic(0, 0, 1, lgc2[1], lgc2[2])
  r1 <- logistic(1, 0, 1, lgc2[1], lgc2[2])
  expect_equal((logistic(c(0.2,0.6), 0, 1, lgc2[1], lgc2[2])-r0)/(r1-r0),
               c(0.2, 0.95),
               tolerance = 0.001)
})

test_that("sigEmax", {
  sgE1 <- guesst(d = c(0.2, 0.6), p = c(0.2, 0.95), "sigEmax")
  expect_equal(sigEmax(c(0.2,0.6), 0, 1, sgE1[1], sgE1[2]),
               c(0.2, 0.95),
               tolerance = 0.001)
})

test_that("sigEmax local", {
  sgE2 <- guesst(d = c(0.2, 0.6), p = c(0.2, 0.95), "sigEmax",
                 local = TRUE, Maxd = 1)
  r1 <- sigEmax(1, 0, 1, sgE2[1], sgE2[2])
  expect_equal(sigEmax(c(0.2,0.6), 0, 1, sgE2[1], sgE2[2])/r1,
               c(0.2,0.95),
               tolerance = 0.001)
})


## test error conditions


test_that("Error conditions for guesst function", {
  
  # Test for invalid percentage values (negative or greater than 1)
  expect_error(guesst(d = 0.5, p = -0.2, model = "emax"), "must have 0 < p <= 1")
  expect_error(guesst(d = 0.5, p = 1.2, model = "emax"), "must have 0 < p <= 1")
  
  # Test for logistic model needing at least two pairs
  expect_error(guesst(d = 0.2, p = 0.5, model = "logistic"), "logistic model needs at least two pairs")
  
  # Test for local version of emax with p <= d/Maxd
  expect_error(guesst(d = 0.3, p = 0.2, model = "emax", local = TRUE, Maxd = 1), "must have p > d/Maxd, for local version")
  
  # Test for exponential model needing p < d/Maxd
  expect_error(guesst(d = 0.8, p = 0.9, model = "exponential", Maxd = 1), "must have p < d/Maxd")
  
  # Test for betaMod model needing scal > dMax
  expect_error(guesst(d = 0.4, p = 0.8, model = "betaMod", dMax = 0.8, scal = 0.8, Maxd = 1), "scal needs to be larger than dMax to calculate guesstimate")
  
  # Test for betaMod model needing dMax <= Maxd
  expect_error(guesst(d = 0.4, p = 0.8, model = "betaMod", dMax = 1.2, scal = 1.5, Maxd = 1), "dose with maximum effect \\(dMax\\) needs to be smaller than maximum dose \\(Maxd\\)")
  
  # Test for sigmoid Emax model needing at least two pairs
  expect_error(guesst(d = 0.2, p = 0.5, model = "sigEmax"), "sigmoid Emax model needs at least two pairs")
  
})
