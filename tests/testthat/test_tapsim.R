library(rapt)

#### readResult ####
test_that("TAPSim results read correctly", {
  res.3225 <- readResult('results_3225.00000001')
  expect_named(res.3225, c("index","id","number","voltage",
                           "startX","startY","startZ",
                           "stopX","stopY","stopZ","tof","probability",
                           "potentialBefore",
                           "fieldBeforeX","fieldBeforeY","fieldBeforeZ",
                           "potentialAfter",
                           "fieldAfterX","fieldAfterY","fieldAfterZ",
                           "normalX","normalY","normalZ",
                           "apexX","apexY","apexZ"))
  expect_type(res.3225, "list")
  expect_s3_class(res.3225, "data.frame")

  res.557 <- readResult('results_557.00000001')
  expect_named(res.557, c("index","id","number","voltage",
                          "startX","startY","startZ",
                          "stopX","stopY","stopZ","tof","probability",
                          "normalX","normalY","normalZ",
                          "apexX","apexY","apexZ"))
  expect_type(res.557, "list")
  expect_s3_class(res.557, "data.frame")
})

#### resultToPOS ####
test_that("TAPSim results start positions convert to a POS", {
  res.3225 <- readResult('results_3225.00000001')
  pos.3225 <- resultToPOS(res.3225)

  expect_named(pos.3225, c("x","y","z","mass"))
  expect_type(pos.3225, "list")
  expect_s3_class(pos.3225, "data.frame")
})

#### resultToDet ####
test_that("TAPSim results end positions convert to detector ppp", {
  res.3225 <- readResult('results_3225.00000001')
  det.3225 <- resultToDet(res.3225)

  expect_type(det.3225, "list")
  expect_s3_class(det.3225, "ppp")
})
